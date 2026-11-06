use std::{
    collections::BTreeMap,
    f64::consts::FRAC_PI_2,
};

type P3 = [f64; 3];

/// Roof faces of a building's points on a grid of `res` cells: the surface height of each cell,
/// each cell's face (0 for none) and the faces' planes, in grid coordinates.
pub struct Faces {
    pub x0: f64,
    pub y0: f64,
    pub w: usize,
    pub ht: usize,
    pub z: Vec<f64>,
    pub label: Vec<usize>,
    pub(crate) planes: Vec<Plane>,
}

pub fn faces(points: &[P3], res: f64) -> Faces {
    let mut q = [f64::MAX, f64::MAX, f64::MIN, f64::MIN];
    for p in points {
        q = [q[0].min(p[0]), q[1].min(p[1]), q[2].max(p[0]), q[3].max(p[1])];
    }
    let (x0, y0) = (q[0] - res, q[1] - res);
    let w = ((q[2] - x0) / res) as usize + 2;
    let ht = ((q[3] - y0) / res) as usize + 2;
    let cell = |p: &P3| ((p[1] - y0) / res) as usize * w + ((p[0] - x0) / res) as usize;

    // The median of a cell's points: flight lines combined leave layers a few cm apart and stray
    // points, and the highest point would make the surface rough.
    let mut zs: Vec<Vec<f64>> = vec![Vec::new(); w * ht];
    for p in points {
        zs[cell(p)].push(p[2]);
    }
    let mut z: Vec<f64> = zs
        .iter_mut()
        .map(|v| {
            if v.is_empty() {
                return f64::NAN;
            }
            v.sort_by(f64::total_cmp);
            v[v.len() / 2]
        })
        .collect();
    let nb = |k: usize| {
        let (i, j) = ((k % w) as isize, (k / w) as isize);
        [(-1, 0), (1, 0), (0, -1), (0, 1)].into_iter().filter_map(move |(di, dj)| {
            let (a, b) = (i + di, j + dj);
            (a >= 0 && b >= 0 && (a as usize) < w && (b as usize) < ht).then(|| b as usize * w + a as usize)
        })
    };

    // Empty cells between points get their neighbours' mean, so the surface has no pits.
    let mask: Vec<bool> = {
        let occ: Vec<bool> = z.iter().map(|v| !v.is_nan()).collect();
        let grown: Vec<bool> = (0..z.len()).map(|k| occ[k] || nb(k).any(|n| occ[n])).collect();
        (0..z.len()).map(|k| grown[k] && (nb(k).all(|n| grown[n]) || occ[k])).collect()
    };
    for _ in 0..4 {
        let prev = z.clone();
        for k in 0..z.len() {
            if mask[k] && prev[k].is_nan() {
                let v: Vec<f64> = nb(k).map(|n| prev[n]).filter(|v| !v.is_nan()).collect();
                if !v.is_empty() {
                    z[k] = v.iter().sum::<f64>() / v.len() as f64;
                }
            }
        }
    }
    let inside: Vec<bool> = (0..z.len()).map(|k| mask[k] && !z[k].is_nan()).collect();

    // Local plane of each cell from its full 3×3 neighbourhood: slopes a, b (per m) and fit RMS.
    let at = |i: isize, j: isize| {
        (i >= 0 && j >= 0 && (i as usize) < w && (j as usize) < ht).then(|| j as usize * w + i as usize)
    };
    let local: Vec<Option<(f64, f64, f64)>> = (0..z.len())
        .map(|k| {
            if !inside[k] {
                return None;
            }
            let (i, j) = ((k % w) as isize, (k / w) as isize);
            let mut v = [[0.0; 3]; 3];
            for dj in -1..=1 {
                for di in -1..=1 {
                    let n = at(i + di, j + dj).filter(|&n| inside[n])?;
                    v[(dj + 1) as usize][(di + 1) as usize] = z[n];
                }
            }
            let (mut sa, mut sb, mut sz) = (0.0, 0.0, 0.0);
            for dj in 0..3 {
                for di in 0..3 {
                    sa += (di as f64 - 1.0) * v[dj][di];
                    sb += (dj as f64 - 1.0) * v[dj][di];
                    sz += v[dj][di];
                }
            }
            let (a, b, c) = (sa / 6.0, sb / 6.0, sz / 9.0);
            let mut rss = 0.0;
            for dj in 0..3 {
                for di in 0..3 {
                    rss += (v[dj][di] - (c + a * (di as f64 - 1.0) + b * (dj as f64 - 1.0))).powi(2);
                }
            }
            Some((a / res, b / res, (rss / 9.0).sqrt()))
        })
        .collect();
    let pos = |k: usize| [((k % w) as f64 + 0.5) * res, ((k / w) as f64 + 0.5) * res];

    // Roof planes grown from the flattest cells while the next cell stays on the plane.
    let mut label = vec![0usize; z.len()];
    let mut planes: Vec<Plane> = vec![Plane::default()];
    let mut order: Vec<usize> = (0..z.len()).filter(|&k| local[k].is_some_and(|l| l.2 < SEED_RMS)).collect();
    order.sort_by(|&a, &b| local[a].unwrap().2.total_cmp(&local[b].unwrap().2));
    for &s in &order {
        if label[s] != 0 {
            continue;
        }
        let r = planes.len();
        let (a, b, _) = local[s].unwrap();
        let ps = pos(s);
        let mut plane = Plane { a, b, c: z[s] - a * ps[0] - b * ps[1], ..Plane::default() };
        let mut fit = Plane::default();
        // The seed's whole 3×3 neighbourhood is on the plane, so the face starts from all of it.
        let (si, sj) = ((s % w) as isize, (s / w) as isize);
        let mut queue = std::collections::VecDeque::new();
        for k in (-1..=1).flat_map(|dj| (-1..=1).filter_map(move |di| at(si + di, sj + dj))) {
            if label[k] == 0 {
                label[k] = r;
                fit.add(pos(k), z[k]);
                queue.push_back(k);
            }
        }
        let mut next_refit = 9;
        while let Some(k) = queue.pop_front() {
            for n in nb(k) {
                if !inside[n] || label[n] != 0 {
                    continue;
                }
                let pn = pos(n);
                if (z[n] - plane.z(pn)).abs() > PLANE_TOL {
                    continue;
                }
                // A face grows only where it already surrounds the cell, never through a one-cell neck
                // such as the corner where two neighbours' coplanar roof faces touch.
                let (ni, nj) = ((n % w) as isize, (n / w) as isize);
                let support = (-1..=1)
                    .flat_map(|dj| (-1..=1).filter_map(move |di| at(ni + di, nj + dj)))
                    .filter(|&m| label[m] == r)
                    .count();
                if support < MIN_SUPPORT {
                    continue;
                }
                if let Some((an, bn, rms)) = local[n] {
                    if rms < SEED_RMS && (an - plane.a).abs() + (bn - plane.b).abs() > SLOPE_TOL {
                        continue;
                    }
                }
                label[n] = r;
                fit.add(pn, z[n]);
                queue.push_back(n);
                if fit.n >= next_refit {
                    if let Some(p) = fit.solve() {
                        plane = Plane { a: p.0, b: p.1, c: p.2, ..Plane::default() };
                    }
                    next_refit *= 2;
                }
            }
        }
        if let Some(p) = fit.solve() {
            plane = Plane { a: p.0, b: p.1, c: p.2, ..Plane::default() };
        }
        plane.n = fit.n;
        planes.push(plane);
    }
    // Planes too small to be a roof face are dissolved; their cells, and ridge and edge cells no
    // plane took, go to the neighbouring plane they fit best.
    for l in label.iter_mut() {
        if *l != 0 && planes[*l].n < MIN_PLANE_CELLS {
            *l = 0;
        }
    }
    loop {
        let mut changed = false;
        for k in 0..z.len() {
            if !inside[k] || label[k] != 0 {
                continue;
            }
            let pk = pos(k);
            if let Some(r) = nb(k).map(|n| label[n]).filter(|&r| r != 0).min_by(|&r, &s| {
                (z[k] - planes[r].z(pk)).abs().total_cmp(&(z[k] - planes[s].z(pk)).abs())
            }) {
                label[k] = r;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }

    // Coplanar faces of two neighbouring roofs can still join through a neck a couple of cells wide.
    // Shrunk by one cell such a face falls apart; each piece then regrows into its share of it.
    let n_faces = planes.len();
    // Each face's cells, so a face's pass costs its own size, not the whole blob's.
    let mut cells_of: Vec<Vec<usize>> = vec![Vec::new(); n_faces];
    for (k, &l) in label.iter().enumerate() {
        if l != 0 {
            cells_of[l].push(k);
        }
    }
    let mut piece = vec![0usize; z.len()];
    let mut owner = vec![0usize; z.len()];
    for r in 1..n_faces {
        let cells = &cells_of[r];
        let core = |k: usize| {
            let (i, j) = ((k % w) as isize, (k / w) as isize);
            (-1..=1).all(|dj| (-1..=1).all(|di| at(i + di, j + dj).is_some_and(|m| label[m] == r)))
        };
        let mut pieces = Vec::new();
        for &s in cells {
            if piece[s] != 0 || !core(s) {
                continue;
            }
            let id = pieces.len() + 1;
            piece[s] = id;
            let mut p = vec![s];
            let mut n = 0;
            while n < p.len() {
                for m in nb(p[n]) {
                    if label[m] == r && piece[m] == 0 && core(m) {
                        piece[m] = id;
                        p.push(m);
                    }
                }
                n += 1;
            }
            pieces.push(p);
        }
        for &k in pieces.iter().flatten() {
            piece[k] = 0;
        }
        let big: Vec<usize> = (0..pieces.len()).filter(|&p| pieces[p].len() >= MIN_PLANE_CELLS).collect();
        if big.len() < 2 {
            continue;
        }
        // The largest piece keeps the face; the others become new faces, all regrown within it.
        let mut queue = std::collections::VecDeque::new();
        for (n, &p) in big.iter().enumerate() {
            let l = if n == 0 { r } else { planes.len() + n - 1 };
            for &k in &pieces[p] {
                owner[k] = l;
                queue.push_back(k);
            }
        }
        while let Some(k) = queue.pop_front() {
            for m in nb(k) {
                if label[m] == r && owner[m] == 0 {
                    owner[m] = owner[k];
                    queue.push_back(m);
                }
            }
        }
        let first_new = planes.len();
        planes.extend((1..big.len()).map(|_| Plane::default()));
        let mut fits: Vec<Plane> = (0..big.len()).map(|_| Plane::default()).collect();
        for &k in cells {
            if owner[k] != 0 {
                label[k] = owner[k];
                let f = if owner[k] == r { 0 } else { owner[k] - first_new + 1 };
                fits[f].add(pos(k), z[k]);
                owner[k] = 0;
            }
        }
        for (f, fit) in fits.iter().enumerate() {
            let l = if f == 0 { r } else { first_new + f - 1 };
            if let Some((a, b, c)) = fit.solve() {
                planes[l] = Plane { a, b, c, n: fit.n, ..Plane::default() };
            }
        }
    }
    Faces { x0, y0, w, ht, z, label, planes }
}

/// Splits the points of one building blob into houses by their roof faces: faces meeting at a ridge
/// or hip are one roof; a valley, step or offset ridge between faces separates them. Parts smaller
/// than `min_part` (m²) join the neighbour they share most boundary with.
pub fn parts(points: &[P3], thetas: &[f64], res: f64, min_part: f64, debug: bool) -> (Vec<Vec<P3>>, Vec<Wall>) {
    let Faces { x0, y0, w, ht, z, mut label, planes } = faces(points, res);
    let cell = |p: &P3| ((p[1] - y0) / res) as usize * w + ((p[0] - x0) / res) as usize;
    let nb = |k: usize| neighbours(k, w, ht);
    let at = |i: isize, j: isize| {
        (i >= 0 && j >= 0 && (i as usize) < w && (j as usize) < ht).then(|| j as usize * w + i as usize)
    };
    let pos = |k: usize| [((k % w) as f64 + 0.5) * res, ((k / w) as f64 + 0.5) * res];

    // Faces of one roof meet at ridges and hips: continuous, and each falling away from the edge.
    // Valleys, steps and offset ridges separate roofs, or wings that are easier joined than split.
    // Per pair: boundary cells that look like one roof, all boundary cells, of those coplanar, and of
    // those continuous; and the boundary itself.
    let mut votes: BTreeMap<(usize, usize), (usize, usize, usize, usize)> = BTreeMap::new();
    let mut edge_pts: BTreeMap<(usize, usize), Vec<[f64; 2]>> = BTreeMap::new();
    for k in 0..z.len() {
        if label[k] == 0 {
            continue;
        }
        for n in nb(k) {
            let (r, s) = (label[k], label[n]);
            if s == 0 || r == s {
                continue;
            }
            let (pk, pn) = (pos(k), pos(n));
            let m = [(pk[0] + pn[0]) / 2.0, (pk[1] + pn[1]) / 2.0];
            let (pr, ps) = (&planes[r], &planes[s]);
            let Meeting { continuous, ridge, flat } = meet(pr, ps, pk, pn, res);
            // Pieces of one face that grew apart.
            let coplanar = (pr.a - ps.a).abs() + (pr.b - ps.b).abs() < SLOPE_TOL;
            edge_pts.entry((r.min(s), r.max(s))).or_default().push(m);
            let e = votes.entry((r.min(s), r.max(s))).or_default();
            e.1 += 1;
            e.3 += continuous as usize;
            if continuous && (ridge || flat) {
                e.0 += 1;
            } else if continuous && coplanar {
                e.2 += 1;
            }
        }
    }
    if debug {
        for (r, p) in planes.iter().enumerate().skip(1) {
            eprintln!("face {r}: slope ({:.2}, {:.2}) cells {}", p.a, p.b, p.n);
        }
        for ((r, s), (same, total, coplanar, _)) in &votes {
            eprintln!("faces {r}-{s}: {same} of {total} boundary cells look like one roof, {coplanar} coplanar");
        }
    }
    // Outline length of each face, in cell edges.
    let mut perim: BTreeMap<usize, usize> = BTreeMap::new();
    for k in 0..z.len() {
        if label[k] != 0 {
            let (i, j) = ((k % w) as isize, (k / w) as isize);
            let edges = [(-1, 0), (1, 0), (0, -1), (0, 1)]
                .iter()
                .filter(|(di, dj)| at(i + di, j + dj).is_none_or(|m| label[m] != label[k]))
                .count();
            *perim.entry(label[k]).or_default() += edges;
        }
    }
    let mut root: Vec<usize> = (0..planes.len()).collect();
    fn find(root: &mut [usize], x: usize) -> usize {
        let mut x = x;
        while root[x] != x {
            root[x] = root[root[x]];
            x = root[x];
        }
        x
    }
    for (&(r, s), &(same, total, coplanar, cont)) in &votes {
        // Neighbours' coplanar faces touch at a corner or a neck; pieces of one face share much of
        // their outline.
        let outline = perim[&r].min(perim[&s]) as f64;
        let pieces = total >= MIN_COPLANAR_EDGE && total as f64 >= COPLANAR_SHARE * outline;
        let same = if pieces { same + coplanar } else { same };
        // A wing meets the main roof in a valley running diagonally across the walls; neighbouring
        // houses meet in one running along them.
        let wing = cont as f64 > SAME_ROOF_SHARE * total as f64
            && total >= MIN_COPLANAR_EDGE
            && diagonal(&edge_pts[&(r, s)], thetas);
        if wing || same as f64 > SAME_ROOF_SHARE * total as f64 {
            if debug {
                let c = |f: usize| {
                    let v: Vec<[f64; 2]> = (0..z.len()).filter(|&k| label[k] == f).map(pos).collect();
                    let n = v.len() as f64;
                    v.iter().fold([x0, y0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n])
                };
                let (cr, cs) = (c(r), c(s));
                eprintln!("join {r} ({:.1}, {:.1}) + {s} ({:.1}, {:.1}): {same}/{total}", cr[0], cr[1], cs[0], cs[1]);
            }
            let (fr, fs) = (find(&mut root, r), find(&mut root, s));
            root[fr] = fs;
        }
    }
    for l in label.iter_mut() {
        if *l != 0 {
            *l = find(&mut root, *l);
        }
    }

    // Slivers, and strips too narrow for a house (side hips, eaves, porch roofs), join the neighbour
    // they share most boundary with. A small roof well below that neighbour's is a shed or garage
    // against the house and stays a building of its own.
    let cell_area = res * res;
    let mut kept: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    loop {
        let mut cells: BTreeMap<usize, Vec<[f64; 2]>> = BTreeMap::new();
        for (k, &l) in label.iter().enumerate() {
            if l != 0 {
                cells.entry(l).or_default().push(pos(k));
            }
        }
        if cells.len() <= 1 {
            break;
        }
        let area: BTreeMap<usize, usize> = cells.iter().map(|(&l, c)| (l, c.len())).collect();
        let small_or_thin = |l: usize| (area[&l] as f64) * cell_area < min_part || width(&cells[&l]) < MIN_WIDTH;
        let Some((&small, _)) =
            area.iter().filter(|&(&l, _)| !kept.contains(&l) && small_or_thin(l)).min_by_key(|&(_, &a)| a)
        else {
            break;
        };
        // Per neighbour: shared edge cells, and the summed height differences across the edge.
        let mut shared: BTreeMap<usize, (usize, f64)> = BTreeMap::new();
        for k in 0..z.len() {
            if label[k] == small {
                for n in nb(k) {
                    if label[n] != 0 && label[n] != small {
                        let e = shared.entry(label[n]).or_default();
                        e.0 += 1;
                        e.1 += z[n] - z[k];
                    }
                }
            }
        }
        let Some((&into, &(edge, rise))) = shared.iter().max_by_key(|&(_, &(c, _))| c) else { break };
        if rise / edge as f64 > SHED_STEP && area[&small] as f64 * cell_area >= MIN_SHED {
            kept.insert(small);
            continue;
        }
        if debug {
            eprintln!("small or thin {small} ({} cells) into {into}", area[&small]);
        }
        for l in label.iter_mut() {
            if *l == small {
                *l = into;
            }
        }
    }

    // Boundary between each pair of neighbouring parts, as points midway between their cells.
    let boundaries = |label: &[usize]| {
        let mut pairs: BTreeMap<(usize, usize), Vec<[f64; 2]>> = BTreeMap::new();
        for k in 0..label.len() {
            if label[k] == 0 {
                continue;
            }
            for n in nb(k) {
                if label[n] != 0 && label[n] > label[k] {
                    let c = |k: usize| [x0 + ((k % w) as f64 + 0.5) * res, y0 + ((k / w) as f64 + 0.5) * res];
                    let (a, b) = (c(k), c(n));
                    pairs.entry((label[k], label[n])).or_default().push([(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0]);
                }
            }
        }
        pairs
    };

    // A party wall runs along the houses' axes; a split oblique to them is a roof split in two.
    // An annex whose roof stays at the house's eaves stays apart whatever its boundary: an
    // L-shaped one fits an oblique line.
    loop {
        let mut zs: BTreeMap<usize, Vec<f64>> = BTreeMap::new();
        for (k, &l) in label.iter().enumerate() {
            if l != 0 && !z[k].is_nan() {
                zs.entry(l).or_default().push(z[k]);
            }
        }
        let pct = |l: usize, q: f64| {
            let mut v = zs.get(&l).cloned().unwrap_or_default();
            v.sort_by(f64::total_cmp);
            if v.is_empty() { f64::NAN } else { v[((v.len() - 1) as f64 * q) as usize] }
        };
        let below = |a: usize, b: usize| pct(a, 0.95) < pct(b, 0.05) + ANNEX_RISE;
        let Some((lp, lq)) = boundaries(&label)
            .into_iter()
            .find(|&((p, q), ref b)| Wall::fit(b).is_some_and(|w| w.axis(thetas).is_none()) && !below(p, q) && !below(q, p))
            .map(|(k, _)| k)
        else {
            break;
        };
        if debug {
            eprintln!("oblique split: {lp} (roof {:.1} to {:.1}) + {lq} ({:.1} to {:.1})", pct(lp, 0.05), pct(lp, 0.95), pct(lq, 0.05), pct(lq, 0.95));
        }
        for l in label.iter_mut() {
            if *l == lq {
                *l = lp;
            }
        }
    }

    let mut point_label: Vec<usize> = points.iter().map(|p| label[cell(p)]).collect();

    // The watershed line between two houses follows the raster; their shared wall is straight.
    let mut walls = Vec::new();
    for ((lp, lq), b) in boundaries(&label) {
        let Some(mut wall) = Wall::fit(&b) else { continue };
        let Some(axis) = wall.axis(thetas) else { continue };
        wall.set_dir(axis, &b);
        // Which side is `lp`: the mean offset of its points near the wall.
        let near = |p: &P3| wall.offset(p).abs() < WALL_REACH && wall.along(p) > wall.t0 - 1.0 && wall.along(p) < wall.t1 + 1.0;
        let side_p: f64 = points.iter().zip(&point_label).filter(|(p, l)| **l == lp && near(p)).map(|(p, _)| wall.offset(p)).sum();
        for (p, l) in points.iter().zip(point_label.iter_mut()) {
            if (*l == lp || *l == lq) && near(p) {
                *l = if (wall.offset(p) > 0.0) == (side_p > 0.0) { lp } else { lq };
            }
        }
        walls.push(wall);
    }

    let mut out: BTreeMap<usize, Vec<P3>> = BTreeMap::new();
    for (p, l) in points.iter().zip(point_label) {
        if l != 0 {
            out.entry(l).or_default().push(*p);
        }
    }
    (out.into_values().collect(), walls)
}

/// How far either side of a fitted wall points are reassigned by side (m).
const WALL_REACH: f64 = 3.0;

/// A wall shared by two houses: the line through `p` along `dir`, from `t0` to `t1` along it.
pub struct Wall {
    pub p: [f64; 2],
    pub dir: [f64; 2],
    pub t0: f64,
    pub t1: f64,
}

impl Wall {
    /// Principal direction of the boundary points; `None` for boundaries too short to be a wall.
    fn fit(b: &[[f64; 2]]) -> Option<Wall> {
        let n = b.len() as f64;
        let c = b.iter().fold([0.0, 0.0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n]);
        let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
        for p in b {
            let (dx, dy) = (p[0] - c[0], p[1] - c[1]);
            sxx += dx * dx;
            sxy += dx * dy;
            syy += dy * dy;
        }
        let ang = 0.5 * (2.0 * sxy).atan2(sxx - syy);
        let mut wall = Wall { p: c, dir: [ang.cos(), ang.sin()], t0: f64::MAX, t1: f64::MIN };
        for p in b {
            let t = wall.along(p);
            wall.t0 = wall.t0.min(t);
            wall.t1 = wall.t1.max(t);
        }
        (wall.t1 - wall.t0 >= MIN_WALL).then_some(wall)
    }

    /// The blob axis (radians) within [`MAX_WALL_SKEW`] of the wall's direction.
    fn axis(&self, thetas: &[f64]) -> Option<f64> {
        let ang = self.dir[1].atan2(self.dir[0]);
        thetas
            .iter()
            .map(|&t| {
                let k = ((ang - t) / FRAC_PI_2).round();
                (t + k * FRAC_PI_2, (ang - t - k * FRAC_PI_2).abs())
            })
            .filter(|&(_, d)| d < MAX_WALL_SKEW.to_radians())
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(a, _)| a)
    }

    /// Turns the wall to `ang` about its boundary's centroid and remeasures its extent.
    fn set_dir(&mut self, ang: f64, b: &[[f64; 2]]) {
        self.dir = [ang.cos(), ang.sin()];
        self.t0 = f64::MAX;
        self.t1 = f64::MIN;
        for p in b {
            let t = self.along(p);
            self.t0 = self.t0.min(t);
            self.t1 = self.t1.max(t);
        }
    }

    pub fn offset(&self, q: &[f64]) -> f64 {
        (q[0] - self.p[0]) * -self.dir[1] + (q[1] - self.p[1]) * self.dir[0]
    }

    pub fn along(&self, q: &[f64]) -> f64 {
        (q[0] - self.p[0]) * self.dir[0] + (q[1] - self.p[1]) * self.dir[1]
    }
}

/// Largest angle between a party wall and the blob's axes (degrees).
const MAX_WALL_SKEW: f64 = 20.0;

/// Shortest boundary treated as a wall between two houses (m).
const MIN_WALL: f64 = 3.0;

/// Cells whose 3×3 neighbourhood fits a plane this well (RMS, m) can seed a roof face.
const SEED_RMS: f64 = 0.1;
/// Neighbours (of 8) a cell must already have in a face to join it.
const MIN_SUPPORT: usize = 3;

/// Largest distance of a cell from the face it joins (m).
const PLANE_TOL: f64 = 0.2;
/// Largest slope difference (sum of both components, per m) of a flat cell joining a face.
const SLOPE_TOL: f64 = 0.35;
/// Faces smaller than this many cells are dissolved into their neighbours.
const MIN_PLANE_CELLS: usize = 20;
/// Two faces are continuous where they meet within this height (m).
const JOIN_TOL: f64 = 0.3;
/// Share of the smaller face's outline that coplanar faces must share to be pieces of one face.
const COPLANAR_SHARE: f64 = 0.25;

/// Shortest boundary (cells) along which coplanar faces are taken for pieces of one face.
const MIN_COPLANAR_EDGE: usize = 6;

/// How much two faces extended over each other must rise above it to meet in a ridge or hip (m).
const CONVEX_EPS: f64 = 0.05;

/// Slopes below this (per m) are flat.
const FLAT: f64 = 0.05;
/// Share of the boundary that must look like one roof for two faces to be joined.
const SAME_ROOF_SHARE: f64 = 0.6;

/// A roof face z = a·x + b·y + c, with least-squares sums for refitting.
#[derive(Default)]
pub(crate) struct Plane {
    pub a: f64,
    pub b: f64,
    c: f64,
    n: usize,
    s: [f64; 9],
}

impl Plane {
    pub fn z(&self, p: [f64; 2]) -> f64 {
        self.a * p[0] + self.b * p[1] + self.c
    }

    fn add(&mut self, p: [f64; 2], z: f64) {
        let (x, y) = (p[0], p[1]);
        for (s, v) in self.s.iter_mut().zip([x, y, z, x * x, x * y, y * y, x * z, y * z, 1.0]) {
            *s += v;
        }
        self.n += 1;
    }

    /// Least-squares plane through the added points.
    fn solve(&self) -> Option<(f64, f64, f64)> {
        let [sx, sy, sz, sxx, sxy, syy, sxz, syz, n] = self.s;
        let m = [[sxx, sxy, sx], [sxy, syy, sy], [sx, sy, n]];
        let r = [sxz, syz, sz];
        let det = |m: [[f64; 3]; 3]| {
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        };
        let d = det(m);
        if d.abs() < 1e-9 {
            return None;
        }
        let col = |c: usize| {
            let mut mm = m;
            for i in 0..3 {
                mm[i][c] = r[i];
            }
            det(mm) / d
        };
        Some((col(0), col(1), col(2)))
    }
}

/// A small part whose roof is this much below its neighbour's along their edge is its own building (m).
const SHED_STEP: f64 = 1.0;
/// A part whose roof tops out less than this (m) above another's eaves is an annex.
const ANNEX_RISE: f64 = 0.5;
/// ... if at least this large (m²).
const MIN_SHED: f64 = 8.0;

/// Narrowest house split off a blob (m).
const MIN_WIDTH: f64 = 3.0;

/// Extent of cell centres across their principal axis.
fn width(c: &[[f64; 2]]) -> f64 {
    let n = c.len() as f64;
    let m = c.iter().fold([0.0, 0.0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n]);
    let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
    for p in c {
        let (dx, dy) = (p[0] - m[0], p[1] - m[1]);
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    let a = 0.5 * (2.0 * sxy).atan2(sxx - syy);
    let nrm = [-a.sin(), a.cos()];
    let (lo, hi) = c.iter().fold((f64::MAX, f64::MIN), |(lo, hi), p| {
        let t = (p[0] - m[0]) * nrm[0] + (p[1] - m[1]) * nrm[1];
        (lo.min(t), hi.max(t))
    });
    hi - lo
}

/// Whether a boundary runs at 45° ± [`DIAGONAL_TOL`] to the building's main axes.
fn diagonal(b: &[[f64; 2]], thetas: &[f64]) -> bool {
    let Some(t) = thetas.first() else { return false };
    let n = b.len() as f64;
    let c = b.iter().fold([0.0, 0.0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n]);
    let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
    for p in b {
        let (dx, dy) = (p[0] - c[0], p[1] - c[1]);
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    let ang = 0.5 * (2.0 * sxy).atan2(sxx - syy);
    let d = (ang - t).rem_euclid(FRAC_PI_2);
    (d - FRAC_PI_2 / 2.0).abs() < DIAGONAL_TOL.to_radians()
}

/// How far from 45° to the walls a valley may run and still join a wing (degrees).
const DIAGONAL_TOL: f64 = 15.0;


/// The 4-neighbours of cell `k` in a `w` × `ht` grid.
pub(crate) fn neighbours(k: usize, w: usize, ht: usize) -> impl Iterator<Item = usize> {
    let (i, j) = ((k % w) as isize, (k / w) as isize);
    [(-1, 0), (1, 0), (0, -1), (0, 1)].into_iter().filter_map(move |(di, dj)| {
        let (a, b) = (i + di, j + dj);
        (a >= 0 && b >= 0 && (a as usize) < w && (b as usize) < ht).then(|| b as usize * w + a as usize)
    })
}

/// How two faces meet across the boundary between neighbouring cells `pk` (face `pr`) and `pn`.
pub(crate) struct Meeting {
    pub continuous: bool,
    /// A ridge or hip, of any direction, rather than a valley.
    pub ridge: bool,
    pub flat: bool,
}

pub(crate) fn meet(pr: &Plane, ps: &Plane, pk: [f64; 2], pn: [f64; 2], res: f64) -> Meeting {
    let u = [(pn[0] - pk[0]) / res, (pn[1] - pk[1]) / res];
    let m = [(pk[0] + pn[0]) / 2.0, (pk[1] + pn[1]) / 2.0];
    let (gr, gs) = (pr.a * u[0] + pr.b * u[1], ps.a * u[0] + ps.b * u[1]);
    Meeting {
        // The line where two faces meet can lie anywhere within a cell; the steeper they differ, the
        // more their heights at the cell boundary may.
        continuous: (pr.z(m) - ps.z(m)).abs() < JOIN_TOL + (gr - gs).abs() * res / 2.0,
        // Over the other face's cell, each face extended rises above the other at a ridge or hip, and
        // drops below at a valley.
        ridge: (pr.z(pn) - ps.z(pn)) + (ps.z(pk) - pr.z(pk)) > CONVEX_EPS,
        flat: pr.a.abs() + pr.b.abs() < FLAT && ps.a.abs() + ps.b.abs() < FLAT,
    }
}
