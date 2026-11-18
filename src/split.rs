use std::{
    collections::BTreeMap,
    f64::consts::{FRAC_PI_2, PI},
};

type P3 = crate::store::Pt;

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
    /// Each cell's intensity, NaN where unknown.
    pub int: Vec<f64>,
}

pub fn faces(points: &[P3], res: f64) -> Faces {
    let min_cells = (MIN_PLANE_AREA / (res * res)).round() as usize;
    let mut q = [f64::MAX, f64::MAX, f64::MIN, f64::MIN];
    for p in points {
        q = [
            q[0].min(p[0]),
            q[1].min(p[1]),
            q[2].max(p[0]),
            q[3].max(p[1]),
        ];
    }
    let (x0, y0) = (q[0] - res, q[1] - res);
    let w = ((q[2] - x0) / res) as usize + 2;
    let ht = ((q[3] - y0) / res) as usize + 2;
    let cell = |p: &P3| ((p[1] - y0) / res) as usize * w + ((p[0] - x0) / res) as usize;

    // The median of a cell's points: flight lines combined leave layers a few cm apart and stray
    // points, and the highest point would make the surface rough.
    let mut zs: Vec<Vec<f64>> = vec![Vec::new(); w * ht];
    let mut is: Vec<Vec<f64>> = vec![Vec::new(); w * ht];
    for p in points {
        zs[cell(p)].push(p[2]);
        if !p[3].is_nan() {
            is[cell(p)].push(p[3]);
        }
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
        [(-1, 0), (1, 0), (0, -1), (0, 1)]
            .into_iter()
            .filter_map(move |(di, dj)| {
                let (a, b) = (i + di, j + dj);
                (a >= 0 && b >= 0 && (a as usize) < w && (b as usize) < ht)
                    .then(|| b as usize * w + a as usize)
            })
    };

    // Empty cells between points get their neighbours' mean, so the surface has no pits.
    let mask: Vec<bool> = {
        let occ: Vec<bool> = z.iter().map(|v| !v.is_nan()).collect();
        let grown: Vec<bool> = (0..z.len())
            .map(|k| occ[k] || nb(k).any(|n| occ[n]))
            .collect();
        (0..z.len())
            .map(|k| grown[k] && (nb(k).all(|n| grown[n]) || occ[k]))
            .collect()
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
    // A few returns per cell are speckled; the median over the 3×3 neighbourhood is the roof's.
    let int: Vec<f64> = (0..z.len())
        .map(|k| {
            let (i, j) = ((k % w) as isize, (k / w) as isize);
            let mut v: Vec<f64> = (-1..=1)
                .flat_map(|dj| (-1..=1).map(move |di| (i + di, j + dj)))
                .filter(|&(a, b)| a >= 0 && b >= 0 && (a as usize) < w && (b as usize) < ht)
                .flat_map(|(a, b)| is[b as usize * w + a as usize].iter().copied())
                .collect();
            if v.is_empty() {
                return f64::NAN;
            }
            v.sort_by(f64::total_cmp);
            v[v.len() / 2]
        })
        .collect();

    // Local plane of each cell from its full 3×3 neighbourhood: slopes a, b (per m) and fit RMS.
    let at = |i: isize, j: isize| {
        (i >= 0 && j >= 0 && (i as usize) < w && (j as usize) < ht)
            .then(|| j as usize * w + i as usize)
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
                    rss +=
                        (v[dj][di] - (c + a * (di as f64 - 1.0) + b * (dj as f64 - 1.0))).powi(2);
                }
            }
            Some((a / res, b / res, (rss / 9.0).sqrt()))
        })
        .collect();
    let pos = |k: usize| [((k % w) as f64 + 0.5) * res, ((k / w) as f64 + 0.5) * res];

    // Roof planes grown from the flattest cells while the next cell stays on the plane.
    let mut label = vec![0usize; z.len()];
    let mut planes: Vec<Plane> = vec![Plane::default()];
    let mut order: Vec<usize> = (0..z.len())
        .filter(|&k| local[k].is_some_and(|l| l.2 < SEED_RMS))
        .collect();
    order.sort_by(|&a, &b| local[a].unwrap().2.total_cmp(&local[b].unwrap().2));
    for &s in &order {
        if label[s] != 0 {
            continue;
        }
        let r = planes.len();
        let (a, b, _) = local[s].unwrap();
        let ps = pos(s);
        let mut plane = Plane {
            a,
            b,
            c: z[s] - a * ps[0] - b * ps[1],
            ..Plane::default()
        };
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
                        plane = Plane {
                            a: p.0,
                            b: p.1,
                            c: p.2,
                            ..Plane::default()
                        };
                    }
                    next_refit *= 2;
                }
            }
        }
        if let Some(p) = fit.solve() {
            plane = Plane {
                a: p.0,
                b: p.1,
                c: p.2,
                ..Plane::default()
            };
        }
        plane.n = fit.n;
        planes.push(plane);
    }
    // Planes too small to be a roof face are dissolved; their cells, and ridge and edge cells no
    // plane took, go to the neighbouring plane they fit best.
    for l in label.iter_mut() {
        if *l != 0 && planes[*l].n < min_cells {
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
            if let Some(r) = nb(k)
                .map(|n| label[n])
                .filter(|&r| r != 0)
                .min_by(|&r, &s| {
                    (z[k] - planes[r].z(pk))
                        .abs()
                        .total_cmp(&(z[k] - planes[s].z(pk)).abs())
                })
            {
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
        let big: Vec<usize> = (0..pieces.len())
            .filter(|&p| pieces[p].len() >= min_cells)
            .collect();
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
                let f = if owner[k] == r {
                    0
                } else {
                    owner[k] - first_new + 1
                };
                fits[f].add(pos(k), z[k]);
                owner[k] = 0;
            }
        }
        for (f, fit) in fits.iter().enumerate() {
            let l = if f == 0 { r } else { first_new + f - 1 };
            if let Some((a, b, c)) = fit.solve() {
                planes[l] = Plane {
                    a,
                    b,
                    c,
                    n: fit.n,
                    ..Plane::default()
                };
            }
        }
    }
    Faces {
        x0,
        y0,
        w,
        ht,
        z,
        label,
        planes,
        int,
    }
}

/// Splits the points of one building blob into houses by their roof faces: faces meeting at a ridge
/// or hip are one roof; a valley, step or offset ridge between faces separates them. Parts smaller
/// than `min_part` (m²) join the neighbour they share most boundary with.
pub fn parts(
    points: &[P3],
    thetas: &[f64],
    res: f64,
    min_part: f64,
    debug: bool,
) -> (Vec<Vec<P3>>, Vec<Wall>) {
    let Faces {
        x0,
        y0,
        w,
        ht,
        z,
        mut label,
        planes,
        int,
    } = faces(points, res);
    let face = label.clone();
    let cell = |p: &P3| ((p[1] - y0) / res) as usize * w + ((p[0] - x0) / res) as usize;
    let nb = |k: usize| neighbours(k, w, ht);
    let at = |i: isize, j: isize| {
        (i >= 0 && j >= 0 && (i as usize) < w && (j as usize) < ht)
            .then(|| j as usize * w + i as usize)
    };
    let pos = |k: usize| [((k % w) as f64 + 0.5) * res, ((k / w) as f64 + 0.5) * res];
    // Median intensity of each label's cells; labels of unlike roofing are not one roof.
    let roofing = |label: &[usize]| {
        let mut v: BTreeMap<usize, Vec<f64>> = BTreeMap::new();
        for (k, &l) in label.iter().enumerate() {
            if l != 0 && !int[k].is_nan() {
                v.entry(l).or_default().push(int[k]);
            }
        }
        v.into_iter()
            .map(|(l, mut v)| {
                v.sort_by(f64::total_cmp);
                (l, v[v.len() / 2])
            })
            .collect::<BTreeMap<usize, f64>>()
    };
    let unlike = |m: &BTreeMap<usize, f64>, a: usize, b: usize| {
        m.get(&a)
            .zip(m.get(&b))
            .is_some_and(|(x, y)| (x - y).abs() > INT_TOL)
    };
    let face_roofing = roofing(&label);

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
            let Meeting {
                continuous,
                ridge,
                flat,
            } = meet(pr, ps, pk, pn, res);
            // Pieces of one face that grew apart.
            let coplanar = (pr.a - ps.a).abs() + (pr.b - ps.b).abs() < SLOPE_TOL;
            edge_pts.entry((r.min(s), r.max(s))).or_default().push(m);
            let e = votes.entry((r.min(s), r.max(s))).or_default();
            e.1 += 1;
            e.3 += continuous as usize;
            if continuous && (ridge || flat) && !pitches_differ(pr, ps) {
                e.0 += 1;
            } else if continuous && coplanar {
                e.2 += 1;
            }
        }
    }
    if debug {
        for (r, p) in planes.iter().enumerate().skip(1) {
            let v: Vec<[f64; 2]> = (0..z.len()).filter(|&k| label[k] == r).map(pos).collect();
            let n = v.len() as f64;
            let c = v
                .iter()
                .fold([x0, y0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n]);
            eprintln!(
                "face {r}: slope ({:.2}, {:.2}) cells {} at {:.1} {:.1}, intensity {:?}",
                p.a,
                p.b,
                p.n,
                c[0],
                c[1],
                face_roofing.get(&r)
            );
        }
        for ((r, s), (same, total, coplanar, _)) in &votes {
            eprintln!(
                "faces {r}-{s}: {same} of {total} boundary cells look like one roof, {coplanar} coplanar"
            );
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
    let mut joins: Vec<(usize, usize, bool)> = Vec::new();
    let mut support: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for (&(r, s), &(same, total, coplanar, cont)) in &votes {
        // Neighbours' coplanar faces touch at a corner or a neck; pieces of one face share much of
        // their outline.
        let outline = perim[&r].min(perim[&s]) as f64;
        let pieces = total >= MIN_COPLANAR_EDGE
            && (total as f64 >= COPLANAR_SHARE * outline || coplanar >= COPLANAR_RUN);
        let same = if pieces && !unlike(&face_roofing, r, s) {
            same + coplanar
        } else {
            same
        };
        let tilt = (planes[r].a - planes[s].a).abs() + (planes[r].b - planes[s].b).abs();
        // Coplanar faces meet along no line.
        let along = tilt < SLOPE_TOL || along_meet(&edge_pts[&(r, s)], &planes[r], &planes[s]);
        // A wing meets the main roof in a valley running diagonally across the walls; neighbouring
        // houses meet in one running along them.
        let wing = cont as f64 > SAME_ROOF_SHARE * total as f64
            && total >= MIN_COPLANAR_EDGE
            && diagonal(&edge_pts[&(r, s)], thetas)
            && along
            && !unlike(&face_roofing, r, s);
        if wing || (along && same >= MIN_JOIN_EDGE && same as f64 > SAME_ROOF_SHARE * total as f64)
        {
            if debug {
                let c = |f: usize| {
                    let v: Vec<[f64; 2]> =
                        (0..z.len()).filter(|&k| label[k] == f).map(pos).collect();
                    let n = v.len() as f64;
                    v.iter()
                        .fold([x0, y0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n])
                };
                let (cr, cs) = (c(r), c(s));
                eprintln!(
                    "join {r} ({:.1}, {:.1}) + {s} ({:.1}, {:.1}): {same}/{total}{}",
                    cr[0],
                    cr[1],
                    cs[0],
                    cs[1],
                    if wing { ", wing" } else { "" }
                );
            }
            // Only a ridge between two roofs is a bridge to judge; pieces of one face are one roof.
            let flat_pair = tilt < PIECE_TOL;
            joins.push((r, s, wing || flat_pair));
        }
        support.insert((r, s), if wing { total } else { same.min(total) });
    }
    // An annex's roof joined to a house's through one short ridge, but stepped from it along the rest
    // of their boundary, is a building of its own.
    let groups = |joins: &[(usize, usize, bool)], skip: usize| {
        let mut root: Vec<usize> = (0..planes.len()).collect();
        for (i, &(r, s, _)) in joins.iter().enumerate() {
            if i != skip {
                let (fr, fs) = (find(&mut root, r), find(&mut root, s));
                root[fr] = fs;
            }
        }
        (0..planes.len())
            .map(|f| find(&mut root, f))
            .collect::<Vec<_>>()
    };
    let mut unjoined = Vec::new();
    loop {
        let mut weakest: Option<(f64, usize, usize, usize)> = None;
        for i in (0..joins.len()).filter(|&i| !joins[i].2) {
            let g = groups(&joins, i);
            let (ga, gb) = (g[joins[i].0], g[joins[i].1]);
            if ga == gb {
                continue;
            }
            let (mut same, mut total) = (0, 0);
            for (&(r, s), &(_, t, _, cont)) in &votes {
                if (g[r] == ga && g[s] == gb) || (g[r] == gb && g[s] == ga) {
                    same += support[&(r, s)].max(cont);
                    total += t;
                }
            }
            let share = same as f64 / total as f64;
            // Cells, eaves and roof top of each side; the smaller must top out at the other's eaves.
            let side = |grp: usize| {
                let mut h: Vec<f64> = (0..z.len())
                    .filter(|&k| label[k] != 0 && g[label[k]] == grp && !z[k].is_nan())
                    .map(|k| z[k])
                    .collect();
                h.sort_by(f64::total_cmp);
                let at = |q: usize| h.get(h.len() * q / 100).copied().unwrap_or(f64::NAN);
                (h.len(), at(5), at(95))
            };
            let (a, b) = (side(ga), side(gb));
            let (small, large) = if a.0 < b.0 { (a, b) } else { (b, a) };
            if debug {
                eprintln!(
                    "bridge {} + {}: share {share:.2}, top {:.1} against eaves {:.1}",
                    joins[i].0, joins[i].1, small.2, large.1
                );
            }
            // A lower end topped well under the house's ridge is judged more readily, where a ridge
            // joins them; one slope running on into it is one roof.
            let (pr, ps) = (&planes[joins[i].0], &planes[joins[i].1]);
            let ridge = pr.a * ps.a + pr.b * ps.b < 0.0;
            let limit = if small.2 < large.1 + ANNEX_RISE {
                BRIDGE_SHARE
            } else if ridge && large.2 - small.2 > OBLIQUE_STEP {
                LOWER_BRIDGE_SHARE
            } else {
                continue;
            };
            if share < limit && weakest.is_none_or(|w| share < w.0) {
                weakest = Some((share, i, same, total));
            }
        }
        let Some((share, i, same, total)) = weakest else {
            break;
        };
        if debug {
            eprintln!(
                "unjoin {} + {}: {same} of {total} boundary cells between their roofs look like one or are continuous ({share:.2})",
                joins[i].0, joins[i].1
            );
        }
        unjoined.extend([joins[i].0, joins[i].1]);
        joins.remove(i);
    }
    for &(r, s, _) in &joins {
        let (fr, fs) = (find(&mut root, r), find(&mut root, s));
        root[fr] = fs;
    }
    for l in label.iter_mut() {
        if *l != 0 {
            *l = find(&mut root, *l);
        }
    }
    let unjoined: Vec<usize> = unjoined.into_iter().map(|f| find(&mut root, f)).collect();
    // Parts holding two faces that meet at a ridge: a gabled roof of their own.
    let ridged: Vec<usize> = joins
        .iter()
        .filter(|&&(r, s, _)| planes[r].a * planes[s].a + planes[r].b * planes[s].b < 0.0)
        .map(|&(r, _, _)| find(&mut root, r))
        .collect();

    // Slivers, and strips too narrow for a house (side hips, eaves, porch roofs), join the neighbour
    // they share most boundary with. A small roof well below any neighbour's is a shed or garage
    // against it and stays a building of its own, though wedged against a garage too.
    let cell_area = res * res;
    // Where along `b`'s side part `a` meets it: the extents of their shared edge and of that side.
    let side = |label: &[usize], a: usize, b: usize| {
        let mut edge = Vec::new();
        for k in (0..label.len()).filter(|&k| label[k] == a) {
            for n in nb(k).filter(|&n| label[n] == b) {
                let (p, q) = (pos(k), pos(n));
                edge.push([(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0]);
            }
        }
        let wall = Wall::fit(&edge)?;
        let (lo, hi) = (0..label.len())
            .filter(|&k| label[k] == b && wall.offset(&pos(k)).abs() < SIDE_REACH)
            .map(|k| wall.along(&pos(k)))
            .fold((f64::MAX, f64::MIN), |(lo, hi), t| (lo.min(t), hi.max(t)));
        Some((wall.t0, wall.t1, lo, hi))
    };
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
        let small_or_thin =
            |l: usize| (area[&l] as f64) * cell_area < min_part || width(&cells[&l]) < MIN_WIDTH;
        let Some((&small, _)) = area
            .iter()
            .filter(|&(&l, _)| !kept.contains(&l) && small_or_thin(l))
            .min_by_key(|&(_, &a)| a)
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
        let Some((&into, _)) = shared.iter().max_by_key(|&(_, &(c, _))| c) else {
            break;
        };
        let total: usize = shared.values().map(|&(c, _)| c).sum();
        let shed = shared
            .values()
            .any(|&(c, rise)| c as f64 >= SHED_EDGE * total as f64 && rise / c as f64 > SHED_STEP);
        // A roof of its own, stepped from the one it touches or with its own ridge, is a building
        // however small.
        let own = (unjoined.contains(&small) || ridged.contains(&small))
            && width(&cells[&small]) >= MIN_WIDTH;
        // A flat roof against a pitched one, well under its top: a garage or annex at the eaves.
        let top = |l: usize| {
            let mut h: Vec<f64> = (0..z.len())
                .filter(|&k| label[k] == l && !z[k].is_nan())
                .map(|k| z[k])
                .collect();
            h.sort_by(f64::total_cmp);
            h.get(h.len() * 95 / 100).copied().unwrap_or(f64::NAN)
        };
        let pitch = |l: usize| planes[l].a.hypot(planes[l].b);
        let flat = pitch(small) < ANNEX_SLOPE
            && pitch(into) > ANNEX_SLOPE
            && top(into) - top(small) > OBLIQUE_STEP
            && area[&small] as f64 * cell_area >= MIN_GARAGE;
        if debug && pitch(small) < ANNEX_SLOPE {
            eprintln!(
                "flat {small}: top {:.1}, {into} pitched {:.2}, top {:.1}",
                top(small),
                pitch(into),
                top(into)
            );
        }
        if (shed || own || flat) && area[&small] as f64 * cell_area >= MIN_SHED {
            if debug {
                eprintln!(
                    "small {small} ({} cells), neighbours {shared:?}: a shed or a roof of its own",
                    area[&small]
                );
            }
            kept.insert(small);
            continue;
        }
        if debug {
            eprintln!(
                "small or thin {small} ({} cells) into {into}, neighbours {shared:?}",
                area[&small]
            );
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
                    let c = |k: usize| {
                        [
                            x0 + ((k % w) as f64 + 0.5) * res,
                            y0 + ((k / w) as f64 + 0.5) * res,
                        ]
                    };
                    let (a, b) = (c(k), c(n));
                    pairs
                        .entry((label[k], label[n]))
                        .or_default()
                        .push([(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0]);
                }
            }
        }
        pairs
    };

    // A party wall runs along the houses' axes; a split oblique to them is a roof split in two.
    // An annex whose roof stays at the house's eaves stays apart whatever its boundary: an
    // L-shaped one fits an oblique line. But one continuing the house's walls and roof extends it, as
    // does a smaller wing pitched like the house, narrower and with a lower ridge.
    loop {
        let mut zs: BTreeMap<usize, Vec<f64>> = BTreeMap::new();
        let mut area: BTreeMap<usize, f64> = BTreeMap::new();
        for (k, &l) in label.iter().enumerate() {
            if l != 0 {
                *area.entry(l).or_default() += cell_area;
                if !z[k].is_nan() {
                    zs.entry(l).or_default().push(z[k]);
                }
            }
        }
        let pct = |l: usize, q: f64| {
            let mut v = zs.get(&l).cloned().unwrap_or_default();
            v.sort_by(f64::total_cmp);
            if v.is_empty() {
                f64::NAN
            } else {
                v[((v.len() - 1) as f64 * q) as usize]
            }
        };
        let below = |a: usize, b: usize| pct(a, 0.95) < pct(b, 0.05) + ANNEX_RISE;
        let bounds = boundaries(&label);
        let part_roofing = roofing(&label);
        let mut neighbours: BTreeMap<usize, usize> = BTreeMap::new();
        for &(p, q) in bounds.keys() {
            *neighbours.entry(p).or_default() += 1;
            *neighbours.entry(q).or_default() += 1;
        }
        // Share of the edge between `a` and `b` along which their faces slope alike: one roof, offset.
        let parallel = |a: usize, b: usize| {
            let (mut same, mut all) = (0, 0);
            for k in (0..label.len()).filter(|&k| label[k] == a) {
                for n in nb(k).filter(|&n| label[n] == b) {
                    let (pr, ps) = (&planes[face[k]], &planes[face[n]]);
                    all += 1;
                    if face[k] != 0
                        && face[n] != 0
                        && pr.a.hypot(pr.b) > ANNEX_SLOPE
                        && (pr.a - ps.a).abs() + (pr.b - ps.b).abs() < SLOPE_TOL
                    {
                        same += 1;
                    }
                }
            }
            all > 0 && same as f64 >= SAME_ROOF_SHARE * all as f64
        };
        let Some((lp, lq)) = bounds
            .into_iter()
            .find(|&((p, q), ref b)| {
                // A wing has no other neighbour.
                let wing = |a: usize, b: usize| {
                    if neighbours[&a] > 1 || unjoined.contains(&a) || unlike(&part_roofing, a, b) {
                        return false;
                    }
                    let Some((t0, t1, lo, hi)) = side(&label, a, b) else {
                        return false;
                    };
                    let flush = t1 - t0 >= FLUSH_SHARE * (hi - lo);
                    if debug {
                        eprintln!(
                            "{a} against {b}: edge {:.1} to {:.1} of side {:.1} to {:.1}",
                            t0 - lo,
                            t1 - lo,
                            0.0,
                            hi - lo
                        );
                    }
                    area[&a] < WING_SHARE * area[&b]
                        && !below(a, b)
                        && pct(a, 0.95) < pct(b, 0.95)
                        && !flush
                        && parallel(a, b)
                };
                let oblique = Wall::fit(b).is_some_and(|w| w.axis(thetas).is_none())
                    && !stepped(b, thetas)
                    && !below(p, q)
                    && !below(q, p)
                    && (pct(p, 0.95) - pct(q, 0.95)).abs() < OBLIQUE_STEP
                    && !unlike(&part_roofing, p, q);
                oblique || wing(p, q) || wing(q, p)
            })
            .map(|(k, _)| k)
        else {
            break;
        };
        if debug {
            eprintln!(
                "oblique split or wing: {lp} (roof {:.1} to {:.1}, intensity {:?}) + {lq} ({:.1} to {:.1}, {:?})",
                pct(lp, 0.05),
                pct(lp, 0.95),
                part_roofing.get(&lp),
                pct(lq, 0.05),
                pct(lq, 0.95),
                part_roofing.get(&lq)
            );
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
        let Some(mut wall) = Wall::fit(&b) else {
            if debug {
                eprintln!(
                    "boundary {lp}-{lq}: {} cells, too short for a wall",
                    b.len()
                );
            }
            continue;
        };
        // A stepped boundary fits a skewed line; the wall is its longest run.
        let (axis, b) = match wall.axis(thetas) {
            Some(axis) => (axis, b),
            None if stepped(&b, thetas) => {
                let Some((axis, run)) = runs(&b, thetas).into_iter().next() else {
                    continue;
                };
                if debug {
                    eprintln!(
                        "boundary {lp}-{lq}: stepped, longest run {} of {} cells at {:.0}°",
                        run.len(),
                        b.len(),
                        axis.to_degrees()
                    );
                }
                // Its position from the run, its extent from the whole boundary.
                let n = run.len() as f64;
                wall.p = run
                    .iter()
                    .fold([0.0, 0.0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n]);
                (axis, b)
            }
            None => {
                if debug {
                    eprintln!(
                        "boundary {lp}-{lq}: not stepped: runs {:?}",
                        runs(&b, thetas)
                            .iter()
                            .map(|(a, r)| (a.to_degrees().round(), r.len()))
                            .collect::<Vec<_>>()
                    );
                    eprintln!(
                        "boundary {lp}-{lq}: at {:.0}° off the axes {thetas:?}",
                        wall.dir[1].atan2(wall.dir[0]).to_degrees()
                    );
                }
                continue;
            }
        };
        wall.set_dir(axis, &b);
        // Which side is `lp`: the mean offset of its points near the wall.
        let near = |p: &P3| {
            wall.offset(p).abs() < WALL_REACH
                && wall.along(p) > wall.t0 - 1.0
                && wall.along(p) < wall.t1 + 1.0
        };
        let side_p: f64 = points
            .iter()
            .zip(&point_label)
            .filter(|(p, l)| **l == lp && near(p))
            .map(|(p, _)| wall.offset(p))
            .sum();
        for (p, l) in points.iter().zip(point_label.iter_mut()) {
            if (*l == lp || *l == lq) && near(p) {
                *l = if (wall.offset(p) > 0.0) == (side_p > 0.0) {
                    lp
                } else {
                    lq
                };
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
        let c = b
            .iter()
            .fold([0.0, 0.0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n]);
        let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
        for p in b {
            let (dx, dy) = (p[0] - c[0], p[1] - c[1]);
            sxx += dx * dx;
            sxy += dx * dy;
            syy += dy * dy;
        }
        let ang = 0.5 * (2.0 * sxy).atan2(sxx - syy);
        let mut wall = Wall {
            p: c,
            dir: [ang.cos(), ang.sin()],
            t0: f64::MAX,
            t1: f64::MIN,
        };
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
/// Faces falling more nearly opposite than this (cosine) are the two sides of one ridge.
const OPPOSITE_COS: f64 = 0.9;
/// Ratio of pitch angles beyond which such faces are of two roofs.
const PITCH_RATIO: f64 = 1.6;
/// Faces or parts whose median intensities differ by more are of unlike roofing.
pub const INT_TOL: f64 = 6000.0;
/// Faces smaller than this (m²) are dissolved into their neighbours.
const MIN_PLANE_AREA: f64 = 5.0;
/// Two faces are continuous where they meet within this height (m).
const JOIN_TOL: f64 = 0.3;
/// Share of the smaller face's outline that coplanar faces must share to be pieces of one face.
const COPLANAR_SHARE: f64 = 0.25;

/// Fewest boundary cells that look like one roof for faces to join: a few where faces touch over a
/// step are noise.
const MIN_JOIN_EDGE: usize = 6;
/// Coplanar faces continuing each other along this many boundary cells are one face, however large.
const COPLANAR_RUN: usize = 20;
/// Shortest boundary (cells) along which coplanar faces are taken for pieces of one face.
const MIN_COPLANAR_EDGE: usize = 6;

/// How much two faces extended over each other must rise above it to meet in a ridge or hip (m).
const CONVEX_EPS: f64 = 0.05;

/// Slopes below this (per m) are flat.
const FLAT: f64 = 0.05;
/// Share of the boundary that must look like one roof for two faces to be joined.
const SAME_ROOF_SHARE: f64 = 0.6;
/// Share of the whole boundary between two roofs joined by a single pair of faces that must look
/// like one roof.
const BRIDGE_SHARE: f64 = 0.3;
/// The same for a part topped more than OBLIQUE_STEP under the other.
const LOWER_BRIDGE_SHARE: f64 = 0.4;
/// Largest slope difference (sum of both components, per m) of two pieces of one face, tighter than
/// [`SLOPE_TOL`] so low roofs falling different ways are not taken for one.
const PIECE_TOL: f64 = 0.15;

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
        for (s, v) in self
            .s
            .iter_mut()
            .zip([x, y, z, x * x, x * y, y * y, x * z, y * z, 1.0])
        {
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
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
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

/// A small part whose roof is this much below a neighbour's along their edge is its own building (m) ...
const SHED_STEP: f64 = 1.0;
/// ... if that edge is at least this share of the part's edges with others.
const SHED_EDGE: f64 = 0.3;
/// A part whose roof tops out less than this (m) above another's eaves is an annex.
const ANNEX_RISE: f64 = 0.5;
/// A part cut off along a hip line reaches nearly as high as the rest; one topped this much lower
/// is a roof of its own.
const OBLIQUE_STEP: f64 = 1.5;
/// Least area of a flat annex (m²); smaller flat patches are porches and terraces.
const MIN_GARAGE: f64 = 15.0;
/// ... if at least this large (m²).
const MIN_SHED: f64 = 8.0;
/// A part sharing this much of its neighbour's side is flush with it, not a wing.
const FLUSH_SHARE: f64 = 0.8;
/// A part's side is its cells this close to the edge it shares with another (m).
const SIDE_REACH: f64 = 1.0;
/// A wing, narrower than the house, is at most this share of its area.
const WING_SHARE: f64 = 0.3;
/// A wing's roof is pitched at least this much (per m).
const ANNEX_SLOPE: f64 = 0.15;

/// Narrowest house split off a blob (m).
const MIN_WIDTH: f64 = 3.0;

/// Extent of cell centres across their principal axis.
fn width(c: &[[f64; 2]]) -> f64 {
    let n = c.len() as f64;
    let m = c
        .iter()
        .fold([0.0, 0.0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n]);
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
    let Some(t) = thetas.first() else {
        return false;
    };
    let d = (direction(b) - t).rem_euclid(FRAC_PI_2);
    (d - FRAC_PI_2 / 2.0).abs() < DIAGONAL_TOL.to_radians()
}

/// Whether a boundary runs along the line where the two faces' planes meet, as a ridge, hip or
/// valley does; faces of two roofs crossing at their ends meet across it.
fn along_meet(b: &[[f64; 2]], pr: &Plane, ps: &Plane) -> bool {
    let dir = direction(b);
    let (sin, cos) = dir.sin_cos();
    let along = b.iter().map(|p| p[0] * cos + p[1] * sin);
    let run = along.clone().fold(f64::MIN, f64::max) - along.fold(f64::MAX, f64::min);
    // A boundary this short has no direction to judge.
    if run < MIN_MEET_RUN {
        return true;
    }
    let line = (pr.a - ps.a).atan2(ps.b - pr.b);
    let d = (dir - line).rem_euclid(PI);
    d.min(PI - d) < MEET_TOL.to_radians()
}

/// The main direction of points `b` (radians).
fn direction(b: &[[f64; 2]]) -> f64 {
    let n = b.len() as f64;
    let c = b
        .iter()
        .fold([0.0, 0.0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n]);
    let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
    for p in b {
        let (dx, dy) = (p[0] - c[0], p[1] - c[1]);
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    0.5 * (2.0 * sxy).atan2(sxx - syy)
}

/// Whether the boundary `b` runs along the walls in at most two straight runs, as between
/// neighbours set back from each other, rather than across them as a valley of one roof.
fn stepped(b: &[[f64; 2]], thetas: &[f64]) -> bool {
    let on_runs: usize = runs(b, thetas).iter().map(|(_, r)| r.len()).sum();
    on_runs as f64 >= STEP_SHARE * b.len() as f64
}

/// The two straight runs along the walls that take most of the boundary `b`, largest first: each
/// one's direction (radians) and points.
fn runs(b: &[[f64; 2]], thetas: &[f64]) -> Vec<(f64, Vec<[f64; 2]>)> {
    let Some(&t) = thetas.first() else {
        return Vec::new();
    };
    let (sin, cos) = t.sin_cos();
    let mut left: Vec<([f64; 2], [f64; 2])> = b
        .iter()
        .map(|p| ([p[0] * cos + p[1] * sin, -p[0] * sin + p[1] * cos], *p))
        .collect();
    let mut out = Vec::new();
    for _ in 0..2 {
        // The line across either axis that takes most points.
        let best = (0..2)
            .flat_map(|ax| left.iter().map(move |(q, _)| (ax, q[ax])))
            .max_by_key(|&(ax, v)| {
                left.iter()
                    .filter(|(q, _)| (q[ax] - v).abs() < STEP_TOL)
                    .count()
            });
        let Some((ax, v)) = best else { break };
        let (on, off): (Vec<_>, Vec<_>) = left
            .into_iter()
            .partition(|(q, _)| (q[ax] - v).abs() < STEP_TOL);
        // Points at one place across the first axis lie on a line along the second.
        out.push((
            if ax == 0 { t + FRAC_PI_2 } else { t },
            on.into_iter().map(|(_, p)| p).collect(),
        ));
        left = off;
    }
    out
}

/// Boundary points within this of a run along the walls lie on it (m) ...
const STEP_TOL: f64 = 0.5;
/// ... and a boundary with this share of its points on two runs is stepped.
const STEP_SHARE: f64 = 0.8;

/// How far from 45° to the walls a valley may run and still join a wing (degrees).
const DIAGONAL_TOL: f64 = 15.0;
/// How far a ridge, hip or valley may run from the line where its faces' planes meet (degrees).
const MEET_TOL: f64 = 30.0;
/// The shortest boundary (m) whose direction is judged against the faces' meeting line.
const MIN_MEET_RUN: f64 = 1.0;

/// The 4-neighbours of cell `k` in a `w` × `ht` grid.
pub(crate) fn neighbours(k: usize, w: usize, ht: usize) -> impl Iterator<Item = usize> {
    let (i, j) = ((k % w) as isize, (k / w) as isize);
    [(-1, 0), (1, 0), (0, -1), (0, 1)]
        .into_iter()
        .filter_map(move |(di, dj)| {
            let (a, b) = (i + di, j + dj);
            (a >= 0 && b >= 0 && (a as usize) < w && (b as usize) < ht)
                .then(|| b as usize * w + a as usize)
        })
}

/// How two faces meet across the boundary between neighbouring cells `pk` (face `pr`) and `pn`.
pub(crate) struct Meeting {
    pub continuous: bool,
    /// A ridge or hip, of any direction, rather than a valley.
    pub ridge: bool,
    pub flat: bool,
}

/// Opposite faces of one gabled roof are pitched alike; a ridge between unlike pitches is where a
/// lower or steeper roof abuts another.
fn pitches_differ(pr: &Plane, ps: &Plane) -> bool {
    let (nr, ns) = (pr.a.hypot(pr.b), ps.a.hypot(ps.b));
    let opposite = pr.a * ps.a + pr.b * ps.b < -OPPOSITE_COS * nr * ns;
    // Compared as angles: steep roofs differ more in slope for the same difference in look.
    let (ar, as_) = (nr.atan(), ns.atan());
    opposite && ar.max(as_) > PITCH_RATIO * ar.min(as_)
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
