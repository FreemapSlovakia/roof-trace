//! Outlines of all the houses of one blob at once, as cells of one right-angled partition: lines
//! along the blob's walls cut it into blocks, and each block goes to one house or to none. Houses
//! that touch then share their walls exactly, and every outline is right-angled.

use crate::outline::Grid;
use crate::store::Pt;

type P = [f64; 2];

/// Raster cell (m).
const RES: f64 = 0.25;
/// A line is a candidate where labels change along at least this much of it (m).
const MIN_RUN: f64 = 1.0;
/// A strip without roof points at least this wide (m) between two houses is a gap between them.
const MIN_GAP: f64 = 0.75;
/// A line stays only if removing it would misassign more than this area (m²).
const LINE_COST: f64 = 1.0;
/// Cost of a corner in an outline, in misassigned m².
const CORNER_COST: f64 = 3.0;
/// Blobs whose raster would need more label counts than this are left to the per-house outlines.
const MAX_COUNTS: usize = 20_000_000;

/// A blob's outlines from the partition, and what they and the houses' own outlines cost.
pub struct Partition {
    /// Per house, its rings, outer first, in map units.
    pub rings: Vec<Option<Vec<Vec<P>>>>,
    /// Misassigned roof (m²) plus `CORNER_COST` per corner, of the partition and of the own outlines.
    pub cost: f64,
    pub own_cost: f64,
}

/// The partition of the houses `parts` for walls along `theta`, measured against their own
/// outlines `own`. `all` are the blob's points, some of them in no part; `cover` is canopy over the
/// roof: what lies under it is not known to be empty, unless bare ground was seen there.
pub fn partition(
    parts: &[&[Pt]],
    own: &[&Vec<Vec<P>>],
    all: &[P],
    cover: &[P],
    bare: &[P],
    theta: f64,
) -> Option<Partition> {
    let n_parts = parts.len();
    let in_parts = parts.iter().flat_map(|p| p.iter());
    let n = parts.iter().map(|p| p.len()).sum::<usize>() as f64;
    if n == 0.0 {
        return None;
    }
    let c = in_parts.fold([0.0, 0.0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n]);
    let (sin, cos) = theta.sin_cos();
    let rot = |p: &[f64]| {
        let (dx, dy) = (p[0] - c[0], p[1] - c[1]);
        [dx * cos + dy * sin, -dx * sin + dy * cos]
    };
    let unrot = |p: P| {
        [
            c[0] + p[0] * cos - p[1] * sin,
            c[1] + p[0] * sin + p[1] * cos,
        ]
    };

    let rotated: Vec<Vec<Pt>> = parts
        .iter()
        .map(|ps| {
            ps.iter()
                .map(|p| {
                    let r = rot(p);
                    [r[0], r[1], p[2], f64::NAN]
                })
                .collect()
        })
        .collect();
    let flat: Vec<Pt> = all
        .iter()
        .map(|p| {
            let r = rot(p);
            [r[0], r[1], 0.0, f64::NAN]
        })
        .collect();
    let mut q = [f64::MAX, f64::MAX, f64::MIN, f64::MIN];
    for p in rotated.iter().flatten().chain(&flat) {
        q = [
            q[0].min(p[0]),
            q[1].min(p[1]),
            q[2].max(p[0]),
            q[3].max(p[1]),
        ];
    }
    let q = [q[0] - 2.0, q[1] - 2.0, q[2] + 2.0, q[3] + 2.0];
    let grid = Grid::new(q, RES, &flat, &[]);
    let (x0, y0, res, w, h, occ) = grid.mask();
    if w * h * (n_parts + 1) > MAX_COUNTS {
        return None;
    }
    let empty = n_parts;
    let unknown = n_parts + 1;
    let cell = |p: &[f64]| {
        let (i, j) = (((p[0] - x0) / res) as usize, ((p[1] - y0) / res) as usize);
        (i < w && j < h).then(|| j * w + i)
    };

    // Each building cell's part: the one with most points in it, else the nearest.
    let mut counts = vec![0u32; w * h * n_parts];
    for (l, ps) in rotated.iter().enumerate() {
        for p in ps {
            if let Some(k) = cell(p) {
                counts[k * n_parts + l] += 1;
            }
        }
    }
    let mut label = vec![empty; w * h];
    let mut queue = std::collections::VecDeque::new();
    for k in 0..w * h {
        let row = &counts[k * n_parts..(k + 1) * n_parts];
        let (l, &m) = row.iter().enumerate().max_by_key(|&(_, &m)| m).unwrap();
        if occ[k] && m > 0 {
            label[k] = l;
            queue.push_back(k);
        }
    }
    let mut seen: Vec<bool> = (0..w * h).map(|k| label[k] != empty).collect();
    while let Some(k) = queue.pop_front() {
        for m in neighbours(k, w, h) {
            if occ[m] && !seen[m] {
                seen[m] = true;
                label[m] = label[k];
                queue.push_back(m);
            }
        }
    }
    // Under canopy a roof may go on unseen; ground seen through it is a gap. Ground is sparse
    // there, so one point clears its neighbours too.
    let mut ground = vec![false; w * h];
    for p in bare {
        if let Some(k) = cell(&rot(p)) {
            ground[k] = true;
            for m in neighbours(k, w, h) {
                ground[m] = true;
            }
        }
    }
    // A strip without roof points between two houses, along the walls and with ground seen in it, is
    // a gap the mask's closing bridged: they stand apart.
    let has_pts: Vec<bool> = (0..w * h)
        .map(|k| {
            counts[k * n_parts..(k + 1) * n_parts]
                .iter()
                .any(|&c| c > 0)
        })
        .collect();
    let min_gap = (MIN_GAP / res).ceil() as usize;
    let mut gap = vec![false; w * h];
    for (len, step, lines, line_step) in [(w, 1, h, w), (h, w, w, 1)] {
        for line in 0..lines {
            let at = |t: usize| line * line_step + t * step;
            let mut t = 0;
            while t < len {
                if has_pts[at(t)] || !occ[at(t)] {
                    t += 1;
                    continue;
                }
                let start = t;
                while t < len && !has_pts[at(t)] && occ[at(t)] {
                    t += 1;
                }
                let ends = (start > 0 && t < len).then(|| (at(start - 1), at(t)));
                let apart =
                    ends.is_some_and(|(a, b)| has_pts[a] && has_pts[b] && label[a] != label[b]);
                if apart && t - start >= min_gap && (start..t).any(|u| ground[at(u)]) {
                    for u in start..t {
                        gap[at(u)] = true;
                    }
                }
            }
        }
    }
    for k in 0..w * h {
        if gap[k] {
            label[k] = empty;
        }
    }
    for p in cover {
        if let Some(k) = cell(&rot(p)) {
            if label[k] == empty && !ground[k] {
                label[k] = unknown;
            }
        }
    }

    // Per label (parts and empty) prefix sums over the raster, for any block's counts at once.
    let nl = n_parts + 1;
    let (pw, ph) = (w + 1, h + 1);
    let mut sum = vec![0u32; nl * pw * ph];
    for j in 0..h {
        for i in 0..w {
            for l in 0..nl {
                let here = (label[j * w + i] == l) as u32;
                let at = |i: usize, j: usize| (l * ph + j) * pw + i;
                sum[at(i + 1, j + 1)] =
                    here + sum[at(i, j + 1)] + sum[at(i + 1, j)] - sum[at(i, j)];
            }
        }
    }
    let count = |l: usize, i0: usize, i1: usize, j0: usize, j1: usize| {
        let at = |i: usize, j: usize| sum[(l * ph + j) * pw + i];
        at(i1, j1) + at(i0, j0) - at(i0, j1) - at(i1, j0)
    };
    // A block's majority label and how many of its known cells it misassigns.
    let block = |i0: usize, i1: usize, j0: usize, j1: usize| {
        let (mut best, mut most, mut known) = (empty, 0u32, 0u32);
        for l in 0..nl {
            let c = count(l, i0, i1, j0, j1);
            known += c;
            if c > most {
                (best, most) = (l, c);
            }
        }
        (best, known - most)
    };

    // Candidate lines: cell boundaries along which the label changes over MIN_RUN.
    let min_run = (MIN_RUN / res).ceil() as usize;
    let known = |l: usize| l != unknown;
    let mut xs: Vec<usize> = (1..w)
        .filter(|&i| {
            (0..h)
                .filter(|&j| {
                    let (a, b) = (label[j * w + i - 1], label[j * w + i]);
                    a != b && known(a) && known(b)
                })
                .count()
                >= min_run
        })
        .collect();
    let mut ys: Vec<usize> = (1..h)
        .filter(|&j| {
            (0..w)
                .filter(|&i| {
                    let (a, b) = (label[(j - 1) * w + i], label[j * w + i]);
                    a != b && known(a) && known(b)
                })
                .count()
                >= min_run
        })
        .collect();
    xs.insert(0, 0);
    xs.push(w);
    ys.insert(0, 0);
    ys.push(h);

    // Lines go, cheapest first, while removing one misassigns less than LINE_COST.
    let line_cost = LINE_COST / (res * res);
    loop {
        let strip_err = |a: &[usize], b: &[usize], k: usize, across: bool| -> f64 {
            // Removing a[k] merges the strips a[k-1]..a[k] and a[k]..a[k+1] in every strip of b.
            (0..b.len() - 1)
                .map(|m| {
                    let (b0, b1) = (b[m], b[m + 1]);
                    let e = |i0: usize, i1: usize| {
                        if across {
                            block(b0, b1, i0, i1).1
                        } else {
                            block(i0, i1, b0, b1).1
                        }
                    };
                    e(a[k - 1], a[k + 1]) as f64
                        - e(a[k - 1], a[k]) as f64
                        - e(a[k], a[k + 1]) as f64
                })
                .sum()
        };
        let best_x = (1..xs.len() - 1)
            .map(|k| (strip_err(&xs, &ys, k, false), k))
            .min_by(|a, b| a.0.total_cmp(&b.0));
        let best_y = (1..ys.len() - 1)
            .map(|k| (strip_err(&ys, &xs, k, true), k))
            .min_by(|a, b| a.0.total_cmp(&b.0));
        match (best_x, best_y) {
            (Some((dx, kx)), by) if dx < line_cost && by.is_none_or(|(dy, _)| dx <= dy) => {
                xs.remove(kx);
            }
            (_, Some((dy, ky))) if dy < line_cost => {
                ys.remove(ky);
            }
            _ => break,
        }
    }

    // Block labels by majority, then improved by moves that lower the misassigned area plus
    // CORNER_COST per corner: a step or bump stays only where it explains enough roof.
    let (nx, ny) = (xs.len() - 1, ys.len() - 1);
    let mut lab: Vec<usize> = (0..nx * ny)
        .map(|b| block(xs[b % nx], xs[b % nx + 1], ys[b / nx], ys[b / nx + 1]).0)
        .collect();
    let miss = |b: usize, l: usize| {
        let (i0, i1, j0, j1) = (xs[b % nx], xs[b % nx + 1], ys[b / nx], ys[b / nx + 1]);
        let total: u32 = (0..nl).map(|m| count(m, i0, i1, j0, j1)).sum();
        (total - count(l, i0, i1, j0, j1)) as f64
    };
    let corner_cost = CORNER_COST / (res * res);
    // Corners of the houses' outlines at grid vertex (vi, vj); beyond the grid is empty.
    let corners_at = |lab: &[usize], vi: usize, vj: usize| {
        let get = |i: isize, j: isize| {
            if i < 0 || j < 0 || i as usize >= nx || j as usize >= ny {
                empty
            } else {
                lab[j as usize * nx + i as usize]
            }
        };
        let (i, j) = (vi as isize, vj as isize);
        let around = [get(i - 1, j - 1), get(i, j - 1), get(i, j), get(i - 1, j)];
        let mut n = 0;
        for (k, &l) in around.iter().enumerate() {
            if l == empty || around[..k].contains(&l) {
                continue;
            }
            let has: Vec<bool> = around.iter().map(|&m| m == l).collect();
            let c = has.iter().filter(|&&v| v).count();
            // One or three of the four blocks, or two diagonally: the outline turns here.
            if c == 1 || c == 3 || (c == 2 && has[0] == has[2]) {
                n += 1;
            }
        }
        n
    };
    // Cost change of giving `blocks` label `to`.
    let delta = |lab: &mut Vec<usize>, blocks: &[usize], to: usize| {
        let mut verts: Vec<(usize, usize)> = blocks
            .iter()
            .flat_map(|&b| {
                let (bi, bj) = (b % nx, b / nx);
                [(bi, bj), (bi + 1, bj), (bi, bj + 1), (bi + 1, bj + 1)]
            })
            .collect();
        verts.sort_unstable();
        verts.dedup();
        let corners = |lab: &[usize]| {
            verts
                .iter()
                .map(|&(i, j)| corners_at(lab, i, j))
                .sum::<usize>() as f64
        };
        let before =
            corners(lab) * corner_cost + blocks.iter().map(|&b| miss(b, lab[b])).sum::<f64>();
        let old: Vec<usize> = blocks.iter().map(|&b| lab[b]).collect();
        for &b in blocks {
            lab[b] = to;
        }
        let after = corners(lab) * corner_cost + blocks.iter().map(|&b| miss(b, to)).sum::<f64>();
        for (&b, &l) in blocks.iter().zip(&old) {
            lab[b] = l;
        }
        after - before
    };
    for _ in 0..50 {
        let mut moves: Vec<(Vec<usize>, usize)> = Vec::new();
        for b in 0..nx * ny {
            let (bi, bj) = (b % nx, b / nx);
            let mut to: Vec<usize> = [
                (bi > 0).then(|| b - 1),
                (bi + 1 < nx).then(|| b + 1),
                (bj > 0).then(|| b - nx),
                (bj + 1 < ny).then(|| b + nx),
            ]
            .into_iter()
            .flatten()
            .map(|m| lab[m])
            .collect();
            to.push(empty);
            moves.extend(
                to.into_iter()
                    .filter(|&l| l != lab[b])
                    .map(|l| (vec![b], l)),
            );
        }
        // Runs of like blocks along a strip, flipped whole to the label beside them all: a step one
        // block at a time would first add corners.
        for (len, n, step, side) in [(ny, nx, nx, 1usize), (nx, ny, 1, nx)] {
            for s in 0..n {
                let at = |k: usize| s * side + k * step;
                let mut k = 0;
                while k < len {
                    let l = lab[at(k)];
                    let mut e = k;
                    while e + 1 < len && lab[at(e + 1)] == l {
                        e += 1;
                    }
                    let run: Vec<usize> = (k..=e).map(at).collect();
                    for beside in [s.checked_sub(1), (s + 1 < n).then_some(s + 1)] {
                        let to = match beside {
                            Some(t) => {
                                let ls: Vec<usize> =
                                    (k..=e).map(|m| lab[t * side + m * step]).collect();
                                if ls.iter().all(|&x| x == ls[0]) {
                                    ls[0]
                                } else {
                                    continue;
                                }
                            }
                            None => empty,
                        };
                        if to != l {
                            moves.push((run.clone(), to));
                        }
                    }
                    k = e + 1;
                }
            }
        }
        let best = moves
            .into_iter()
            .map(|(blocks, to)| (delta(&mut lab, &blocks, to), blocks, to))
            .min_by(|a, b| a.0.total_cmp(&b.0));
        match best {
            Some((d, blocks, to)) if d < -1e-9 => {
                for b in blocks {
                    lab[b] = to;
                }
            }
            _ => break,
        }
    }

    let misassigned: u32 = (0..nx * ny)
        .map(|b| {
            let (i0, i1, j0, j1) = (xs[b % nx], xs[b % nx + 1], ys[b / nx], ys[b / nx + 1]);
            let total: u32 = (0..nl).map(|l| count(l, i0, i1, j0, j1)).sum();
            total - count(lab[b], i0, i1, j0, j1)
        })
        .sum();
    let corners: usize = (0..=nx)
        .flat_map(|i| (0..=ny).map(move |j| (i, j)))
        .map(|(i, j)| corners_at(&lab, i, j))
        .sum();
    let cost = misassigned as f64 * res * res + CORNER_COST * corners as f64;

    // The own outlines on the same raster; a cell two of them hold is misassigned by one.
    let own_rot: Vec<Vec<Vec<P>>> = own
        .iter()
        .map(|rings| {
            rings
                .iter()
                .map(|r| r.iter().map(|p| rot(p)).collect())
                .collect()
        })
        .collect();
    let mut own_miss = 0u32;
    for k in 0..w * h {
        if label[k] == unknown {
            continue;
        }
        let p = [
            x0 + ((k % w) as f64 + 0.5) * res,
            y0 + ((k / w) as f64 + 0.5) * res,
        ];
        let holders: Vec<usize> = (0..own_rot.len())
            .filter(|&l| inside(&own_rot[l][0], p) && !own_rot[l][1..].iter().any(|h| inside(h, p)))
            .collect();
        own_miss += match holders[..] {
            [] => (label[k] != empty) as u32,
            [l] => (l != label[k]) as u32,
            _ => 1,
        };
    }
    let own_corners: usize = own
        .iter()
        .flat_map(|rings| rings.iter())
        .map(Vec::len)
        .sum();
    let own_cost = own_miss as f64 * res * res + CORNER_COST * own_corners as f64;

    let rings = (0..n_parts)
        .map(|l| {
            let rings = trace(&lab, nx, ny, l);
            (!rings.is_empty()).then(|| {
                rings
                    .into_iter()
                    .map(|r| {
                        r.into_iter()
                            .map(|(i, j)| unrot([x0 + xs[i] as f64 * res, y0 + ys[j] as f64 * res]))
                            .collect()
                    })
                    .collect()
            })
        })
        .collect();
    Some(Partition {
        rings,
        cost,
        own_cost,
    })
}

/// Whether `p` lies inside ring `r`.
fn inside(r: &[P], p: P) -> bool {
    let n = r.len();
    let mut ins = false;
    for i in 0..n {
        let (a, b) = (r[i], r[(i + 1) % n]);
        if (a[1] > p[1]) != (b[1] > p[1])
            && p[0] < a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0])
        {
            ins = !ins;
        }
    }
    ins
}

fn neighbours(k: usize, w: usize, h: usize) -> impl Iterator<Item = usize> {
    let (i, j) = (k % w, k / w);
    [
        (i > 0).then(|| k - 1),
        (i + 1 < w).then(|| k + 1),
        (j > 0).then(|| k - w),
        (j + 1 < h).then(|| k + w),
    ]
    .into_iter()
    .flatten()
}

/// Rings around the largest 4-connected region of blocks labelled `l`, outer first, as block
/// corners; holes follow, clockwise.
fn trace(lab: &[usize], nx: usize, ny: usize, l: usize) -> Vec<Vec<(usize, usize)>> {
    // The largest region only: a stray block elsewhere is noise.
    let mut region = vec![usize::MAX; lab.len()];
    let mut sizes = Vec::new();
    for s in 0..lab.len() {
        if lab[s] != l || region[s] != usize::MAX {
            continue;
        }
        let id = sizes.len();
        region[s] = id;
        let mut stack = vec![s];
        let mut size = 0;
        while let Some(k) = stack.pop() {
            size += 1;
            for m in neighbours(k, nx, ny) {
                if lab[m] == l && region[m] == usize::MAX {
                    region[m] = id;
                    stack.push(m);
                }
            }
        }
        sizes.push(size);
    }
    let Some(id) = (0..sizes.len()).max_by_key(|&r| sizes[r]) else {
        return Vec::new();
    };
    let inside = |i: isize, j: isize| {
        i >= 0
            && j >= 0
            && (i as usize) < nx
            && (j as usize) < ny
            && region[j as usize * nx + i as usize] == id
    };
    let mut out: std::collections::BTreeMap<(isize, isize), Vec<(isize, isize)>> =
        std::collections::BTreeMap::new();
    for k in (0..lab.len()).filter(|&k| region[k] == id) {
        let (i, j) = ((k % nx) as isize, (k / nx) as isize);
        if !inside(i, j - 1) {
            out.entry((i, j)).or_default().push((1, 0));
        }
        if !inside(i + 1, j) {
            out.entry((i + 1, j)).or_default().push((0, 1));
        }
        if !inside(i, j + 1) {
            out.entry((i + 1, j + 1)).or_default().push((-1, 0));
        }
        if !inside(i - 1, j) {
            out.entry((i, j + 1)).or_default().push((0, -1));
        }
    }
    let mut rings: Vec<Vec<(usize, usize)>> = Vec::new();
    while let Some((&start, _)) = out
        .iter()
        .find(|(_, v)| v.len() == 1)
        .or_else(|| out.iter().find(|(_, v)| !v.is_empty()))
    {
        let mut v = start;
        let mut dir = out.get_mut(&v).unwrap().pop().unwrap();
        let mut corners = vec![v];
        loop {
            v = (v.0 + dir.0, v.1 + dir.1);
            if v == start {
                break;
            }
            let edges = out.get_mut(&v).unwrap();
            // Turning right at a corner two blocks share diagonally keeps them apart.
            let right = (dir.1, -dir.0);
            let pick = edges
                .iter()
                .position(|&d| d == right)
                .or_else(|| edges.iter().position(|&d| d == dir))
                .unwrap_or(0);
            let next = edges.swap_remove(pick);
            if next != dir {
                corners.push(v);
            }
            dir = next;
        }
        // The walk may start mid-edge.
        let m = corners.len();
        let corners: Vec<(usize, usize)> = (0..m)
            .filter(|&k| {
                let (a, b, c) = (corners[(k + m - 1) % m], corners[k], corners[(k + 1) % m]);
                (b.0 - a.0) * (c.1 - b.1) != (b.1 - a.1) * (c.0 - b.0)
            })
            .map(|k| (corners[k].0 as usize, corners[k].1 as usize))
            .collect();
        if corners.len() >= 4 {
            rings.push(corners);
        }
    }
    // Outer ring first: the largest in area.
    let area = |r: &Vec<(usize, usize)>| {
        let m = r.len();
        (0..m)
            .map(|k| {
                r[k].0 as f64 * r[(k + 1) % m].1 as f64 - r[(k + 1) % m].0 as f64 * r[k].1 as f64
            })
            .sum::<f64>()
            .abs()
    };
    rings.sort_by(|a, b| area(b).total_cmp(&area(a)));
    rings
}
