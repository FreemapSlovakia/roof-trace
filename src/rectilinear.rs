use crate::outline::{Grid, ring_area};

type P = [f64; 2];

/// Outline of a building whose walls follow the directions `theta` and `theta` + 90°: its points
/// are rasterized in a grid rotated to that direction, so walls fall on cell edges, and steps
/// shorter than `jog` are removed. Cells under `cover` (canopy over the roof) are filled within the
/// house's rows and columns. Returns the rings, outer first, open, in map units, and the mask area.
pub fn outline(points: &[P], cover: &[P], theta: f64, res: f64, jog: f64, min_area: f64) -> Option<(Vec<Vec<P>>, f64)> {
    let n = points.len() as f64;
    let c = points.iter().fold([0.0, 0.0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n]);
    let (sin, cos) = theta.sin_cos();
    let to_rot = |p: &P| {
        let (dx, dy) = (p[0] - c[0], p[1] - c[1]);
        [dx * cos + dy * sin, -dx * sin + dy * cos, 0.0]
    };
    let rot: Vec<[f64; 3]> = points.iter().map(to_rot).collect();
    let rot_cover: Vec<[f64; 2]> = cover.iter().map(|p| {
        let r = to_rot(p);
        [r[0], r[1]]
    }).collect();
    let mut q = [f64::MAX, f64::MAX, f64::MIN, f64::MIN];
    for p in &rot {
        q = [q[0].min(p[0]), q[1].min(p[1]), q[2].max(p[0]), q[3].max(p[1])];
    }
    // Margin keeps the outline off the grid edge.
    let q = [q[0] - 2.0, q[1] - 2.0, q[2] + 2.0, q[3] + 2.0];
    let mut g = Grid::new(q, res, &rot, &[]);
    g.complete(&rot_cover);
    let b = g.buildings(min_area).into_iter().max_by(|a, b| a.area.total_cmp(&b.area))?;
    let rings: Vec<Vec<P>> = b
        .rings
        .iter()
        .filter_map(|r| simplify(r, jog))
        .map(|r| r.into_iter().map(|[x, y]| [c[0] + x * cos - y * sin, c[1] + x * sin + y * cos]).collect())
        .collect();
    if rings.is_empty() || ring_area(&rings[0]) <= 0.0 {
        return None;
    }
    Some((rings, b.area))
}

/// Removes steps shorter than `jog` from a ring of alternating horizontal and vertical edges,
/// each time merging the two parallel edges either side of the shortest step.
fn simplify(corners: &[P], jog: f64) -> Option<Vec<P>> {
    // The tracer may start mid-edge; a point between collinear edges would break the alternation.
    let n = corners.len();
    let flat = |a: P, b: P| (a[1] - b[1]).abs() < 1e-9;
    let corners: Vec<P> = (0..n)
        .filter(|&i| flat(corners[(i + n - 1) % n], corners[i]) != flat(corners[i], corners[(i + 1) % n]))
        .map(|i| corners[i])
        .collect();
    let n = corners.len();
    if n < 4 || n % 2 == 1 {
        return None;
    }
    // Edge j runs from corner j to j + 1 at a fixed y (horizontal) or x.
    let mut horiz: Vec<bool> = (0..n).map(|j| (corners[j][1] - corners[(j + 1) % n][1]).abs() < 1e-9).collect();
    let mut level: Vec<f64> = (0..n).map(|j| if horiz[j] { corners[j][1] } else { corners[j][0] }).collect();
    let len = |level: &[f64], j: usize| {
        let m = level.len();
        (level[(j + 1) % m] - level[(j + m - 1) % m]).abs()
    };
    loop {
        let m = level.len();
        if m <= 4 {
            break;
        }
        let (j, l) = (0..m).map(|j| (j, len(&level, j))).min_by(|a, b| a.1.total_cmp(&b.1)).unwrap();
        if l >= jog {
            break;
        }
        let (a, b) = ((j + m - 1) % m, (j + 1) % m);
        let (la, lb) = (len(&level, a), len(&level, b));
        let merged = if la + lb > 0.0 { (level[a] * la + level[b] * lb) / (la + lb) } else { (level[a] + level[b]) / 2.0 };
        level[a] = merged;
        // Drop the step and the edge merged into `a`, higher index first.
        let (hi, lo) = if j > b { (j, b) } else { (b, j) };
        for k in [hi, lo] {
            level.remove(k);
            horiz.remove(k);
        }
    }
    let m = level.len();
    let out: Vec<P> = (0..m)
        .map(|j| {
            let next = level[(j + 1) % m];
            if horiz[j] { [next, level[j]] } else { [level[j], next] }
        })
        .collect();
    (ring_area(&out).abs() > 0.0).then_some(out)
}

/// Fraction of `points` outside the outer ring (holes ignored).
pub fn outside(points: &[P], ring: &[P]) -> f64 {
    let n = ring.len();
    let out = points
        .iter()
        .filter(|p| {
            let mut inside = false;
            for i in 0..n {
                let (a, b) = (ring[i], ring[(i + 1) % n]);
                if (a[1] > p[1]) != (b[1] > p[1]) && p[0] < a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]) {
                    inside = !inside;
                }
            }
            !inside
        })
        .count();
    out as f64 / points.len() as f64
}

/// Angles tried around each candidate direction: a degree or two off turns a long wall into steps.
const SEARCH_DEG: f64 = 3.0;
const SEARCH_STEP_DEG: f64 = 0.25;

/// Largest share of the building's points left outside an accepted outline.
const MAX_OUTSIDE: f64 = 0.05;

/// The right-angled outline with the fewest corners over directions near `thetas` that still fits
/// the points and its mask area; `None` when none fits, as for round or oblique buildings.
pub fn best(points: &[P], cover: &[P], thetas: &[f64], res: f64, jog: f64, min_area: f64, max_area_change: f64) -> Option<Vec<Vec<P>>> {
    let steps = (SEARCH_DEG / SEARCH_STEP_DEG) as i32;
    let corners = |r: &Vec<Vec<P>>| r.iter().map(Vec::len).sum::<usize>();
    let fit = |cover: &[P]| {
        thetas
            .iter()
            .flat_map(|&t| (-steps..=steps).map(move |k| t + (k as f64 * SEARCH_STEP_DEG).to_radians()))
            .filter_map(|t| outline(points, cover, t, res, jog, min_area))
            .filter(|(r, area)| (ring_area(&r[0]) / area - 1.0).abs() < max_area_change && outside(points, &r[0]) < MAX_OUTSIDE)
            .map(|(r, _)| r)
            .min_by_key(corners)
    };
    let plain = fit(&[]);
    if cover.is_empty() {
        return plain;
    }
    // Canopy beside a house is no evidence of roof; only completion that simplifies the outline is.
    match (plain, fit(cover)) {
        (Some(p), Some(c)) => Some(if corners(&c) < corners(&p) { c } else { p }),
        (p, c) => p.or(c),
    }
}
