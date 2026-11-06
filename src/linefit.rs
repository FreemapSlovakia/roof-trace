use std::f64::consts::FRAC_PI_2;

type P = [f64; 2];

/// A wall fitted to a run of boundary points.
struct Line {
    pts: Vec<P>,
    c: P,
    dir: P,
    /// Snapped to the building's main axes.
    ortho: bool,
}

impl Line {
    fn fit(pts: Vec<P>) -> Line {
        let n = pts.len() as f64;
        let c = pts.iter().fold([0.0, 0.0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n]);
        let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
        for p in &pts {
            let (dx, dy) = (p[0] - c[0], p[1] - c[1]);
            sxx += dx * dx;
            sxy += dx * dy;
            syy += dy * dy;
        }
        let ang = 0.5 * (2.0 * sxy).atan2(sxx - syy);
        Line { pts, c, dir: [ang.cos(), ang.sin()], ortho: false }
    }

    fn offset(&self, q: P) -> f64 {
        (q[0] - self.c[0]) * -self.dir[1] + (q[1] - self.c[1]) * self.dir[0]
    }

    fn along(&self, q: P) -> f64 {
        (q[0] - self.c[0]) * self.dir[0] + (q[1] - self.c[1]) * self.dir[1]
    }

    fn len(&self) -> f64 {
        let (lo, hi) = self.pts.iter().fold((f64::MAX, f64::MIN), |(lo, hi), &p| {
            let t = self.along(p);
            (lo.min(t), hi.max(t))
        });
        hi - lo
    }

    fn rms(&self) -> f64 {
        (self.pts.iter().map(|&p| self.offset(p).powi(2)).sum::<f64>() / self.pts.len() as f64).sqrt()
    }

    fn angle(&self) -> f64 {
        self.dir[1].atan2(self.dir[0])
    }

    fn project(&self, q: P) -> P {
        let t = self.along(q);
        [self.c[0] + t * self.dir[0], self.c[1] + t * self.dir[1]]
    }
}

/// Runs of boundary points straighter than this are one wall (m).
const SPLIT_TOL: f64 = 0.5;
/// Neighbouring walls closer in direction than this are tried as one (degrees).
const MERGE_ANGLE: f64 = 12.0;
/// ... and merged when the points still fit one line this well (RMS, m).
const MERGE_RMS: f64 = 0.3;
/// Walls shorter than this are raster noise at corners (m).
const MIN_WALL: f64 = 0.8;
/// Walls within this angle of the main axes are made exactly parallel or perpendicular (degrees).
const ORTHO_TOL: f64 = 4.0;
/// Walls at least this long decide whether a house is square (m).
const MAIN_WALL: f64 = 3.0;
/// Walls off the main axes shorter than this are dropped (m).
const MIN_SKEWED_WALL: f64 = 2.0;
/// Parallel neighbouring walls offset less than this are one wall (m).
const PARALLEL_MERGE: f64 = 0.8;

/// Fits straight walls to a traced ring (cell-edge corners, `res` apart at most) and snaps those
/// close to the building's main axes; `hidden` tells boundary under tree crowns. Returns the ring and whether every wall was snapped.
pub fn ring(r: &[P], res: f64, hidden: &dyn Fn(&P) -> bool) -> Option<(Vec<P>, bool)> {
    // Midpoints of the unit cell edges: the staircase sampled evenly.
    let n = r.len();
    let mut dense = Vec::new();
    for i in 0..n {
        let (a, b) = (r[i], r[(i + 1) % n]);
        let k = ((dist(a, b) / res).round() as usize).max(1);
        for s in 0..k {
            let t = (s as f64 + 0.5) / k as f64;
            dense.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
        }
    }
    let m = dense.len();
    if m < 8 {
        return None;
    }

    // Break points of a Douglas-Peucker pass split the boundary into straight runs.
    let far = (1..m).max_by(|&a, &b| dist(dense[0], dense[a]).total_cmp(&dist(dense[0], dense[b]))).unwrap();
    let mut breaks = vec![0];
    dp(&dense, 0, far, &mut breaks);
    breaks.push(far);
    let ext: Vec<P> = dense.iter().chain(std::iter::once(&dense[0])).copied().collect();
    dp(&ext, far, m, &mut breaks);
    breaks.sort_unstable();
    breaks.dedup();
    let mut lines: Vec<Line> = (0..breaks.len())
        .map(|i| {
            let (a, b) = (breaks[i], if i + 1 < breaks.len() { breaks[i + 1] } else { m });
            Line::fit((a..b).map(|k| dense[k % m]).collect())
        })
        .filter(|l| l.pts.len() >= 2)
        .collect();

    // Merge nearly collinear neighbours while their points still fit one line.
    loop {
        let k = lines.len();
        if k < 4 {
            return None;
        }
        let Some((i, merged)) = (0..k).find_map(|i| {
            let (a, b) = (&lines[i], &lines[(i + 1) % k]);
            if angle_diff(a.angle(), b.angle()) > MERGE_ANGLE.to_radians() {
                return None;
            }
            let l = Line::fit(a.pts.iter().chain(&b.pts).copied().collect());
            (l.rms() < MERGE_RMS).then_some((i, l))
        }) else {
            break;
        };
        lines[i] = merged;
        lines.remove((i + 1) % k);
    }

    // Short walls are corner noise; their neighbours meet instead.
    loop {
        let k = lines.len();
        if k <= 4 {
            break;
        }
        let Some(i) = (0..k).filter(|&i| lines[i].len() < MIN_WALL).min_by(|&a, &b| lines[a].len().total_cmp(&lines[b].len()))
        else {
            break;
        };
        lines.remove(i);
    }

    // Main direction: length-weighted mean of wall angles ×4, so walls 90° apart agree.
    let (mut sx, mut sy) = (0.0, 0.0);
    for l in &lines {
        let len = l.len();
        sx += len * (4.0 * l.angle()).cos();
        sy += len * (4.0 * l.angle()).sin();
    }
    let theta = sy.atan2(sx) / 4.0;
    for l in lines.iter_mut() {
        let k = ((l.angle() - theta) / FRAC_PI_2).round();
        let d = l.angle() - theta - k * FRAC_PI_2;
        if d.abs() < ORTHO_TOL.to_radians() {
            let a = theta + k * FRAC_PI_2;
            l.dir = [a.cos(), a.sin()];
            l.ortho = true;
        }
    }

    // A short wall off the main axes is a cut corner from the raster, and one along a tree crown
    // follows the tree: neither is a skewed wall.
    let shaded = |l: &Line| l.pts.iter().filter(|p| hidden(p)).count() * 2 > l.pts.len();
    loop {
        let k = lines.len();
        if k <= 4 {
            break;
        }
        let Some(i) = (0..k).filter(|&i| !lines[i].ortho && (lines[i].len() < MIN_SKEWED_WALL || shaded(&lines[i]))).min_by(|&a, &b| lines[a].len().total_cmp(&lines[b].len()))
        else {
            break;
        };
        lines.remove(i);
    }

    // Parallel neighbours only slightly offset are one wall with a raster jog.
    loop {
        let k = lines.len();
        if k <= 4 {
            break;
        }
        let Some(i) = (0..k).find(|&i| {
            let (a, b) = (&lines[i], &lines[(i + 1) % k]);
            angle_diff(a.angle(), b.angle()) < 1e-6 && a.offset(b.c).abs() < PARALLEL_MERGE
        }) else {
            break;
        };
        let j = (i + 1) % k;
        let (a, b) = (&lines[i], &lines[j]);
        let pts: Vec<P> = a.pts.iter().chain(&b.pts).copied().collect();
        let n = pts.len() as f64;
        let c = pts.iter().fold([0.0, 0.0], |s, p| [s[0] + p[0] / n, s[1] + p[1] / n]);
        let merged = Line { pts, c, dir: a.dir, ortho: a.ortho && b.ortho };
        lines[i] = merged;
        lines.remove(j);
    }

    let k = lines.len();
    let mut out = Vec::with_capacity(k * 2);
    for i in 0..k {
        let (a, b) = (&lines[i], &lines[(i + 1) % k]);
        let joint = *a.pts.last().unwrap();
        let cross = a.dir[0] * b.dir[1] - a.dir[1] * b.dir[0];
        match intersect(a, b) {
            // Nearly parallel walls meet far away; a short step between them is truer.
            Some(x) if cross.abs() > 0.25 && dist(x, joint) < 3.0 => out.push(x),
            _ => {
                out.push(a.project(joint));
                out.push(b.project(joint));
            }
        }
    }
    // Only real walls tell whether the house is square; short ones are often cut corners.
    let ortho = lines.iter().filter(|l| l.len() >= MAIN_WALL).all(|l| l.ortho);
    (out.len() >= 3).then_some((out, ortho))
}

fn dp(pts: &[P], lo: usize, hi: usize, out: &mut Vec<usize>) {
    if hi <= lo + 1 {
        return;
    }
    let (a, b) = (pts[lo], pts[hi.min(pts.len() - 1)]);
    let (k, d) = (lo + 1..hi).map(|i| (i, seg_dist(pts[i], a, b))).max_by(|x, y| x.1.total_cmp(&y.1)).unwrap();
    if d > SPLIT_TOL {
        dp(pts, lo, k, out);
        out.push(k);
        dp(pts, k, hi, out);
    }
}

fn intersect(a: &Line, b: &Line) -> Option<P> {
    let den = a.dir[0] * b.dir[1] - a.dir[1] * b.dir[0];
    if den.abs() < 1e-9 {
        return None;
    }
    let t = ((b.c[0] - a.c[0]) * b.dir[1] - (b.c[1] - a.c[1]) * b.dir[0]) / den;
    Some([a.c[0] + t * a.dir[0], a.c[1] + t * a.dir[1]])
}

/// Difference of two line directions, modulo 180°.
fn angle_diff(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(std::f64::consts::PI);
    d.min(std::f64::consts::PI - d)
}

fn seg_dist(p: P, a: P, b: P) -> f64 {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let l2 = dx * dx + dy * dy;
    let t = if l2 == 0.0 { 0.0 } else { (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2).clamp(0.0, 1.0) };
    dist(p, [a[0] + t * dx, a[1] + t * dy])
}

fn dist(a: P, b: P) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1])
}
