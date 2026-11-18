//! Houses built side by side under parallel roofs meet in a valley running straight along their
//! walls; one house's valleys run diagonally, from its wings.

use crate::store::Pt;

const RES: f64 = 0.5;
/// How far to each side the roof must rise from the valley (m), and by how much.
const REACH: f64 = 1.5;
const RISE: f64 = 0.3;
/// Share of the line's cells that must lie in the valley, and its least length (m).
const SHARE: f64 = 0.6;
const MIN_LENGTH: f64 = 6.0;
/// Narrowest house on either side (m).
const MIN_SIDE: f64 = 3.0;

/// The points split along the deepest straight valley along `theta` or across it, if any.
pub fn cut(points: &[Pt], theta: f64) -> Option<(Vec<Pt>, Vec<Pt>)> {
    [theta, theta + std::f64::consts::FRAC_PI_2]
        .into_iter()
        .filter_map(|t| line(points, t).map(|(score, at)| (score, t, at)))
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .map(|(_, t, at)| {
            let (s, c) = t.sin_cos();
            points.iter().partition(|p| -p[0] * s + p[1] * c < at)
        })
}

/// The best valley line parallel to direction `t`: its share of valley cells and its offset.
fn line(points: &[Pt], t: f64) -> Option<(f64, f64)> {
    let (s, c) = t.sin_cos();
    // u along the line, v across it.
    let uv: Vec<[f64; 3]> = points
        .iter()
        .map(|p| [p[0] * c + p[1] * s, -p[0] * s + p[1] * c, p[2]])
        .collect();
    let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
    for p in &uv {
        lo = [lo[0].min(p[0]), lo[1].min(p[1])];
        hi = [hi[0].max(p[0]), hi[1].max(p[1])];
    }
    let (nu, nv) = (
        ((hi[0] - lo[0]) / RES) as usize + 1,
        ((hi[1] - lo[1]) / RES) as usize + 1,
    );
    if (nu as f64) * RES < MIN_LENGTH || (nv as f64) * RES < 2.0 * MIN_SIDE || nu * nv > 4_000_000 {
        return None;
    }
    let mut z = vec![f64::NAN; nu * nv];
    for p in &uv {
        let k = ((p[1] - lo[1]) / RES) as usize * nu + ((p[0] - lo[0]) / RES) as usize;
        if !(z[k] >= p[2]) {
            z[k] = p[2];
        }
    }
    let reach = (REACH / RES).round() as usize;
    let side = (MIN_SIDE / RES).round() as usize;
    let mut best: Option<(f64, f64)> = None;
    for j in side.max(reach)..nv.saturating_sub(side.max(reach)) {
        let (mut on, mut valley) = (0usize, 0usize);
        for i in 0..nu {
            let here = z[j * nu + i];
            if here.is_nan() {
                continue;
            }
            on += 1;
            // The lowest of a few cells across, so the valley may wander a cell.
            let low = (j - 1..=j + 1)
                .map(|jj| z[jj * nu + i])
                .filter(|v| !v.is_nan())
                .fold(here, f64::min);
            let (a, b) = (z[(j - reach) * nu + i], z[(j + reach) * nu + i]);
            valley += (a > low + RISE && b > low + RISE) as usize;
        }
        let share = valley as f64 / on.max(1) as f64;
        if (valley as f64) * RES >= MIN_LENGTH
            && share >= SHARE
            && best.is_none_or(|(s, _)| share > s)
        {
            best = Some((share, lo[1] + (j as f64 + 0.5) * RES));
        }
    }
    best
}
