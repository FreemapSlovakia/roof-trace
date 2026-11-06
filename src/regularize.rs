use std::{
    f64::consts::FRAC_PI_2,
    sync::atomic::{AtomicUsize, Ordering::Relaxed},
};

use crate::outline::ring_area;

type P = [f64; 2];

/// A boundary edge as an infinite line through `p`.
struct Line {
    p: P,
    dir: P,
    len: f64,
    /// Orientation ×2 + axis (0 or 1) when snapped to one of the building's orientations.
    axis: Option<u8>,
    /// Original corner where the edge ends.
    end: P,
}

/// Simplifies a closed ring and snaps its edges to the building's two main directions.
/// Returns `None` for rings that collapse.
pub fn ring(r: &[P], tol: f64, snap: f64, merge: f64, debug: bool) -> Option<Vec<P>> {
    let s = simplify_closed(r, tol);
    if debug {
        eprintln!("raw {} corners, simplified {}: {s:?}", r.len(), s.len());
    }
    if s.len() < 3 {
        return None;
    }
    let n = s.len();

    let edges: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let (a, b) = (s[i], s[(i + 1) % n]);
            (dist(a, b), (b[1] - a[1]).atan2(b[0] - a[0]))
        })
        .collect();
    let thetas = orientations(&edges);
    if debug {
        let deg: Vec<f64> = thetas.iter().map(|t| t.to_degrees()).collect();
        eprintln!("orientations {deg:?}");
        for (l, a) in &edges {
            eprintln!("  edge {l:.2} m at {:.1}°", a.to_degrees());
        }
    }

    let mut lines: Vec<Line> = (0..n)
        .map(|i| {
            let (a, b) = (s[i], s[(i + 1) % n]);
            let (len, ang) = edges[i];
            let best = thetas
                .iter()
                .enumerate()
                .map(|(t, &theta)| {
                    let k = ((ang - theta) / FRAC_PI_2).round();
                    (t, theta, k, ang - theta - k * FRAC_PI_2)
                })
                .min_by(|x, y| x.3.abs().total_cmp(&y.3.abs()));
            let (ang, axis) = match best {
                Some((t, theta, k, d)) if d.abs() < snap => {
                    (theta + k * FRAC_PI_2, Some(t as u8 * 2 + k.rem_euclid(2.0) as u8))
                }
                _ => (ang, None),
            };
            Line { p: mid(a, b), dir: [ang.cos(), ang.sin()], len, axis, end: b }
        })
        .collect();

    // Short edges that fit no main direction are raster noise at corners; their neighbours meet instead.
    loop {
        let m = lines.len();
        let Some(i) = (0..m)
            .filter(|&i| lines[i].axis.is_none() && lines[i].len < MIN_FREE_EDGE)
            .min_by(|&a, &b| lines[a].len.total_cmp(&lines[b].len))
        else {
            break;
        };
        if m <= 4 {
            break;
        }
        let prev = (i + m - 1) % m;
        lines[prev].end = lines[i].end;
        lines.remove(i);
    }

    // Merge parallel edges that are only slightly offset, whether adjacent or joined by a short jog:
    // the steps a raster outline leaves along a wall.
    loop {
        let m = lines.len();
        if m < 4 {
            return None;
        }
        let parallel = |a: &Line, b: &Line| a.axis.is_some() && a.axis == b.axis && offset(a, b.p).abs() < merge;
        let Some((i, jog)) = (0..m).find_map(|i| {
            let (a, b, c) = (&lines[i], &lines[(i + 1) % m], &lines[(i + 2) % m]);
            if parallel(a, b) {
                Some((i, false))
            } else if m > 4 && b.len < 2.0 * merge && parallel(a, c) {
                Some((i, true))
            } else {
                None
            }
        }) else {
            break;
        };
        if jog {
            lines.remove((i + 1) % m);
        }
        let m = lines.len();
        let i = if jog && (i + 1) % (m + 1) < i { i - 1 } else { i };
        let j = (i + 1) % m;
        let (a, b) = (&lines[i], &lines[j]);
        let w = a.len + b.len;
        let merged = Line {
            p: [(a.p[0] * a.len + b.p[0] * b.len) / w, (a.p[1] * a.len + b.p[1] * b.len) / w],
            dir: a.dir,
            len: w,
            axis: a.axis,
            end: b.end,
        };
        lines[i] = merged;
        lines.remove(j);
    }

    if debug {
        for l in &lines {
            eprintln!("  line {:.2} m axis {:?} dir {:.1}°", l.len, l.axis, l.dir[1].atan2(l.dir[0]).to_degrees());
        }
    }
    let m = lines.len();
    let mut out = Vec::with_capacity(m * 2);
    for i in 0..m {
        let (a, b) = (&lines[i], &lines[(i + 1) % m]);
        let cross = a.dir[0] * b.dir[1] - a.dir[1] * b.dir[0];
        let x = if cross.abs() > 0.2 { intersect(a, b) } else { None };
        match x {
            // Far-off intersections come from nearly parallel edges; a step between them is truer.
            Some(x) if dist(x, a.end) < 3.0 => out.push(x),
            _ => {
                out.push(project(a, a.end));
                out.push(project(b, a.end));
            }
        }
    }
    let out = clean(out);
    if out.len() < 3 {
        return None;
    }
    // Guard against a regularization that turned the ring inside out or lost much of it.
    let (a0, a1) = (ring_area(r), ring_area(&out));
    if a0.signum() != a1.signum() || (a1 / a0 - 1.0).abs() > 0.25 {
        FALLBACKS.fetch_add(1, Relaxed);
        return Some(s);
    }
    Some(out)
}

/// Unsnapped edges shorter than this are dropped (m).
const MIN_FREE_EDGE: f64 = 1.0;

/// A minor orientation needs this share of the outline's length; merged neighbouring houses often have two or three.
const MIN_ORIENTATION_SHARE: f64 = 0.15;

/// Dominant directions (radians, modulo 90°) of edges given as (length, angle).
fn orientations(edges: &[(f64, f64)]) -> Vec<f64> {
    let total: f64 = edges.iter().map(|e| e.0).sum();
    let mut hist = [0.0; 90];
    for &(len, ang) in edges {
        let b = (ang.to_degrees().rem_euclid(90.0) as usize).min(89);
        for d in -3..=3_isize {
            hist[(b as isize + d).rem_euclid(90) as usize] += len * (4 - d.abs()) as f64;
        }
    }
    let mut out = Vec::new();
    while out.len() < 3 {
        let (peak, &w) = hist.iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1)).unwrap();
        if w <= 0.0 || (!out.is_empty() && w < MIN_ORIENTATION_SHARE * total * 4.0) {
            break;
        }
        // Refine with the length-weighted mean of nearby edges, angles ×4 so 90° apart agree.
        let peak_r = (peak as f64 + 0.5).to_radians();
        let (mut sx, mut sy) = (0.0, 0.0);
        for &(len, ang) in edges {
            let d = ((ang - peak_r) / FRAC_PI_2 - ((ang - peak_r) / FRAC_PI_2).round()) * FRAC_PI_2;
            if d.abs() < 10f64.to_radians() {
                sx += len * (4.0 * ang).cos();
                sy += len * (4.0 * ang).sin();
            }
        }
        out.push(sy.atan2(sx) / 4.0);
        for d in -15..=15_isize {
            hist[(peak as isize + d).rem_euclid(90) as usize] = 0.0;
        }
    }
    out
}

/// Rings whose regularization was rejected and kept only simplified.
pub static FALLBACKS: AtomicUsize = AtomicUsize::new(0);

fn simplify_closed(r: &[P], tol: f64) -> Vec<P> {
    let far = (1..r.len()).max_by(|&a, &b| dist(r[0], r[a]).total_cmp(&dist(r[0], r[b]))).unwrap_or(0);
    let mut a = r[..=far].to_vec();
    let mut b = r[far..].to_vec();
    b.push(r[0]);
    a = dp(&a, tol);
    b = dp(&b, tol);
    a.pop();
    b.pop();
    a.extend(b);
    a
}

/// Douglas-Peucker on an open polyline.
fn dp(pts: &[P], tol: f64) -> Vec<P> {
    if pts.len() < 3 {
        return pts.to_vec();
    }
    let (a, b) = (pts[0], pts[pts.len() - 1]);
    let (k, d) = (1..pts.len() - 1)
        .map(|i| (i, seg_dist(pts[i], a, b)))
        .max_by(|x, y| x.1.total_cmp(&y.1))
        .unwrap();
    if d <= tol {
        return vec![a, b];
    }
    let mut left = dp(&pts[..=k], tol);
    left.pop();
    left.extend(dp(&pts[k..], tol));
    left
}

/// Drops repeated and collinear corners.
fn clean(mut v: Vec<P>) -> Vec<P> {
    loop {
        let n = v.len();
        if n < 3 {
            return v;
        }
        let Some(i) = (0..n).find(|&i| {
            let (a, b, c) = (v[(i + n - 1) % n], v[i], v[(i + 1) % n]);
            dist(a, b) < 0.05 || seg_dist(b, a, c) < 0.05
        }) else {
            return v;
        };
        v.remove(i);
    }
}

fn intersect(a: &Line, b: &Line) -> Option<P> {
    let den = a.dir[0] * b.dir[1] - a.dir[1] * b.dir[0];
    if den.abs() < 1e-9 {
        return None;
    }
    let t = ((b.p[0] - a.p[0]) * b.dir[1] - (b.p[1] - a.p[1]) * b.dir[0]) / den;
    Some([a.p[0] + t * a.dir[0], a.p[1] + t * a.dir[1]])
}

fn project(l: &Line, q: P) -> P {
    let t = (q[0] - l.p[0]) * l.dir[0] + (q[1] - l.p[1]) * l.dir[1];
    [l.p[0] + t * l.dir[0], l.p[1] + t * l.dir[1]]
}

/// Signed distance of `q` from the line.
fn offset(l: &Line, q: P) -> f64 {
    (q[0] - l.p[0]) * -l.dir[1] + (q[1] - l.p[1]) * l.dir[0]
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

fn mid(a: P, b: P) -> P {
    [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0]
}

/// Dominant wall directions of a ring (radians, modulo 90°), strongest first.
/// Direction of the long side of the smallest rectangle around `r`: the walls' direction even
/// where the traced edges are too ragged to vote for it.
pub fn box_orientation(r: &[P]) -> Option<f64> {
    use geo::MinimumRotatedRect;
    let ring: Vec<geo::Coord<f64>> = r.iter().map(|&[x, y]| geo::Coord { x, y }).collect();
    let b = geo::Polygon::new(geo::LineString(ring), vec![]).minimum_rotated_rect()?;
    let c = &b.exterior().0;
    let side = |i: usize| (c[i + 1].x - c[i].x, c[i + 1].y - c[i].y);
    let (a, b) = (side(0), side(1));
    let (dx, dy) = if a.0.hypot(a.1) >= b.0.hypot(b.1) { a } else { b };
    Some(dy.atan2(dx))
}

pub fn wall_orientations(r: &[P], tol: f64) -> Vec<f64> {
    let s = simplify_closed(r, tol);
    let n = s.len();
    let edges: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let (a, b) = (s[i], s[(i + 1) % n]);
            (dist(a, b), (b[1] - a[1]).atan2(b[0] - a[0]))
        })
        .collect();
    orientations(&edges)
}
