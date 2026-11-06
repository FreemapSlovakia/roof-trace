use std::collections::BTreeMap;

use crate::split;

type P3 = [f64; 3];

/// Grid for finding roof faces (m), as for splitting.
const RES: f64 = 0.5;
/// Slopes smaller than this share of the sloped roof are dormers or noise, not its shape ...
const MIN_SLOPE_SHARE: f64 = 0.08;
/// ... unless together they are this much of it.
const MAX_MINOR: f64 = 0.25;
/// Faces looking within this many degrees of each other are one slope.
const SAME_ASPECT: f64 = 30.0;
/// Faces of one slope differing more than this in pitch (rise per m) make a mansard or gambrel.
const PITCH_TOL: f64 = 0.3;
/// Share of flat roof that makes it a flat roof ...
const FLAT_ROOF: f64 = 0.85;
/// ... and above which, short of that, the roof is mixed.
const MAX_FLAT_PART: f64 = 0.25;
/// Slope (rise per m) under which a face is flat: about 6°.
const FLAT_SLOPE: f64 = 0.1;
/// Share of the edges between slopes that must be ridges or hips for the roof to be simple.
const CONVEX_SHARE: f64 = 0.65;
/// Faces whose slopes point apart at least this much (degrees) are opposite.
const OPPOSITE: f64 = 160.0;
/// Pairs of opposite faces this close to right angles (degrees) form a hipped or pyramidal roof.
const SQUARE_TOL: f64 = 20.0;
/// End faces with at least this share of the roof make a pyramid rather than a hipped roof.
const PYRAMID_SHARE: f64 = 0.4;
/// Share of the cells with points that faces must cover for the shape to be known.
const MIN_COVER: f64 = 0.7;
/// Points this close below the top (m) are the top edge ...
const TOP_BAND: f64 = 0.3;
/// ... which on a skillion lies within this share of the footprint from its uphill side.
const TOP_EDGE: f64 = 0.2;
/// Footprints at least this much longer than wide have a ridge direction worth tagging.
const ELONGATED: f64 = 1.2;

/// What OSM can be told about one house's roof.
pub struct Roof {
    /// `roof:shape`, when the roof is one of the simple shapes.
    pub shape: Option<&'static str>,
    /// The ridge runs across the footprint's long side (`roof:orientation=across`).
    pub across: bool,
    /// Ridge height without chimneys and antennas (m above sea level).
    pub top: f64,
    /// Height of the eaves (m above sea level).
    pub eave: f64,
}

/// Shape and heights of one house's roof from its points. A roof with valleys (wings, cross gables)
/// or faces at odd angles gets no shape: none is better than a wrong one.
pub fn roof(points: &[P3], debug: bool) -> Roof {
    let mut z: Vec<f64> = points.iter().map(|p| p[2]).collect();
    z.sort_by(f64::total_cmp);
    let pct = |q: f64| z[((z.len() - 1) as f64 * q) as usize];
    let (top, eave) = (pct(0.99), pct(0.03));
    let (shape, across) = shape(points, debug).map_or((None, false), |(s, a)| (Some(s), a));
    Roof { shape, across, top, eave }
}

fn shape(points: &[P3], debug: bool) -> Option<(&'static str, bool)> {
    let f = split::faces(points, RES);
    let pos = |k: usize| [((k % f.w) as f64 + 0.5) * RES, ((k / f.w) as f64 + 0.5) * RES];
    let mut cells: BTreeMap<usize, usize> = BTreeMap::new();
    for &l in f.label.iter().filter(|&&l| l != 0) {
        *cells.entry(l).or_default() += 1;
    }
    let total: usize = cells.values().sum();
    if total == 0 {
        return None;
    }
    // Where the faces leave much of the roof out, as under a tree, its shape is unknown.
    let seen: std::collections::BTreeSet<(i64, i64)> =
        points.iter().map(|p| ((p[0] / RES).floor() as i64, (p[1] / RES).floor() as i64)).collect();
    if (total as f64) < MIN_COVER * seen.len() as f64 {
        return None;
    }
    let slope = |l: usize| f.planes[l].a.hypot(f.planes[l].b);
    // The way a face looks: down its slope.
    let aspect = |l: usize| (-f.planes[l].b).atan2(-f.planes[l].a).to_degrees();
    let apart = |a: f64, b: f64| {
        let d = (a - b).rem_euclid(360.0);
        d.min(360.0 - d)
    };
    let flat: usize = cells.iter().filter(|&(&l, _)| slope(l) < FLAT_SLOPE).map(|(_, &n)| n).sum();
    if flat as f64 >= FLAT_ROOF * total as f64 {
        return Some(("flat", false));
    }
    if flat as f64 > MAX_FLAT_PART * total as f64 {
        return None;
    }

    // A roof slope often comes as several faces; those looking the same way are one slope.
    let mut faces: Vec<usize> = cells.keys().copied().filter(|&l| slope(l) >= FLAT_SLOPE).collect();
    faces.sort_by_key(|l| std::cmp::Reverse(cells[l]));
    let sloped: usize = faces.iter().map(|l| cells[l]).sum();
    // Per slope: its faces, cells, and aspect of its largest face.
    let mut slopes: Vec<(Vec<usize>, usize, f64)> = Vec::new();
    for &l in &faces {
        match slopes.iter_mut().find(|s| apart(s.2, aspect(l)) < SAME_ASPECT) {
            Some(s) => {
                s.0.push(l);
                s.1 += cells[&l];
            }
            None => slopes.push((vec![l], cells[&l], aspect(l))),
        }
    }
    let (main, minor): (Vec<_>, Vec<_>) = slopes.into_iter().partition(|s| s.1 as f64 >= MIN_SLOPE_SHARE * sloped as f64);
    let minor_cells: usize = minor.iter().map(|s| s.1).sum();
    let slope_of = |l: usize| main.iter().position(|s| s.0.contains(&l));

    // Between different slopes, a simple roof has ridges and hips only; valleys mean wings.
    let (mut convex, mut edges) = (0, 0);
    for k in 0..f.label.len() {
        let Some(sr) = slope_of(f.label[k]) else { continue };
        for n in split::neighbours(k, f.w, f.ht) {
            let Some(sn) = slope_of(f.label[n]) else { continue };
            if sn <= sr {
                continue;
            }
            let m = split::meet(&f.planes[f.label[k]], &f.planes[f.label[n]], pos(k), pos(n), RES);
            convex += m.ridge as usize;
            edges += 1;
        }
    }
    // Faces of one slope at clearly different pitches: a mansard or gambrel roof.
    let broken = main.iter().any(|s| {
        let big: Vec<f64> = s.0.iter().filter(|l| cells[l] as f64 >= MIN_SLOPE_SHARE * sloped as f64 / 2.0).map(|&l| slope(l)).collect();
        big.iter().any(|a| big.iter().any(|b| (a - b).abs() > PITCH_TOL))
    });
    if debug {
        for s in &main {
            eprintln!("roof slope: aspect {:.0}, {} of {sloped} cells, faces {:?}", s.2, s.1, s.0);
        }
        eprintln!("roof: {minor_cells} cells in minor slopes, {convex} of {edges} edges convex, broken pitch {broken}");
    }
    if minor_cells as f64 > MAX_MINOR * sloped as f64 || broken || (edges > 0 && (convex as f64) < CONVEX_SHARE * edges as f64) {
        return None;
    }

    let big = main.iter().max_by_key(|s| s.1)?;
    // The ridge lies across the slope of the largest face; is it across the footprint too?
    let across = across(points, big.2 + 90.0);
    match main.len() {
        1 if top_at_edge(points, big.2) => Some(("skillion", false)),
        2 if apart(main[0].2, main[1].2) >= OPPOSITE => Some(("gabled", across)),
        4 => {
            // Two pairs of opposite slopes at right angles: the larger pair carries the ridge.
            let (pair, ends): (Vec<_>, Vec<_>) = main.iter().partition(|s| apart(s.2, big.2) >= OPPOSITE || s.2 == big.2);
            if pair.len() != 2 || apart(ends[0].2, ends[1].2) < OPPOSITE || (apart(pair[0].2, ends[0].2) - 90.0).abs() > SQUARE_TOL {
                return None;
            }
            let share = (ends[0].1 + ends[1].1) as f64 / sloped as f64;
            Some(if share >= PYRAMID_SHARE { ("pyramidal", false) } else { ("hipped", across) })
        }
        _ => None,
    }
}

/// Whether a ridge in direction `ridge` (degrees) runs across an elongated footprint.
fn across(points: &[P3], ridge: f64) -> bool {
    let n = points.len() as f64;
    let c = points.iter().fold([0.0, 0.0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n]);
    let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
    for p in points {
        let (dx, dy) = (p[0] - c[0], p[1] - c[1]);
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    // Principal axes of the footprint: the long one's direction and the two spreads.
    let long = 0.5 * (2.0 * sxy).atan2(sxx - syy);
    let d = ((sxx - syy).powi(2) / 4.0 + sxy * sxy).sqrt();
    let (l1, l2) = ((sxx + syy) / 2.0 + d, (sxx + syy) / 2.0 - d);
    if l2 <= 0.0 || (l1 / l2).sqrt() < ELONGATED {
        return false;
    }
    let off = (ridge.to_radians() - long).rem_euclid(std::f64::consts::PI);
    off.min(std::f64::consts::PI - off) > std::f64::consts::FRAC_PI_4
}

/// Whether the highest points lie along the uphill edge of the footprint, as on a skillion roof,
/// rather than in its middle, as on a gable roof whose other side was not seen.
fn top_at_edge(points: &[P3], aspect: f64) -> bool {
    let (s, c) = aspect.to_radians().sin_cos();
    let u = |p: &P3| p[0] * c + p[1] * s;
    let (lo, hi) = points.iter().fold((f64::MAX, f64::MIN), |(lo, hi), p| (lo.min(u(p)), hi.max(u(p))));
    let mut z: Vec<f64> = points.iter().map(|p| p[2]).collect();
    z.sort_by(f64::total_cmp);
    let top = z[(z.len() - 1) * 99 / 100];
    let ridge: Vec<f64> = points.iter().filter(|p| p[2] >= top - TOP_BAND).map(u).collect();
    let mean = ridge.iter().sum::<f64>() / ridge.len() as f64;
    hi > lo && (mean - lo) / (hi - lo) < TOP_EDGE
}
