use geo::{Area, BooleanOps, BoundingRect, Centroid, Coord, LineString, MultiPolygon, Polygon, Rect};

/// Overlaps larger than this share of the smaller house are not two neighbours meeting.
const MAX_SHARE: f64 = 0.3;

/// Neighbouring houses fitted apart can overlap along their common wall. Each overlap is cut by
/// the line along its length, and each house keeps its own side, so both end on one wall.
pub fn resolve(polys: &mut [Polygon<f64>]) -> usize {
    let rects: Vec<Option<Rect<f64>>> = polys.iter().map(|p| p.bounding_rect()).collect();
    let mut fixed = 0;
    for i in 0..polys.len() {
        for j in i + 1..polys.len() {
            let (Some(a), Some(b)) = (rects[i], rects[j]) else { continue };
            if a.max().x < b.min().x || b.max().x < a.min().x || a.max().y < b.min().y || b.max().y < a.min().y {
                continue;
            }
            let inter = polys[i].intersection(&polys[j]);
            let area = inter.unsigned_area();
            if area < 0.05 || area > MAX_SHARE * polys[i].unsigned_area().min(polys[j].unsigned_area()) {
                continue;
            }
            let pts: Vec<Coord<f64>> = inter.iter().flat_map(|p| p.exterior().coords().copied()).collect();
            let (c, d) = principal(&pts);
            let n = Coord { x: -d.y, y: d.x };
            let Some(ci) = polys[i].centroid() else { continue };
            let side = ((ci.x() - c.x) * n.x + (ci.y() - c.y) * n.y).signum();
            if side == 0.0 {
                continue;
            }
            let cut_i = inter.intersection(&half_plane(c, d, n, -side));
            let cut_j = inter.intersection(&half_plane(c, d, n, side));
            polys[i] = largest(polys[i].difference(&cut_i)).unwrap_or_else(|| polys[i].clone());
            polys[j] = largest(polys[j].difference(&cut_j)).unwrap_or_else(|| polys[j].clone());
            fixed += 1;
        }
    }
    fixed
}

/// Centroid and unit direction of the longest spread of `pts`.
fn principal(pts: &[Coord<f64>]) -> (Coord<f64>, Coord<f64>) {
    let n = pts.len() as f64;
    let c = pts.iter().fold(Coord { x: 0.0, y: 0.0 }, |a, p| Coord { x: a.x + p.x / n, y: a.y + p.y / n });
    let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
    for p in pts {
        let (dx, dy) = (p.x - c.x, p.y - c.y);
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    let a = 0.5 * (2.0 * sxy).atan2(sxx - syy);
    (c, Coord { x: a.cos(), y: a.sin() })
}

/// The side `s` (±1 along `n`) of the line through `c` along `d`, as a polygon large enough for a house.
fn half_plane(c: Coord<f64>, d: Coord<f64>, n: Coord<f64>, s: f64) -> MultiPolygon<f64> {
    const L: f64 = 200.0;
    let p = |u: f64, v: f64| Coord { x: c.x + d.x * u + n.x * s * v, y: c.y + d.y * u + n.y * s * v };
    MultiPolygon(vec![Polygon::new(LineString(vec![p(-L, 0.0), p(L, 0.0), p(L, L), p(-L, L), p(-L, 0.0)]), vec![])])
}

/// The largest piece, without the sliver holes a cut leaves where the two outlines nearly met.
fn largest(m: MultiPolygon<f64>) -> Option<Polygon<f64>> {
    let p = m.into_iter().max_by(|a, b| a.unsigned_area().total_cmp(&b.unsigned_area()))?;
    let (ext, holes) = p.into_inner();
    let holes = holes.into_iter().filter(|h| Polygon::new(h.clone(), vec![]).unsigned_area() >= 1.0).collect();
    Some(Polygon::new(ext, holes))
}

/// The polygon without self-intersections, such as a loop left where snapping folded a short wall
/// back on itself.
pub fn clean(p: &Polygon<f64>) -> Polygon<f64> {
    largest(MultiPolygon(vec![p.clone()]).union(&MultiPolygon(vec![]))).unwrap_or_else(|| p.clone())
}

/// Corners this close (m) to the line through their neighbours are not corners.
const STRAIGHT: f64 = 0.05;
/// Edges shorter than this (m) are left over from moving corners onto walls, not walls.
const MIN_EDGE: f64 = 0.25;

/// Drops repeated, straight-through and spike vertices that cutting and cleaning leave behind.
pub fn tidy(p: &Polygon<f64>) -> Polygon<f64> {
    let ring = |l: &LineString<f64>| {
        let mut v: Vec<Coord<f64>> = l.0.clone();
        v.pop();
        loop {
            let n = v.len();
            if n <= 3 {
                break;
            }
            let Some(i) = (0..n).find(|&i| {
                let (a, b, c) = (v[(i + n - 1) % n], v[i], v[(i + 1) % n]);
                let (dx, dy) = (c.x - a.x, c.y - a.y);
                let len = dx.hypot(dy);
                let off = if len < STRAIGHT { 0.0 } else { ((b.x - a.x) * dy - (b.y - a.y) * dx).abs() / len };
                // Straight through, out and straight back (a spike without width), or barely past the last corner.
                off < STRAIGHT || (b.x - a.x).hypot(b.y - a.y) < MIN_EDGE
            }) else {
                break;
            };
            v.remove(i);
        }
        v.push(v[0]);
        LineString(v)
    };
    let t = Polygon::new(ring(p.exterior()), p.interiors().iter().map(ring).collect());
    // Cutting a corner where the outline is narrow can make it cross itself.
    if geo::Validation::is_valid(&t) { t } else { p.clone() }
}
