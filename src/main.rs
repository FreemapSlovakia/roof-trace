mod linefit;
mod outline;
mod overlap;
mod rectilinear;
mod regularize;
mod roof;
mod split;
mod store;

use std::{
    fs::File,
    path::PathBuf,
    sync::{
        Mutex,
        atomic::{AtomicUsize, Ordering},
    },
};

use anyhow::{Result, bail};
use clap::Parser;
use rayon::prelude::*;
use rusqlite::Connection;
use serde_json::json;

use store::Cloud;

/// Traces regularized roof outlines from building points (class 6) of lidar flight lines.
///
/// With `--index` it reads the flight lines over the area directly; with `--bin` it sorts them into
/// a tile store instead, which `--store` then traces tile by tile.
#[derive(Parser)]
struct Args {
    /// Flight-line index (`laz_index` table of `laz2geotiff index`).
    #[arg(long)]
    index: Option<PathBuf>,

    /// Tile store to trace, made with `--bin`.
    #[arg(long, conflicts_with = "index")]
    store: Option<PathBuf>,

    /// Sort the flight lines' points around the area into this tile store; resumes where it stopped.
    #[arg(long, requires = "index")]
    bin: Option<PathBuf>,

    /// Only flight lines whose path matches this SQL `LIKE` pattern, such as `%LOT26%`.
    #[arg(long, default_value = "%")]
    files: String,

    /// Area as `xmin,ymin,xmax,ymax` in the data's CRS.
    #[arg(long, value_delimiter = ',', allow_hyphen_values = true)]
    bbox: Vec<f64>,

    /// Output GeoJSON.
    #[arg(long, required_unless_present = "bin")]
    out: Option<PathBuf>,

    /// Mask cell size (m).
    #[arg(long, default_value_t = 0.25)]
    res: f64,

    /// Smallest building kept (m²).
    #[arg(long, default_value_t = 10.0)]
    min_area: f64,

    /// Douglas-Peucker tolerance before regularization (m).
    #[arg(long, default_value_t = 0.4)]
    simplify: f64,

    /// Edges within this angle of the building's main axes are snapped to them (degrees).
    #[arg(long, default_value_t = 20.0)]
    snap_angle: f64,

    /// Parallel edges closer than this, adjacent or across a short jog, are merged into one (m).
    #[arg(long, default_value_t = 0.8)]
    merge_offset: f64,

    /// Keep attached houses together.
    #[arg(long)]
    no_split: bool,

    /// Steps shorter than this are removed from right-angled outlines (m).
    #[arg(long, default_value_t = 1.0)]
    jog: f64,

    /// Also write the highest building point per cell as an ESRI ASCII grid, for inspection.
    #[arg(long)]
    dsm: Option<PathBuf>,

    /// Check the result against expected cases: lines `same|split x1 y1 x2 y2`, `#` comments.
    #[arg(long)]
    check: Option<PathBuf>,

    /// EPSG code written to the output: S-JTSK [JTSK03] / Krovak East North, the lidar's CRS.
    #[arg(long, default_value_t = 8353)]
    epsg: u32,
}

/// Smallest house split off a blob (m²).
const MIN_PART: f64 = 30.0;

/// Corners this close to a shared wall are moved onto it (m).
const WALL_SNAP: f64 = 0.75;

/// Cell size of the roof surface used for splitting (m).
const SPLIT_RES: f64 = 0.5;

/// Largest relative area change accepted from the right-angled outline.
const MAX_AREA_CHANGE: f64 = 0.1;

/// Largest share of a house's points left outside an accepted outline.
const MAX_OUTSIDE: f64 = 0.05;

/// Mean deviation from the traced boundary one more outline corner is worth (m).
const CORNER_COST: f64 = 0.03;
/// Smallest hole kept as a courtyard (m²).
const MIN_COURT: f64 = 25.0;
/// Building points per m² above which a hole is roof, not a courtyard.
const COURT_DENSITY: f64 = 5.0;
/// Footprints filling their rectangle less than this get no pitched `roof:shape`.
const SIMPLE_FILL: f64 = 0.85;
/// Outlines filling less of their smallest enclosing rectangle than this, with more corners than
/// [`REVIEW_CORNERS`], are flagged for review.
const REVIEW_FILL: f64 = 0.7;
const REVIEW_CORNERS: usize = 8;

/// Share of a blob's mask cells with points below which its roof is sparse ...
const SPARSE_ROOF: f64 = 0.6;
/// ... and split on this grid (m).
const SPARSE_SPLIT_RES: f64 = 1.0;
/// Percentile of a part's point heights below which points are walls, not roof.
const EAVE_PCT: usize = 5;
/// Degrees by which the bounding box's direction must differ from the walls' to be tried as well ...
const BOX_APART: f64 = 5.0;
/// ... once the walls' outline has more corners than this, or none fits.
const BOX_CORNERS: usize = 24;
/// Margin read around the area so buildings crossing its edge are traced whole.
const MARGIN: f64 = 40.0;
/// Farthest (m) a tile's read grows past it to take in a building complex cut off at the margin.
const MAX_REACH: f64 = 300.0;

fn main() -> Result<()> {
    let args = Args::parse();
    let [x0, y0, x1, y1] = args.bbox[..] else { bail!("--bbox needs 4 values") };
    let area = [x0, y0, x1, y1];
    let q = [x0 - MARGIN, y0 - MARGIN, x1 + MARGIN, y1 + MARGIN];

    if let Some(path) = &args.bin {
        // Far enough out for a tile at the edge to grow its read over a building complex.
        let r = MAX_REACH + MARGIN;
        return bin(args.index.as_ref().unwrap(), path, [x0 - r, y0 - r, x1 + r, y1 + r], &args.files);
    }
    let (polys, props) = if let Some(path) = &args.store {
        if let Some(dsm) = &args.dsm {
            write_dsm(dsm, area, 0.25, &store::load(&Connection::open(path)?, area)?.pts)?;
        }
        let tiles = store::tiles(&Connection::open(path)?, area)?;
        eprintln!("{} tiles", tiles.len());
        let done = AtomicUsize::new(0);
        let out = Mutex::new(Vec::new());
        tiles.par_iter().try_for_each(|&(tx, ty)| -> Result<()> {
            let t = [tx as f64 * store::TILE, ty as f64 * store::TILE, (tx + 1) as f64 * store::TILE, (ty + 1) as f64 * store::TILE];
            // Only the part of the tile inside the area, so adjacent runs do not overlap.
            let t = [t[0].max(x0), t[1].max(y0), t[2].min(x1), t[3].min(y1)];
            let conn = Connection::open(path)?;
            let mut q = [t[0] - MARGIN, t[1] - MARGIN, t[2] + MARGIN, t[3] + MARGIN];
            // A building complex the margin cuts off would be traced in pieces, differently by each tile.
            let (cloud, blobs) = loop {
                let mut cloud = store::load(&conn, q)?;
                cloud.sort();
                let blobs = outline::Grid::new(q, args.res, &cloud.pts, &cloud.ground).buildings(args.min_area);
                let g = grown(q, t, &blobs, args.res);
                if g == q {
                    break (cloud, blobs);
                }
                q = g;
            };
            let r = trace_area(&cloud, blobs, t, &args);
            out.lock().unwrap().push(((tx, ty), r));
            let n = done.fetch_add(1, Ordering::Relaxed) + 1;
            if n % 100 == 0 || n == tiles.len() {
                eprintln!("{n} of {} tiles traced", tiles.len());
            }
            Ok(())
        })?;
        let mut out = out.into_inner().unwrap();
        out.sort_by_key(|(t, _)| *t);
        out.into_iter().map(|(_, r)| r).fold((Vec::new(), Vec::new()), |mut a, (p, v)| {
            a.0.extend(p);
            a.1.extend(v);
            a
        })
    } else {
        let Some(index) = &args.index else { bail!("--index or --store needed") };
        let files = flight_lines(index, q, &args.files)?;
        eprintln!("{} flight lines", files.len());
        let clouds: Vec<Cloud> = files.par_iter().map(|f| store::read_file(f, q)).collect::<Result<_>>()?;
        let mut cloud = Cloud::default();
        for c in clouds {
            cloud.append(c);
        }
        cloud.sort();
        eprintln!("{} building points, {} ground points, {} vegetation points", cloud.pts.len(), cloud.ground.len(), cloud.veg.len());
        if let Some(path) = &args.dsm {
            write_dsm(path, area, 0.25, &cloud.pts)?;
        }
        let blobs = outline::Grid::new(q, args.res, &cloud.pts, &cloud.ground).buildings(args.min_area);
        trace_area(&cloud, blobs, area, &args)
    };
    eprintln!("{} buildings", polys.len());
    if let Some(path) = &args.check {
        check(path, &polys)?;
    }
    let features: Vec<serde_json::Value> = polys
        .iter()
        .zip(props)
        .map(|(p, props)| {
            let ring = |l: &geo::LineString<f64>| l.coords().map(|c| [c.x, c.y]).collect::<Vec<_>>();
            let coords: Vec<Vec<[f64; 2]>> = std::iter::once(p.exterior()).chain(p.interiors()).map(ring).collect();
            json!({ "type": "Feature", "properties": props, "geometry": { "type": "Polygon", "coordinates": coords } })
        })
        .collect();
    let fc = json!({
        "type": "FeatureCollection",
        "crs": { "type": "name", "properties": { "name": format!("urn:ogc:def:crs:EPSG::{}", args.epsg) } },
        "features": features,
    });
    serde_json::to_writer(std::io::BufWriter::new(File::create(args.out.as_ref().unwrap())?), &fc)?;
    Ok(())
}

/// Flight lines whose extent meets `q`.
fn flight_lines(index: &PathBuf, q: [f64; 4], like: &str) -> Result<Vec<String>> {
    let conn = Connection::open(index)?;
    let files = conn
        .prepare("SELECT file FROM laz_index WHERE max_x > ?1 AND min_x < ?3 AND max_y > ?2 AND min_y < ?4 AND file LIKE ?5 ORDER BY file")?
        .query_map(rusqlite::params![q[0], q[1], q[2], q[3], like], |r| r.get(0))?
        .collect::<rusqlite::Result<_>>()?;
    Ok(files)
}

/// Sorts the points inside `q` of the flight lines over it into the tile store at `path`, skipping
/// lines already there.
fn bin(index: &PathBuf, path: &PathBuf, q: [f64; 4], like: &str) -> Result<()> {
    let conn = Connection::open(path)?;
    store::create(&conn)?;
    let files: Vec<String> =
        flight_lines(index, q, like)?.into_iter().filter(|f| !store::is_binned(&conn, f).unwrap_or(false)).collect();
    eprintln!("{} flight lines to bin", files.len());
    let conn = Mutex::new(conn);
    let done = AtomicUsize::new(0);
    files.par_iter().try_for_each(|f| -> Result<()> {
        let chunks = store::bin(f, q)?;
        store::save(&mut conn.lock().unwrap(), f, &chunks)?;
        let bytes: usize = chunks.iter().map(|(_, d)| d.len()).sum();
        eprintln!("{} of {} binned: {f} ({} MB)", done.fetch_add(1, Ordering::Relaxed) + 1, files.len(), bytes >> 20);
        Ok(())
    })
}

/// The read box `q` around tile `t` grown to take in the blobs near the tile that reach its edge.
fn grown(q: [f64; 4], t: [f64; 4], blobs: &[outline::Building], res: f64) -> [f64; 4] {
    let near = [t[0] - MARGIN / 2.0, t[1] - MARGIN / 2.0, t[2] + MARGIN / 2.0, t[3] + MARGIN / 2.0];
    let mut g = q;
    for b in blobs {
        let r = b.rings[0].iter().fold([f64::MAX, f64::MAX, f64::MIN, f64::MIN], |r, p| [r[0].min(p[0]), r[1].min(p[1]), r[2].max(p[0]), r[3].max(p[1])]);
        let cut = r[0] < q[0] + 2.0 * res || r[1] < q[1] + 2.0 * res || r[2] > q[2] - 2.0 * res || r[3] > q[3] - 2.0 * res;
        if cut && r[0] < near[2] && r[2] > near[0] && r[1] < near[3] && r[3] > near[1] {
            g = [g[0].min(r[0] - MARGIN), g[1].min(r[1] - MARGIN), g[2].max(r[2] + MARGIN), g[3].max(r[3] + MARGIN)];
        }
    }
    [g[0].max(t[0] - MAX_REACH), g[1].max(t[1] - MAX_REACH), g[2].min(t[2] + MAX_REACH), g[3].min(t[3] + MAX_REACH)]
}

/// Those of `blobs` whose mask's centroid lies in `area`, traced with `cloud`'s ground and
/// vegetation: their outlines and properties.
fn trace_area(cloud: &Cloud, blobs: Vec<outline::Building>, area: [f64; 4], args: &Args) -> (Vec<geo::Polygon<f64>>, Vec<serde_json::Value>) {
    let [x0, y0, x1, y1] = area;
    let veg = VegIndex::new(&cloud.veg);
    let ground_z = GroundIndex::new(&cloud.ground);
    let debug_at: Option<Vec<f64>> =
        std::env::var("ROOF_DEBUG").ok().map(|v| v.split(',').map(|s| s.parse().unwrap()).collect());
    let (mut polys, mut props) = (Vec::new(), Vec::new());
    for (blob, b) in blobs.into_iter().enumerate() {
        let (cx, cy) = b.centroid();
        if cx < x0 || cx >= x1 || cy < y0 || cy >= y1 {
            continue;
        }
        let debug = debug_at.as_ref().is_some_and(|d| (cx - d[0]).hypot(cy - d[1]) < d.get(2).copied().unwrap_or(15.0));
        let all: Vec<[f64; 2]> = b.points.iter().map(|p| [p[0], p[1]]).collect();
        let (parts, walls) = if args.no_split {
            (vec![b.points], Vec::new())
        } else {
            let thetas = regularize::wall_orientations(&b.rings[0], args.simplify);
            // Roof faces fitted to a few points per cell are noise; a sparse roof is split coarser.
            let seen: std::collections::BTreeSet<(i64, i64)> =
                b.points.iter().map(|p| ((p[0] / args.res).floor() as i64, (p[1] / args.res).floor() as i64)).collect();
            let share = seen.len() as f64 * args.res * args.res / b.area;
            let res = if share < SPARSE_ROOF { SPARSE_SPLIT_RES } else { SPLIT_RES };
            if debug {
                eprintln!("{:.0}% of the mask's cells have points; split at {res} m", share * 100.0);
            }
            let (parts, walls) = split::parts(&b.points, &thetas, res, MIN_PART, debug);
            // Too few faces to split by is one building, not none.
            if parts.is_empty() { (vec![b.points], walls) } else { (parts, walls) }
        };
        let mut traced: Vec<_> = parts
            .into_iter()
            .filter_map(|p| {
                let t = trace(&p, &veg, &args, debug)?;
                Some((p, t))
            })
            .collect();
        // Holes filled and parts merged can each reveal more of the other, so both repeat until settled.
        loop {
            fill_roof_holes(&mut traced, &all);
            let Some((i, j)) = misplit(&traced, debug) else { break };
            let (pj, _) = traced.swap_remove(j);
            let mut pts = std::mem::take(&mut traced[i].0);
            pts.extend(pj);
            match trace(&pts, &veg, &args, debug) {
                Some(t) => traced[i] = (pts, t),
                None => {
                    traced.swap_remove(i);
                }
            }
        }
        let first = polys.len();
        let mut tops = Vec::new();
        for (part, (mut rings, rect, pb)) in traced {
            // Small holes are gaps in the points or terraces, not courtyards.
            let outer_ring = rings.remove(0);
            rings.retain(|h| outline::ring_area(h).abs() >= MIN_COURT);
            rings.insert(0, outer_ring);
            let fill = rect_fill(&rings[0]);
            let review = fill < REVIEW_FILL && rings[0].len() > REVIEW_CORNERS;
            // Houses traced apart meet on their shared walls only once their corners there are moved onto them.
            for v in rings.iter_mut().flatten() {
                for wall in &walls {
                    let t = wall.along(v);
                    if wall.offset(v).abs() < WALL_SNAP && t > wall.t0 - 2.0 && t < wall.t1 + 2.0 {
                        *v = [wall.p[0] + t * wall.dir[0], wall.p[1] + t * wall.dir[1]];
                    }
                }
            }
            let ring = |r: Vec<[f64; 2]>| {
                let mut c: Vec<geo::Coord<f64>> = r.into_iter().map(|[x, y]| geo::Coord { x, y }).collect();
                c.push(c[0]);
                geo::LineString(c)
            };
            let mut rings = rings.into_iter();
            let outer = ring(rings.next().unwrap());
            let poly = geo::Polygon::new(outer, rings.map(ring).collect());
            let roof = roof::roof(&part, debug);
            let mut tags = serde_json::Map::new();
            tags.insert("building".into(), "yes".into());
            if let Some(g) = ground_z.lowest_around(&poly) {
                tags.insert("height".into(), json!(((roof.top - g) * 10.0).round() / 10.0));
            }
            // A pitched roof over an L or notched footprint has wings the faces alone may not show.
            let shape = roof.shape.filter(|&s| s == "flat" || fill >= SIMPLE_FILL);
            if let Some(s) = shape {
                tags.insert("roof:shape".into(), s.into());
                if s != "flat" {
                    tags.insert("roof:height".into(), json!(((roof.top - roof.eave) * 10.0).round() / 10.0));
                }
                if roof.across {
                    tags.insert("roof:orientation".into(), "across".into());
                }
            }
            polys.push(poly);
            tags.extend(json!({
                "area": (pb.area * 10.0).round() / 10.0,
                "z_max": (pb.z_max * 10.0).round() / 10.0,
                "blob": blob,
                "orthogonal": rect,
                "rect_fill": (fill * 100.0).round() / 100.0,
                "review": review,
            }).as_object().unwrap().clone());
            tops.push(roof.top);
            props.push(serde_json::Value::Object(tags));
        }
        half_roofs(&polys[first..], &mut props[first..], &tops);
    }
    for p in polys.iter_mut() {
        *p = overlap::clean(p);
    }
    overlap::resolve(&mut polys);
    for p in polys.iter_mut() {
        *p = overlap::tidy(&overlap::clean(p));
    }
    (polys, props)
}

/// Outline of one house from its points: right-angled where that fits, else a free outline.
/// Returns the rings, whether they are right-angled, and the house's mask.
fn trace(points: &[[f64; 3]], veg: &VegIndex, args: &Args, debug: bool) -> Option<(Vec<Vec<[f64; 2]>>, bool, outline::Building)> {
    let mut q = [f64::MAX, f64::MAX, f64::MIN, f64::MIN];
    for p in points {
        q = [q[0].min(p[0]), q[1].min(p[1]), q[2].max(p[0]), q[3].max(p[1])];
    }
    let q = [q[0] - 2.0, q[1] - 2.0, q[2] + 2.0, q[3] + 2.0];
    let b = outline::Grid::new(q, args.res, points, &[])
        .buildings(args.min_area)
        .into_iter()
        .max_by(|a, b| a.area.total_cmp(&b.area))?;
    let flat: Vec<[f64; 2]> = b.points.iter().map(|p| [p[0], p[1]]).collect();
    // Only canopy at roof height can hide part of the roof; a few points down the walls must not
    // bring the garden's shrubs in.
    let mut z: Vec<f64> = points.iter().map(|p| p[2]).collect();
    z.sort_by(f64::total_cmp);
    let eave = z[z.len() * EAVE_PCT / 100];
    let cover: Vec<[f64; 2]> = veg.query(q).filter(|p| p[2] > eave - 0.5).map(|p| [p[0], p[1]]).collect();
    // Right-angled houses get the cleaner rotated tracing; round and oblique ones keep a free outline.
    let thetas = regularize::wall_orientations(&b.rings[0], args.simplify);
    let rect_at = |t: &[f64]| rectilinear::best(&flat, &cover, t, args.res, args.jog, args.min_area, MAX_AREA_CHANGE);
    let rect = rect_at(&thetas);
    // Ragged edges can outvote the walls with their grid steps; the bounding box still has them.
    let boxed = regularize::box_orientation(&b.rings[0]).filter(|&t| {
        let apart = |a: f64| ((a - t) / std::f64::consts::FRAC_PI_2 - ((a - t) / std::f64::consts::FRAC_PI_2).round()).abs() * 90.0;
        thetas.iter().all(|&a| apart(a) > BOX_APART)
    });
    let jagged = rect.as_ref().is_none_or(|r| r.iter().map(Vec::len).sum::<usize>() > BOX_CORNERS);
    let boxed = boxed.filter(|_| jagged).and_then(|t| rect_at(&[t]));
    if debug {
        eprintln!("rectilinear: {}, at the box's direction: {}", rect.is_some(), boxed.is_some());
    }
    // Boundary under canopy follows the tree, not the wall, so it does not judge the outlines.
    let canopy: std::collections::BTreeSet<(i64, i64)> = cover.iter().map(|p| (p[0].floor() as i64, p[1].floor() as i64)).collect();
    let near_canopy = |p: &[f64; 2]| {
        let (x, y) = (p[0].floor() as i64, p[1].floor() as i64);
        (-1..=1).any(|dx| (-1..=1).any(|dy| canopy.contains(&(x + dx, y + dy))))
    };
    let open: Vec<[f64; 2]> = b.rings[0].iter().filter(|p| !near_canopy(p)).copied().collect();
    let judged = if open.len() * 2 > b.rings[0].len() { &open } else { &b.rings[0] };
    // The general line fit keeps walls that are not quite at right angles, where forcing them
    // would leave steps; the simpler of the two outlines wins.
    let fits = |r: &Vec<[f64; 2]>| {
        (outline::ring_area(r) / b.area - 1.0).abs() < MAX_AREA_CHANGE && rectilinear::outside(&flat, r) < MAX_OUTSIDE
    };
    // A skewed wall beside a tree is dropped as the tree's, unless the outline needs it: shrubs at
    // eave height stand along real walls too.
    let open_air = |_: &[f64; 2]| false;
    let general = [&near_canopy as &dyn Fn(&[f64; 2]) -> bool, &open_air].into_iter().find_map(|hidden| {
        let (r, ortho) = linefit::ring(&b.rings[0], args.res, hidden).filter(|(r, _)| fits(r))?;
        let holes = b.rings[1..].iter().filter_map(|h| linefit::ring(h, args.res, hidden).map(|(h, _)| h));
        Some((std::iter::once(r).chain(holes).collect::<Vec<_>>(), ortho))
    });
    let corners = |r: &Vec<Vec<[f64; 2]>>| r.iter().map(Vec::len).sum::<usize>();
    if debug {
        eprintln!("canopy: {} cover points, {} of {} boundary corners clear", cover.len(), open.len(), b.rings[0].len());
    }
    let dev = |o: &Vec<Vec<[f64; 2]>>| deviation(judged, &o[0]);
    let score = |o: &Vec<Vec<[f64; 2]>>| dev(o) + CORNER_COST * corners(o) as f64;
    if debug {
        if let Some(r) = &boxed {
            eprintln!("rect at the box's direction: {} corners, deviation {:.2}", corners(r), dev(r));
        }
    }
    let rect = match (rect, boxed) {
        (Some(r), Some(x)) => Some(if score(&x) < score(&r) { x } else { r }),
        (r, x) => r.or(x),
    };
    if debug {
        if let Some(r) = &rect {
            eprintln!("rect: {} corners, deviation {:.2}", corners(r), dev(r));
        }
        if let Some((g, o)) = &general {
            eprintln!("general (ortho {o}): {} corners, deviation {:.2}", corners(g), dev(g));
        }
    }
    // The outline closest to the traced boundary wins, each corner costing a little: a skewed house
    // keeps its skew and a round apse its facets, but a wall bent by a tree is squared.
    match (rect, general) {
        (Some(r), Some((g, ortho))) => {
            // Skew costs a corner: most houses are square, and a traced edge leans more often than a wall.
            let skew = if ortho { 0.0 } else { CORNER_COST };
            return Some(if score(&g) + skew < score(&r) { (g, ortho, b) } else { (r, true, b) });
        }
        (Some(r), None) => return Some((r, true, b)),
        (None, Some((g, ortho))) => return Some((g, ortho, b)),
        (None, None) => {}
    }
    let reg =
        |r: &Vec<[f64; 2]>| regularize::ring(r, args.simplify, args.snap_angle.to_radians(), args.merge_offset, debug);
    let outer = b.rings.first().and_then(reg)?;
    let rings = std::iter::once(outer).chain(b.rings[1..].iter().filter_map(reg)).collect();
    Some((rings, false, b))
}

/// Vegetation points bucketed for lookups around one house.
struct VegIndex(std::collections::BTreeMap<(i64, i64), Vec<[f64; 3]>>);

const VEG_CELL: f64 = 10.0;

impl VegIndex {
    fn new(pts: &[[f64; 3]]) -> Self {
        let mut m: std::collections::BTreeMap<(i64, i64), Vec<[f64; 3]>> = std::collections::BTreeMap::new();
        for p in pts {
            m.entry(((p[0] / VEG_CELL).floor() as i64, (p[1] / VEG_CELL).floor() as i64)).or_default().push(*p);
        }
        VegIndex(m)
    }

    fn query(&self, q: [f64; 4]) -> impl Iterator<Item = &[f64; 3]> {
        let (i0, j0) = ((q[0] / VEG_CELL).floor() as i64, (q[1] / VEG_CELL).floor() as i64);
        let (i1, j1) = ((q[2] / VEG_CELL).floor() as i64, (q[3] / VEG_CELL).floor() as i64);
        (i0..=i1)
            .flat_map(move |i| (j0..=j1).map(move |j| (i, j)))
            .filter_map(|k| self.0.get(&k))
            .flatten()
            .filter(move |p| p[0] >= q[0] && p[0] <= q[2] && p[1] >= q[1] && p[1] <= q[3])
    }
}

fn write_dsm(path: &std::path::Path, q: [f64; 4], res: f64, pts: &[[f64; 3]]) -> Result<()> {
    use std::io::Write;
    let w = ((q[2] - q[0]) / res).ceil() as usize;
    let h = ((q[3] - q[1]) / res).ceil() as usize;
    let mut z = vec![f64::NAN; w * h];
    for p in pts {
        if p[0] < q[0] || p[1] < q[1] {
            continue;
        }
        let (i, j) = (((p[0] - q[0]) / res) as usize, ((q[3] - p[1]) / res) as usize);
        if i < w && j < h && !(z[j * w + i] >= p[2]) {
            z[j * w + i] = p[2];
        }
    }
    let mut f = std::io::BufWriter::new(File::create(path)?);
    writeln!(f, "ncols {w}\nnrows {h}\nxllcorner {}\nyllcorner {}\ncellsize {res}\nNODATA_value -9999", q[0], q[1])?;
    for row in z.chunks(w) {
        let line: Vec<String> = row.iter().map(|v| if v.is_nan() { "-9999".into() } else { format!("{v:.2}") }).collect();
        writeln!(f, "{}", line.join(" "))?;
    }
    Ok(())
}

/// Reports, for each expected case, whether its two points fall in one building or in two.
fn check(path: &std::path::Path, polys: &[geo::Polygon<f64>]) -> Result<()> {
    use geo::Contains;
    let find = |x: f64, y: f64| polys.iter().position(|p| p.contains(&geo::Point::new(x, y)));
    let (mut pass, mut total) = (0, 0);
    for line in std::fs::read_to_string(path)?.lines() {
        let (case, note) = line.split_once('#').map_or((line, ""), |(c, n)| (c, n.trim()));
        let f: Vec<&str> = case.split_whitespace().collect();
        if f.len() != 5 {
            continue;
        }
        let v: Vec<f64> = f[1..].iter().map(|s| s.parse()).collect::<Result<_, _>>()?;
        let (a, b) = (find(v[0], v[1]), find(v[2], v[3]));
        let ok = match f[0] {
            "same" => a.is_some() && a == b,
            "split" => a.is_some() && b.is_some() && a != b,
            other => bail!("unknown case kind {other}"),
        };
        total += 1;
        pass += ok as usize;
        eprintln!("{} {} {note} ({a:?}, {b:?})", if ok { "PASS" } else { "FAIL" }, f[0]);
    }
    eprintln!("{pass} of {total} cases pass");
    Ok(())
}

/// Mean distance of a traced ring's corners from an outline.
fn deviation(traced: &[[f64; 2]], outline: &[[f64; 2]]) -> f64 {
    let n = outline.len();
    let seg = |p: [f64; 2], a: [f64; 2], b: [f64; 2]| {
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let l2 = dx * dx + dy * dy;
        let t = if l2 == 0.0 { 0.0 } else { (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2).clamp(0.0, 1.0) };
        (p[0] - a[0] - t * dx).hypot(p[1] - a[1] - t * dy)
    };
    let sum: f64 = traced
        .iter()
        .map(|&p| (0..n).map(|i| seg(p, outline[i], outline[(i + 1) % n])).fold(f64::MAX, f64::min))
        .sum();
    sum / traced.len() as f64
}

/// Area of a ring over the area of its smallest enclosing rectangle.
fn rect_fill(r: &[[f64; 2]]) -> f64 {
    let n = r.len();
    let area = outline::ring_area(r).abs();
    let best = (0..n)
        .map(|i| {
            let (a, b) = (r[i], r[(i + 1) % n]);
            let ang = (b[1] - a[1]).atan2(b[0] - a[0]);
            let (s, c) = ang.sin_cos();
            let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
            for p in r {
                let (u, v) = (p[0] * c + p[1] * s, -p[0] * s + p[1] * c);
                lo = [lo[0].min(u), lo[1].min(v)];
                hi = [hi[0].max(u), hi[1].max(v)];
            }
            (hi[0] - lo[0]) * (hi[1] - lo[1])
        })
        .fold(f64::MAX, f64::min);
    if best > 0.0 { area / best } else { 1.0 }
}

/// Parts sharing this much of the smaller one's area are one roof traced twice.
const MAX_PART_OVERLAP: f64 = 0.5;
/// Two parts that together fill their bounding rectangle this well, better than either alone, are
/// one roof the split cut along a ridge or around a wing...
const NOTCH_FILL: f64 = 0.85;
/// ... when that rectangle is little bigger than the larger part's: the smaller fills its notch.
const NOTCH_GROWTH: f64 = 1.15;
/// Parts narrower than this (m) ...
const MAX_STRIP: f64 = 4.0;
/// ... rejoin a part as high whose rectangle they complete this well.
const STRIP_FILL: f64 = 0.9;
/// A part whose roof rises less than this (m) is flat enough to be an annex.
const ANNEX_RISE: f64 = 0.6;
/// Halves meeting at the ridge need less: the ridge along their edge already tells.
const RIDGE_FILL: f64 = 0.8;
/// Roofs stepping down more than this (m) are separate buildings, such as a lean-to.
const SHED_STEP: f64 = 1.0;
/// Halves of a roof whose points within `EDGE_BAND` (m) of each other lie at most `RIDGE_DROP` (m) below
/// the ridge were split along it.
const EDGE_BAND: f64 = 0.75;
const RIDGE_DROP: f64 = 0.5;

/// A pair of parts the roof split should not have separated.
fn misplit(traced: &[Part], debug: bool) -> Option<(usize, usize)> {
    use geo::{Area, BooleanOps, Distance, Euclidean};
    let poly = |r: &[[f64; 2]]| {
        let mut c: Vec<geo::Coord<f64>> = r.iter().map(|&[x, y]| geo::Coord { x, y }).collect();
        c.push(c[0]);
        geo::Polygon::new(geo::LineString(c), vec![])
    };
    let rings: Vec<&[[f64; 2]]> = traced.iter().map(|(_, t)| &t.0[0][..]).collect();
    let polys: Vec<_> = rings.iter().map(|r| poly(r)).collect();
    let areas: Vec<f64> = polys.iter().map(|p| p.unsigned_area()).collect();
    let rects: Vec<(f64, f64)> = rings.iter().map(|r| min_rect(&[r])).collect();
    let tops: Vec<f64> = traced.iter().map(|(_, t)| t.2.z_max).collect();
    let eaves: Vec<f64> = traced
        .iter()
        .map(|(p, _)| {
            let mut z: Vec<f64> = p.iter().map(|p| p[2]).collect();
            z.sort_by(f64::total_cmp);
            z[z.len() * 5 / 100]
        })
        .collect();
    // The ridge without chimneys and antennas.
    let ridges: Vec<f64> = traced
        .iter()
        .map(|(p, _)| {
            let mut z: Vec<f64> = p.iter().map(|p| p[2]).collect();
            z.sort_by(f64::total_cmp);
            z[z.len() * 95 / 100]
        })
        .collect();
    for i in 0..polys.len() {
        for j in i + 1..polys.len() {
            let (a, b) = (&polys[i], &polys[j]);
            let inter = a.intersection(b).unsigned_area();
            if inter > MAX_PART_OVERLAP * areas[i].min(areas[j]) {
                return Some((i, j));
            }
            if Euclidean.distance(a, b) > 0.5 {
                continue;
            }
            let (rect, _) = min_rect(&[rings[i], rings[j]]);
            let fill = (areas[i] + areas[j] - inter) / rect;
            let (small, big) = if areas[i] > areas[j] { (j, i) } else { (i, j) };
            let (fi, fj) = (areas[i] / rects[i].0, areas[j] / rects[j].0);
            let fb = areas[big] / rects[big].0;
            let big = rects[big].0;
            if debug {
                eprintln!("misplit {i} {j}: fills {fi:.2} {fj:.2}, together {fill:.2}, rect growth {:.2}, ridges {:.2} {:.2}, eaves {:.2} {:.2}", rect / big, ridges[i], ridges[j], eaves[i], eaves[j]);
            }
            // A flat roof well below the other is a garage or porch built on, even where it fills a notch.
            let (lo, hi) = if ridges[i] < ridges[j] { (i, j) } else { (j, i) };
            let annex = ridges[lo] - eaves[lo] < ANNEX_RISE && ridges[hi] - ridges[lo] > SHED_STEP;
            let notch = fill >= NOTCH_FILL && fill > fb + 0.05 && rect < NOTCH_GROWTH * big && !annex;
            // Interlocking within the larger part's rectangle, but leaving another part's corner out.
            let interlock = rect < 1.02 * big && fi.max(fj) < REVIEW_FILL && fill > fi.max(fj) + 0.05;
            // A narrow strip along a part, as high as it: a hip end or eave cut off, not a lean-to.
            let strip = rects[small].1 < MAX_STRIP && (tops[i] - tops[j]).abs() < SHED_STEP && fill >= STRIP_FILL;
            // Two halves of one roof split along its ridge: the points along their common edge are the highest.
            let ridge = fill >= RIDGE_FILL && (ridges[i] - ridges[j]).abs() < SHED_STEP && {
                let near = |pts: &[[f64; 3]], other: &geo::Polygon<f64>| {
                    pts.iter().filter(|p| Euclidean.distance(&geo::Point::new(p[0], p[1]), other) < EDGE_BAND).map(|p| p[2]).collect::<Vec<_>>()
                };
                let mut z = near(&traced[i].0, b);
                z.extend(near(&traced[j].0, a));
                z.sort_by(f64::total_cmp);
                if debug && !z.is_empty() {
                    eprintln!("misplit {i} {j}: {} points along the edge, median {:.2}, ridges {:.2} {:.2}", z.len(), z[z.len() / 2], ridges[i], ridges[j]);
                }
                z.len() >= 20 && z[z.len() / 2] > ridges[i].max(ridges[j]) - RIDGE_DROP
            };
            if debug && ridge {
                eprintln!("misplit {i} {j}: split along the ridge");
            }
            if notch || interlock || strip || ridge {
                return Some((i, j));
            }
        }
    }
    None
}

/// Area and short side of the smallest rectangle around `rings` aligned to one of their edges.
fn min_rect(rings: &[&[[f64; 2]]]) -> (f64, f64) {
    let pts = || rings.iter().flat_map(|r| r.iter());
    rings
        .iter()
        .flat_map(|r| (0..r.len()).map(move |i| (r[i], r[(i + 1) % r.len()])))
        .map(|(a, b)| {
            let (s, c) = (b[1] - a[1]).atan2(b[0] - a[0]).sin_cos();
            let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
            for p in pts() {
                let (u, v) = (p[0] * c + p[1] * s, -p[0] * s + p[1] * c);
                lo = [lo[0].min(u), lo[1].min(v)];
                hi = [hi[0].max(u), hi[1].max(v)];
            }
            let (w, h) = (hi[0] - lo[0], hi[1] - lo[1]);
            (w * h, w.min(h))
        })
        .fold((f64::MAX, 0.0), |a, b| if b.0 < a.0 { b } else { a })
}

/// A part traced from its points, with its outline rings, whether they are right-angled, and its mask.
type Part = (Vec<[f64; 3]>, (Vec<Vec<[f64; 2]>>, bool, outline::Building));

/// Fills holes full of roof points: another part of the same roof enclosed by this one, not a
/// courtyard. Parts inside a filled hole are dropped.
fn fill_roof_holes(traced: &mut Vec<Part>, all: &[[f64; 2]]) {
    let mut absorbed = vec![false; traced.len()];
    for i in 0..traced.len() {
        let mut k = 1;
        while k < traced[i].1.0.len() {
            let h = &traced[i].1.0[k];
            let inside = (1.0 - rectilinear::outside(all, h)) * all.len() as f64;
            if inside / outline::ring_area(h).abs() < COURT_DENSITY {
                k += 1;
                continue;
            }
            let h = traced[i].1.0.remove(k);
            for (j, t) in traced.iter().enumerate() {
                let (x, y) = t.1.2.centroid();
                if j != i && rectilinear::outside(&[[x, y]], &h) < 0.5 {
                    absorbed[j] = true;
                }
            }
        }
    }
    let mut k = 0;
    traced.retain(|_| {
        k += 1;
        !absorbed[k - 1]
    });
}

/// Mean ground height per square metre.
struct GroundIndex(std::collections::HashMap<(i64, i64), (f64, u32)>);

/// Ground this close to a house (m) is where it stands.
const GROUND_REACH: f64 = 3.0;

impl GroundIndex {
    fn new(ground: &[[f64; 3]]) -> Self {
        let mut m: std::collections::HashMap<(i64, i64), (f64, u32)> = std::collections::HashMap::new();
        for p in ground {
            let e = m.entry((p[0].floor() as i64, p[1].floor() as i64)).or_default();
            e.0 += p[2];
            e.1 += 1;
        }
        GroundIndex(m)
    }

    /// The lowest ground next to a house: what OSM measures `height` from.
    fn lowest_around(&self, poly: &geo::Polygon<f64>) -> Option<f64> {
        use geo::{BoundingRect, Distance, Euclidean};
        let r = poly.bounding_rect()?;
        let (x0, y0) = ((r.min().x - GROUND_REACH).floor() as i64, (r.min().y - GROUND_REACH).floor() as i64);
        let (x1, y1) = ((r.max().x + GROUND_REACH).ceil() as i64, (r.max().y + GROUND_REACH).ceil() as i64);
        let mut z: Vec<f64> = (x0..x1)
            .flat_map(|x| (y0..y1).map(move |y| (x, y)))
            .filter_map(|c| {
                let &(sum, n) = self.0.get(&c)?;
                let d = Euclidean.distance(&geo::Point::new(c.0 as f64 + 0.5, c.1 as f64 + 0.5), poly);
                (d > 0.0 && d <= GROUND_REACH).then_some(sum / n as f64)
            })
            .collect();
        if z.is_empty() {
            return None;
        }
        // A low percentile, not the minimum: a ditch or a stray point is not where the house stands.
        z.sort_by(f64::total_cmp);
        Some(z[z.len() / 10])
    }
}

/// A one-sloped part against another as high is half of a roof split along its ridge, not a lean-to:
/// its shape is the whole roof's, which this part alone does not show.
fn half_roofs(polys: &[geo::Polygon<f64>], props: &mut [serde_json::Value], tops: &[f64]) {
    use geo::{Distance, Euclidean};
    for i in 0..polys.len() {
        if props[i]["roof:shape"] != "skillion" {
            continue;
        }
        let half = (0..polys.len())
            .any(|j| j != i && (tops[i] - tops[j]).abs() < SHED_STEP && Euclidean.distance(&polys[i], &polys[j]) < 0.5);
        if half {
            let p = props[i].as_object_mut().unwrap();
            p.remove("roof:shape");
            p.remove("roof:height");
        }
    }
}
