//! Turns traced outlines into OSM files for review in JOSM: buildings OSM lacks go to `new.osm`,
//! those overlapping an existing OSM building to `overlap.osm`. Neighbouring outlines share nodes.

use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::PathBuf,
};

use anyhow::{Context, Result};
use clap::Parser;
use geo::{Area, BooleanOps, BoundingRect, Coord, LineString, Polygon, Validation};
use proj::Proj;
use serde_json::Value;

#[derive(Parser)]
struct Args {
    /// Existing OSM buildings as GeoJSON lines (`osmium export -f geojsonseq`), in WGS 84.
    #[arg(long)]
    osm: Option<PathBuf>,

    /// Directory for `new.osm` and `overlap.osm`.
    #[arg(long)]
    out: PathBuf,

    /// Traced outlines (`roof-trace --out`), in S-JTSK [JTSK03].
    inputs: Vec<PathBuf>,
}

/// Corners of different buildings closer than this (m) become one node, and a corner this close to
/// another building's wall becomes a node of that wall too.
const SNAP: f64 = 0.4;
/// Nodes closer than this (m) to the line through their neighbours add nothing to the outline.
const STRAIGHT: f64 = 0.05;
/// Share of either building's area that an overlap must reach for a traced building to be in OSM already.
const SAME: f64 = 0.1;
/// Keys written to OSM; the rest are the tracer's diagnostics.
const KEYS: [&str; 5] = ["building", "height", "roof:shape", "roof:height", "roof:orientation"];
/// The survey every tag and the geometry come from.
const SOURCE: &str = "LLS: ÚGKK SR";
/// Cell of the spatial hashes (m).
const CELL: f64 = 50.0;

type P = [f64; 2];

struct Building {
    rings: Vec<Vec<P>>,
    tags: Vec<(String, String)>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let mut buildings = Vec::new();
    for path in &args.inputs {
        let fc: Value = serde_json::from_reader(BufReader::new(File::open(path).with_context(|| path.display().to_string())?))?;
        for f in fc["features"].as_array().into_iter().flatten() {
            let rings: Vec<Vec<P>> = f["geometry"]["coordinates"]
                .as_array()
                .into_iter()
                .flatten()
                .map(|r| {
                    let mut r: Vec<P> = r.as_array().unwrap().iter().map(|c| [c[0].as_f64().unwrap(), c[1].as_f64().unwrap()]).collect();
                    r.pop();
                    r
                })
                .collect();
            let props = &f["properties"];
            let tags = KEYS
                .iter()
                .filter_map(|&k| {
                    let v = &props[k];
                    let v = v.as_str().map(str::to_owned).or_else(|| v.as_f64().map(|n| n.to_string()))?;
                    Some((k.to_owned(), v))
                })
                .chain(std::iter::once(("source".to_owned(), SOURCE.to_owned())))
                .collect();
            buildings.push(Building { rings, tags });
        }
    }
    eprintln!("{} traced buildings", buildings.len());

    let to_wgs = Proj::new_known_crs("EPSG:8353", "EPSG:4326", None)?;
    let existing = match &args.osm {
        Some(path) => existing(path, &buildings)?,
        None => Vec::new(),
    };
    eprintln!("{} OSM buildings around them", existing.len());
    let index = Grid::new(existing.iter().map(|p| p.bounding_rect().unwrap()));
    let (mut new, mut overlap) = (Vec::new(), Vec::new());
    for b in buildings {
        let p = polygon(&b.rings);
        let r = p.bounding_rect().unwrap();
        let area = p.unsigned_area();
        let known = index.near(r).any(|i| {
            let o = &existing[i];
            let inter = p.intersection(o).unsigned_area();
            inter > SAME * area || inter > SAME * o.unsigned_area()
        });
        if known { overlap.push(b) } else { new.push(b) }
    }
    std::fs::create_dir_all(&args.out)?;
    for (name, set) in [("new.osm", new), ("overlap.osm", overlap)] {
        let (nodes, ways) = share_nodes(set);
        write(&args.out.join(name), &nodes, &ways, &to_wgs)?;
        eprintln!("{name}: {} buildings, {} nodes", ways.len(), nodes.len());
    }
    Ok(())
}

fn polygon(rings: &[Vec<P>]) -> Polygon<f64> {
    let ls = |r: &Vec<P>| {
        let mut c: Vec<Coord<f64>> = r.iter().map(|&[x, y]| Coord { x, y }).collect();
        c.push(c[0]);
        LineString(c)
    };
    Polygon::new(ls(&rings[0]), rings[1..].iter().map(ls).collect())
}

/// The OSM buildings in `path` near any of `buildings`, in the tracer's CRS.
fn existing(path: &PathBuf, buildings: &[Building]) -> Result<Vec<Polygon<f64>>> {
    let to_jtsk = Proj::new_known_crs("EPSG:4326", "EPSG:8353", None)?;
    let traced = Grid::new(buildings.iter().map(|b| polygon(&b.rings).bounding_rect().unwrap()));
    let mut out = Vec::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let line = line.trim_start_matches('\u{1e}');
        if line.is_empty() {
            continue;
        }
        let f: Value = serde_json::from_str(line)?;
        let g = &f["geometry"];
        let polys: Vec<&Value> = match g["type"].as_str() {
            Some("Polygon") => vec![&g["coordinates"]],
            Some("MultiPolygon") => g["coordinates"].as_array().unwrap().iter().collect(),
            _ => continue,
        };
        for p in polys {
            let rings: Vec<Vec<P>> = p
                .as_array()
                .unwrap()
                .iter()
                .map(|r| {
                    let mut r: Vec<P> = r
                        .as_array()
                        .unwrap()
                        .iter()
                        .filter_map(|c| to_jtsk.convert((c[0].as_f64()?, c[1].as_f64()?)).ok().map(|(x, y)| [x, y]))
                        .collect();
                    r.pop();
                    r
                })
                .filter(|r| r.len() >= 3)
                .collect();
            if rings.is_empty() {
                continue;
            }
            let p = polygon(&rings);
            if traced.near(p.bounding_rect().unwrap()).next().is_some() {
                out.push(p);
            }
        }
    }
    Ok(out)
}

/// Rectangles bucketed by the cells they cover.
struct Grid(HashMap<(i64, i64), Vec<usize>>);

impl Grid {
    fn new(rects: impl Iterator<Item = geo::Rect<f64>>) -> Self {
        let mut m: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        for (i, r) in rects.enumerate() {
            for k in cells(r) {
                m.entry(k).or_default().push(i);
            }
        }
        Grid(m)
    }

    /// Rectangles sharing a cell with `r`, each once.
    fn near(&self, r: geo::Rect<f64>) -> impl Iterator<Item = usize> + '_ {
        let mut ids: Vec<usize> = cells(r).filter_map(|k| self.0.get(&k)).flatten().copied().collect();
        ids.sort_unstable();
        ids.dedup();
        ids.into_iter()
    }
}

fn cells(r: geo::Rect<f64>) -> impl Iterator<Item = (i64, i64)> {
    let c = |v: f64| (v / CELL).floor() as i64;
    let (x0, y0, x1, y1) = (c(r.min().x), c(r.min().y), c(r.max().x), c(r.max().y));
    (x0..=x1).flat_map(move |x| (y0..=y1).map(move |y| (x, y)))
}

/// A building as rings of node indices, with its tags.
struct Way {
    rings: Vec<Vec<usize>>,
    tags: Vec<(String, String)>,
}

/// Nodes and rings of `buildings`, with corners of neighbours merged and corners on a neighbour's
/// wall added to it, so attached buildings share their walls. A building the sharing would make
/// invalid keeps its own nodes.
fn share_nodes(buildings: Vec<Building>) -> (Vec<P>, Vec<Way>) {
    let mut nodes: Vec<P> = Vec::new();
    let mut owner: Vec<usize> = Vec::new();
    let mut at: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    let key = |p: P| ((p[0] / SNAP).floor() as i64, (p[1] / SNAP).floor() as i64);
    let mut ways: Vec<Way> = Vec::new();
    for (bi, b) in buildings.iter().enumerate() {
        let mut used: Vec<usize> = Vec::new();
        let rings = b
            .rings
            .iter()
            .map(|r| {
                r.iter()
                    .map(|&p| {
                        let (kx, ky) = key(p);
                        let found = (kx - 1..=kx + 1)
                            .flat_map(|x| (ky - 1..=ky + 1).map(move |y| (x, y)))
                            .filter_map(|k| at.get(&k))
                            .flatten()
                            .copied()
                            .filter(|&n| owner[n] != bi && !used.contains(&n) && dist(nodes[n], p) < SNAP)
                            .min_by(|&a, &b| dist(nodes[a], p).total_cmp(&dist(nodes[b], p)));
                        let n = found.unwrap_or_else(|| {
                            nodes.push(p);
                            owner.push(bi);
                            at.entry(key(p)).or_default().push(nodes.len() - 1);
                            nodes.len() - 1
                        });
                        used.push(n);
                        n
                    })
                    .collect()
            })
            .collect();
        ways.push(Way { rings, tags: b.tags.clone() });
    }

    // Corners lying on another building's wall join that wall.
    let edges_of = |w: &Way, nodes: &[P]| -> Vec<geo::Rect<f64>> {
        w.rings
            .iter()
            .flat_map(|r| (0..r.len()).map(move |i| (r[i], r[(i + 1) % r.len()])))
            .map(|(a, b)| geo::Rect::new(Coord { x: nodes[a][0], y: nodes[a][1] }, Coord { x: nodes[b][0], y: nodes[b][1] }))
            .collect()
    };
    let way_index = Grid::new(ways.iter().map(|w| {
        let r = edges_of(w, &nodes).into_iter().reduce(|a, b| {
            geo::Rect::new(
                Coord { x: a.min().x.min(b.min().x), y: a.min().y.min(b.min().y) },
                Coord { x: a.max().x.max(b.max().x), y: a.max().y.max(b.max().y) },
            )
        });
        let r = r.unwrap();
        geo::Rect::new(r.min() - Coord { x: SNAP, y: SNAP }, r.max() + Coord { x: SNAP, y: SNAP })
    }));
    let originals: Vec<Vec<Vec<usize>>> = ways.iter().map(|w| w.rings.clone()).collect();
    for wi in 0..ways.len() {
        let bbox = {
            let pts: Vec<P> = originals[wi].iter().flatten().map(|&n| nodes[n]).collect();
            let (lo, hi) = pts.iter().fold(([f64::MAX; 2], [f64::MIN; 2]), |(lo, hi), p| {
                ([lo[0].min(p[0]), lo[1].min(p[1])], [hi[0].max(p[0]), hi[1].max(p[1])])
            });
            geo::Rect::new(Coord { x: lo[0] - SNAP, y: lo[1] - SNAP }, Coord { x: hi[0] + SNAP, y: hi[1] + SNAP })
        };
        let candidates: Vec<usize> = way_index
            .near(bbox)
            .filter(|&o| o != wi)
            .flat_map(|o| originals[o].iter().flatten().copied())
            .collect();
        for ring in ways[wi].rings.iter_mut() {
            let m = ring.len();
            // Each corner joins only its nearest wall, and none beside one of the ring's own corners,
            // where it would fold the wall back.
            let mut on: Vec<Vec<(f64, usize)>> = vec![Vec::new(); m];
            let mut seen = candidates.clone();
            seen.sort_unstable();
            seen.dedup();
            for n in seen {
                if ring.contains(&n) || ring.iter().any(|&r| dist(nodes[r], nodes[n]) < SNAP) {
                    continue;
                }
                let best = (0..m)
                    .filter_map(|i| {
                        let (t, d) = along(nodes[ring[i]], nodes[ring[(i + 1) % m]], nodes[n]);
                        (t > 0.0 && t < 1.0 && d < SNAP).then_some((i, t, d))
                    })
                    .min_by(|x, y| x.2.total_cmp(&y.2));
                if let Some((i, t, _)) = best {
                    on[i].push((t, n));
                }
            }
            let mut out = Vec::with_capacity(m);
            for (i, mut on) in on.into_iter().enumerate() {
                out.push(ring[i]);
                on.sort_by(|x, y| x.0.total_cmp(&y.0));
                out.extend(on.into_iter().map(|(_, n)| n));
            }
            *ring = out;
        }
    }

    // Sharing moves corners by up to SNAP; where that breaks a building, it keeps its own corners.
    let mut own = 0;
    for (wi, w) in ways.iter_mut().enumerate() {
        for r in w.rings.iter_mut() {
            tidy(r, &nodes);
        }
        let shared = polygon(&w.rings.iter().map(|r| r.iter().map(|&n| nodes[n]).collect()).collect::<Vec<_>>());
        if shared.is_valid() {
            continue;
        }
        own += 1;
        w.rings = buildings[wi]
            .rings
            .iter()
            .map(|r| {
                r.iter()
                    .map(|&p| {
                        nodes.push(p);
                        nodes.len() - 1
                    })
                    .collect()
            })
            .collect();
    }
    if own > 0 {
        eprintln!("{own} buildings keep their own corners");
    }
    // Only nodes some ring still uses.
    let mut map = vec![usize::MAX; nodes.len()];
    let mut kept = Vec::new();
    for w in ways.iter_mut() {
        for r in w.rings.iter_mut() {
            for n in r.iter_mut() {
                if map[*n] == usize::MAX {
                    map[*n] = kept.len();
                    kept.push(nodes[*n]);
                }
                *n = map[*n];
            }
        }
    }
    (kept, ways)
}

/// Drops nodes that merging left straight between their neighbours, on a spike without width or on
/// top of the previous node.
fn tidy(r: &mut Vec<usize>, nodes: &[P]) {
    while r.len() > 3 {
        let n = r.len();
        let Some(i) = (0..n).find(|&i| {
            let (a, b, c) = (nodes[r[(i + n - 1) % n]], nodes[r[i]], nodes[r[(i + 1) % n]]);
            let (_, off) = along(a, c, b);
            dist(a, c) < STRAIGHT || off < STRAIGHT || dist(a, b) < STRAIGHT
        }) else {
            return;
        };
        r.remove(i);
    }
}

fn dist(a: P, b: P) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1])
}

/// Position of `p` along the segment `a`–`b` (0 at `a`, 1 at `b`) and its distance from the line.
fn along(a: P, b: P, p: P) -> (f64, f64) {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let l2 = dx * dx + dy * dy;
    if l2 == 0.0 {
        return (-1.0, f64::MAX);
    }
    let t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2;
    (t, ((p[0] - a[0]) * dy - (p[1] - a[1]) * dx).abs() / l2.sqrt())
}

fn write(path: &PathBuf, nodes: &[P], ways: &[Way], to_wgs: &Proj) -> Result<()> {
    let mut o = BufWriter::new(File::create(path)?);
    writeln!(o, "<?xml version='1.0' encoding='UTF-8'?>")?;
    writeln!(o, "<osm version='0.6' generator='roof-trace' upload='never'>")?;
    for (i, &[x, y]) in nodes.iter().enumerate() {
        let (lon, lat) = to_wgs.convert((x, y))?;
        writeln!(o, "  <node id='-{}' lat='{lat:.7}' lon='{lon:.7}'/>", i + 1)?;
    }
    let mut id = 0;
    let mut next = || {
        id += 1;
        id
    };
    let tag = |o: &mut BufWriter<File>, tags: &[(String, String)]| -> Result<()> {
        for (k, v) in tags {
            writeln!(o, "    <tag k='{k}' v='{}'/>", v.replace('&', "&amp;").replace('\'', "&apos;").replace('<', "&lt;"))?;
        }
        Ok(())
    };
    for w in ways {
        let ids: Vec<i64> = w
            .rings
            .iter()
            .map(|r| -> Result<i64> {
                let id = next();
                writeln!(o, "  <way id='-{id}'>")?;
                for &n in r.iter().chain(std::iter::once(&r[0])) {
                    writeln!(o, "    <nd ref='-{}'/>", n + 1)?;
                }
                if w.rings.len() == 1 {
                    tag(&mut o, &w.tags)?;
                }
                writeln!(o, "  </way>")?;
                Ok(id)
            })
            .collect::<Result<_>>()?;
        if ids.len() > 1 {
            writeln!(o, "  <relation id='-{}'>", next())?;
            for (k, id) in ids.iter().enumerate() {
                writeln!(o, "    <member type='way' ref='-{id}' role='{}'/>", if k == 0 { "outer" } else { "inner" })?;
            }
            writeln!(o, "    <tag k='type' v='multipolygon'/>")?;
            tag(&mut o, &w.tags)?;
            writeln!(o, "  </relation>")?;
        }
    }
    writeln!(o, "</osm>")?;
    Ok(())
}
