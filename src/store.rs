use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{BufReader, Seek, SeekFrom},
};

use anyhow::{Result, bail, ensure};
use laz::{DecompressionSelection, LasZipDecompressor, LazVlr};
use rusqlite::{Connection, params};

const GROUND: u8 = 2;
const MEDIUM_VEGETATION: u8 = 4;
const HIGH_VEGETATION: u8 = 5;
const BUILDING: u8 = 6;

/// Side of a stored tile (m).
pub const TILE: f64 = 500.0;
/// Ground and vegetation are kept only in cells of this size (m) with building points in them or
/// next to them: where a roof's outline, height and canopy cover are judged.
const NEAR_CELL: f64 = 5.0;
/// Ground is thinned to its lowest point per cell of this size (m): fine enough to keep the gap
/// between two houses open.
const GROUND_CELL: f64 = 0.25;
/// Vegetation is thinned to its highest point per cell of this size (m).
const VEG_CELL: f64 = 0.5;

/// The points tracing needs: building, ground, and medium and high vegetation.
#[derive(Default)]
pub struct Cloud {
    pub pts: Vec<[f64; 3]>,
    pub ground: Vec<[f64; 3]>,
    pub veg: Vec<[f64; 3]>,
}

impl Cloud {
    pub fn append(&mut self, mut other: Cloud) {
        self.pts.append(&mut other.pts);
        self.ground.append(&mut other.ground);
        self.veg.append(&mut other.veg);
    }

    /// Files and tiles arrive in any order; the roof split depends on point order, so fix it.
    pub fn sort(&mut self) {
        for v in [&mut self.pts, &mut self.ground, &mut self.veg] {
            v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        }
    }
}

/// The points of the LAZ file `path` inside `q`.
pub fn read_file(path: &str, q: [f64; 4]) -> Result<Cloud> {
    let mut c = Cloud::default();
    scan(path, &[BUILDING, GROUND, MEDIUM_VEGETATION, HIGH_VEGETATION], |class, p| {
        if p[0] < q[0] || p[0] > q[2] || p[1] < q[1] || p[1] > q[3] {
            return;
        }
        match class {
            BUILDING => c.pts.push(p),
            GROUND => c.ground.push(p),
            _ => c.veg.push(p),
        }
    })?;
    Ok(c)
}

/// Calls `f` with the class and position of each point of the LAZ file `path` in one of `classes`,
/// decompressing it as a stream.
fn scan(path: &str, classes: &[u8], mut f: impl FnMut(u8, [f64; 3])) -> Result<()> {
    let mut r = BufReader::with_capacity(1 << 20, File::open(path)?);
    let h = las::raw::Header::read_from(&mut r)?;
    r.seek(SeekFrom::Start(h.header_size as u64))?;
    let mut vlr = None;
    for _ in 0..h.number_of_variable_length_records {
        let v = las::raw::Vlr::read_from(&mut r, false)?;
        if v.record_id == 22204 && v.user_id.starts_with(b"laszip encoded") {
            vlr = Some(LazVlr::from_buffer(&v.data)?);
        }
    }
    let Some(vlr) = vlr else { bail!("{path}: not a LAZ file") };
    let format = h.point_data_record_format & 0x3f;
    ensure!(format <= 10, "{path}: point format {format} unsupported");
    let legacy = format < 6;
    let n = h.large_file.map(|l| l.number_of_point_records).unwrap_or(h.number_of_point_records as u64);
    r.seek(SeekFrom::Start(h.offset_to_point_data as u64))?;

    // Only the layered formats 6+ can skip fields.
    let sel = if legacy {
        DecompressionSelection::all()
    } else {
        DecompressionSelection(DecompressionSelection::Z | DecompressionSelection::CLASSIFICATION)
    };
    let mut dec = LasZipDecompressor::selective(r, vlr, sel)?;
    let len = h.point_data_record_length as usize;
    // Formats 0-5 keep the class in byte 15's low 5 bits, 6+ in byte 16.
    let class = |r: &[u8]| if legacy { r[15] & 0x1f } else { r[16] };
    let i32_at = |r: &[u8], o: usize| i32::from_le_bytes(r[o..o + 4].try_into().unwrap());

    const BATCH: u64 = 50_000;
    let mut buf = vec![0u8; BATCH as usize * len];
    let mut done = 0;
    while done < n {
        let k = (n - done).min(BATCH) as usize;
        dec.decompress_many(&mut buf[..k * len])?;
        done += k as u64;
        for r in buf[..k * len].chunks_exact(len) {
            let c = class(r);
            if !classes.contains(&c) {
                continue;
            }
            let x = i32_at(r, 0) as f64 * h.x_scale_factor + h.x_offset;
            let y = i32_at(r, 4) as f64 * h.y_scale_factor + h.y_offset;
            f(c, [x, y, i32_at(r, 8) as f64 * h.z_scale_factor + h.z_offset]);
        }
    }
    Ok(())
}

pub fn create(conn: &Connection) -> Result<()> {
    conn.execute_batch(
        "PRAGMA journal_mode = WAL;
         PRAGMA synchronous = OFF;
         CREATE TABLE IF NOT EXISTS chunk (tx INTEGER, ty INTEGER, data BLOB);
         CREATE INDEX IF NOT EXISTS chunk_tile ON chunk (tx, ty);
         CREATE TABLE IF NOT EXISTS binned (file TEXT PRIMARY KEY);",
    )?;
    Ok(())
}

pub fn is_binned(conn: &Connection, file: &str) -> Result<bool> {
    Ok(conn.query_row("SELECT COUNT(*) FROM binned WHERE file = ?1", [file], |r| r.get::<_, u32>(0))? > 0)
}

/// Sorts the points of one LAZ file inside `clip` into tiles, keeping ground and vegetation only
/// near buildings and thinned, and returns the tiles' chunks.
pub fn bin(path: &str, clip: [f64; 4]) -> Result<Vec<((i64, i64), Vec<u8>)>> {
    let inside = |p: &[f64; 3]| p[0] >= clip[0] && p[0] <= clip[2] && p[1] >= clip[1] && p[1] <= clip[3];
    let cell = |p: &[f64; 3], s: f64| ((p[0] / s).floor() as i64, (p[1] / s).floor() as i64);
    // Two passes, so a file's ground and vegetation never have to be held whole: buildings first,
    // then the rest near them.
    let mut pts = Vec::new();
    scan(path, &[BUILDING], |_, p| if inside(&p) { pts.push(p) })?;
    let built: BTreeSet<(i64, i64)> = pts.iter().map(|p| cell(p, NEAR_CELL)).collect();
    let near = |p: &[f64; 3]| {
        let (i, j) = cell(p, NEAR_CELL);
        (-1..=1).any(|di| (-1..=1).any(|dj| built.contains(&(i + di, j + dj))))
    };
    let mut ground: BTreeMap<(i64, i64), [f64; 3]> = BTreeMap::new();
    let mut veg: BTreeMap<(i64, i64), [f64; 3]> = BTreeMap::new();
    if !built.is_empty() {
        scan(path, &[GROUND, MEDIUM_VEGETATION, HIGH_VEGETATION], |class, p| {
            if !inside(&p) || !near(&p) {
                return;
            }
            if class == GROUND {
                ground.entry(cell(&p, GROUND_CELL)).and_modify(|q| if p[2] < q[2] { *q = p }).or_insert(p);
            } else {
                veg.entry(cell(&p, VEG_CELL)).and_modify(|q| if p[2] > q[2] { *q = p }).or_insert(p);
            }
        })?;
    }
    let (ground, veg): (Vec<_>, Vec<_>) = (ground.into_values().collect(), veg.into_values().collect());

    let mut tiles: BTreeMap<(i64, i64), Vec<u8>> = BTreeMap::new();
    for (pts, class) in [(&pts, BUILDING), (&ground, GROUND), (&veg, HIGH_VEGETATION)] {
        for p in pts {
            let t = cell(p, TILE);
            let buf = tiles.entry(t).or_default();
            // Offsets from the tile's corner fit f32 to well under a millimetre.
            buf.extend_from_slice(&((p[0] - t.0 as f64 * TILE) as f32).to_le_bytes());
            buf.extend_from_slice(&((p[1] - t.1 as f64 * TILE) as f32).to_le_bytes());
            buf.extend_from_slice(&(p[2] as f32).to_le_bytes());
            buf.push(class);
        }
    }
    Ok(tiles.into_iter().collect())
}

pub fn save(conn: &mut Connection, file: &str, chunks: &[((i64, i64), Vec<u8>)]) -> Result<()> {
    let tx = conn.transaction()?;
    for ((x, y), data) in chunks {
        tx.execute("INSERT INTO chunk (tx, ty, data) VALUES (?1, ?2, ?3)", params![x, y, data])?;
    }
    tx.execute("INSERT INTO binned (file) VALUES (?1)", [file])?;
    tx.commit()?;
    Ok(())
}

/// Stored points inside `q`.
pub fn load(conn: &Connection, q: [f64; 4]) -> Result<Cloud> {
    let t = |v: f64| (v / TILE).floor() as i64;
    let mut stmt = conn.prepare_cached("SELECT tx, ty, data FROM chunk WHERE tx BETWEEN ?1 AND ?2 AND ty BETWEEN ?3 AND ?4")?;
    let mut rows = stmt.query(params![t(q[0]), t(q[2]), t(q[1]), t(q[3])])?;
    let mut c = Cloud::default();
    while let Some(r) = rows.next()? {
        let (x0, y0) = (r.get::<_, i64>(0)? as f64 * TILE, r.get::<_, i64>(1)? as f64 * TILE);
        let data = r.get_ref(2)?.as_blob()?;
        for p in data.chunks_exact(13) {
            let f = |o: usize| f32::from_le_bytes(p[o..o + 4].try_into().unwrap()) as f64;
            let pt = [x0 + f(0), y0 + f(4), f(8)];
            if pt[0] < q[0] || pt[0] > q[2] || pt[1] < q[1] || pt[1] > q[3] {
                continue;
            }
            match p[12] {
                BUILDING => c.pts.push(pt),
                GROUND => c.ground.push(pt),
                _ => c.veg.push(pt),
            }
        }
    }
    Ok(c)
}

/// Tiles holding building points, as stored chunks show.
pub fn tiles(conn: &Connection, q: [f64; 4]) -> Result<Vec<(i64, i64)>> {
    let t = |v: f64| (v / TILE).floor() as i64;
    let mut stmt = conn.prepare("SELECT DISTINCT tx, ty FROM chunk WHERE tx BETWEEN ?1 AND ?2 AND ty BETWEEN ?3 AND ?4 ORDER BY tx, ty")?;
    let rows = stmt.query_map(params![t(q[0]), t(q[2] - 1e-6), t(q[1]), t(q[3] - 1e-6)], |r| Ok((r.get(0)?, r.get(1)?)))?;
    Ok(rows.collect::<rusqlite::Result<_>>()?)
}
