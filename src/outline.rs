use std::collections::BTreeMap;

/// Raster mask of building points.
pub struct Grid {
    x0: f64,
    y0: f64,
    res: f64,
    w: usize,
    h: usize,
    occ: Vec<bool>,
    z: Vec<f32>,
    pts: Vec<[f64; 3]>,
}

pub struct Building {
    /// Outer ring first (counter-clockwise), then holes (clockwise); open, in map units.
    pub rings: Vec<Vec<[f64; 2]>>,
    pub area: f64,
    pub z_max: f64,
    /// The building points that fell into it.
    pub points: Vec<[f64; 3]>,
    sum: (f64, f64, usize),
}

impl Building {
    pub fn centroid(&self) -> (f64, f64) {
        (self.sum.0 / self.sum.2 as f64, self.sum.1 / self.sum.2 as f64)
    }
}

/// Holes in the mask up to this size (m²) are gaps between points, not courtyards.
const GAP_AREA: f64 = 4.0;
/// Where less than this share of cells has points (dark roofs return few), the 3×3 closing leaves
/// gaps between them; there the closing reaches this many point spacings instead.
const SPARSE: f64 = 0.5;
const SPARSE_REACH: f64 = 1.5;
/// Cells (m) over which point density is judged.
const DENSITY_CELL: f64 = 2.0;

impl Grid {
    pub fn new(q: [f64; 4], res: f64, pts: &[[f64; 3]], ground: &[[f64; 3]]) -> Self {
        let w = ((q[2] - q[0]) / res).ceil() as usize;
        let h = ((q[3] - q[1]) / res).ceil() as usize;
        let mut g = Grid { x0: q[0], y0: q[1], res, w, h, occ: vec![false; w * h],
            z: vec![f32::MIN; w * h],
            pts: pts.to_vec(),
        };
        for p in pts {
            let i = ((p[0] - g.x0) / res) as usize;
            let j = ((p[1] - g.y0) / res) as usize;
            if i < w && j < h {
                let k = j * w + i;
                g.occ[k] = true;
                g.z[k] = g.z[k].max(p[2] as f32);
            }
        }
        let mut bare = vec![false; w * h];
        for p in ground {
            let i = ((p[0] - g.x0) / res) as usize;
            let j = ((p[1] - g.y0) / res) as usize;
            if i < w && j < h {
                bare[j * w + i] = true;
            }
        }
        // Closing bridges the spacing between points without growing the outline. It never fills a
        // cell where the ground was seen, so a narrow passage between two houses stays open.
        let reach = g.reach();
        let mut closed: BTreeMap<usize, Vec<bool>> = BTreeMap::new();
        for &r in &reach {
            closed.entry(r).or_insert_with(|| if r == 1 { g.morph(&g.morph(&g.occ, true), false) } else { g.close_round(r) });
        }
        for k in 0..g.occ.len() {
            g.occ[k] |= closed[&reach[k]][k] && !bare[k];
        }
        g.fill_gaps();
        g
    }

    /// Per cell, the closing radius (cells) that bridges the spacing of the points around it.
    fn reach(&self) -> Vec<usize> {
        let (cw, ch) = (((self.w as f64 * self.res) / DENSITY_CELL).ceil() as usize, ((self.h as f64 * self.res) / DENSITY_CELL).ceil() as usize);
        // Occupied cells rather than points: scan lines bunch points up without covering more roof.
        let mut n = vec![0usize; cw * ch];
        for k in (0..self.occ.len()).filter(|&k| self.occ[k]) {
            let (x, y) = ((k % self.w) as f64 * self.res, (k / self.w) as f64 * self.res);
            n[((y / DENSITY_CELL) as usize).min(ch - 1) * cw + ((x / DENSITY_CELL) as usize).min(cw - 1)] += 1;
        }
        // The densest neighbouring cell, so that a roof's edge cells, half off the roof, count as
        // dense as its middle.
        let radius: Vec<usize> = (0..cw * ch)
            .map(|k| {
                let (i, j) = ((k % cw) as isize, (k / cw) as isize);
                let most = (-1..=1)
                    .flat_map(|dj| (-1..=1).map(move |di| (i + di, j + dj)))
                    .filter(|&(a, b)| a >= 0 && b >= 0 && (a as usize) < cw && (b as usize) < ch)
                    .map(|(a, b)| n[b as usize * cw + a as usize])
                    .max()
                    .unwrap_or(0);
                let share = most as f64 * self.res * self.res / (DENSITY_CELL * DENSITY_CELL);
                if share == 0.0 || share >= SPARSE {
                    1
                } else {
                    // Points about res/√share apart.
                    ((SPARSE_REACH / share.sqrt()).ceil() as usize).clamp(1, 8)
                }
            })
            .collect();
        (0..self.w * self.h)
            .map(|k| {
                let (x, y) = ((k % self.w) as f64 * self.res, (k / self.w) as f64 * self.res);
                radius[((y / DENSITY_CELL) as usize).min(ch - 1) * cw + ((x / DENSITY_CELL) as usize).min(cw - 1)]
            })
            .collect()
    }

    /// Fills cells under `cover` that share a row and a column with the building: the part of a roof
    /// hidden under a tree crown. Meant for a grid aligned with one house's walls.
    pub fn complete(&mut self, cover: &[[f64; 2]]) {
        let mut row = vec![false; self.h];
        let mut col = vec![false; self.w];
        for k in 0..self.occ.len() {
            if self.occ[k] {
                row[k / self.w] = true;
                col[k % self.w] = true;
            }
        }
        for p in cover {
            let i = ((p[0] - self.x0) / self.res) as usize;
            let j = ((p[1] - self.y0) / self.res) as usize;
            if i < self.w && j < self.h && row[j] && col[i] {
                self.occ[j * self.w + i] = true;
            }
        }
        let closed = self.morph(&self.morph(&self.occ, true), false);
        for k in 0..self.occ.len() {
            self.occ[k] |= closed[k];
        }
        self.fill_gaps();
    }

    /// 3×3 dilation (`grow`) or erosion.
    fn morph(&self, m: &[bool], grow: bool) -> Vec<bool> {
        let (w, h) = (self.w as isize, self.h as isize);
        let mut out = vec![!grow; m.len()];
        for j in 0..h {
            for i in 0..w {
                let mut any = false;
                let mut all = true;
                for dj in -1..=1 {
                    for di in -1..=1 {
                        let (a, b) = (i + di, j + dj);
                        let v = a >= 0 && b >= 0 && a < w && b < h && m[(b * w + a) as usize];
                        any |= v;
                        all &= v;
                    }
                }
                out[(j * w + i) as usize] = if grow { any } else { all };
            }
        }
        out
    }

    /// Closing by a disc of radius `r` cells; a square would leave steps along the grid.
    fn close_round(&self, r: usize) -> Vec<bool> {
        let r2 = (r * r) as f64;
        let grown: Vec<bool> = edt(&self.occ, self.w, self.h).iter().map(|&d| d <= r2).collect();
        let outside: Vec<bool> = grown.iter().map(|v| !v).collect();
        edt(&outside, self.w, self.h).iter().map(|&d| d > r2).collect()
    }

    fn fill_gaps(&mut self) {
        let empty: Vec<bool> = self.occ.iter().map(|v| !v).collect();
        let max = (GAP_AREA / (self.res * self.res)) as usize;
        for comp in self.components(&empty) {
            let border = comp.iter().any(|&k| {
                let (i, j) = (k % self.w, k / self.w);
                i == 0 || j == 0 || i == self.w - 1 || j == self.h - 1
            });
            if !border && comp.len() <= max {
                for k in comp {
                    self.occ[k] = true;
                }
            }
        }
    }

    /// 4-connected components of `m`.
    fn components(&self, m: &[bool]) -> Vec<Vec<usize>> {
        let mut seen = vec![false; m.len()];
        let mut out = Vec::new();
        for start in 0..m.len() {
            if !m[start] || seen[start] {
                continue;
            }
            seen[start] = true;
            let mut comp = vec![start];
            let mut n = 0;
            while n < comp.len() {
                let k = comp[n];
                n += 1;
                let (i, j) = (k % self.w, k / self.w);
                let mut push = |k2: usize| {
                    if m[k2] && !seen[k2] {
                        seen[k2] = true;
                        comp.push(k2);
                    }
                };
                if i > 0 {
                    push(k - 1);
                }
                if i + 1 < self.w {
                    push(k + 1);
                }
                if j > 0 {
                    push(k - self.w);
                }
                if j + 1 < self.h {
                    push(k + self.w);
                }
            }
            out.push(comp);
        }
        out
    }

    pub fn buildings(&self, min_area: f64) -> Vec<Building> {
        let cell = self.res * self.res;
        let mut label = vec![usize::MAX; self.occ.len()];
        let comps: Vec<Vec<usize>> =
            self.components(&self.occ).into_iter().filter(|c| c.len() as f64 * cell >= min_area).collect();
        for (n, c) in comps.iter().enumerate() {
            for &k in c {
                label[k] = n;
            }
        }
        let mut points = vec![Vec::new(); comps.len()];
        for p in &self.pts {
            let i = ((p[0] - self.x0) / self.res) as usize;
            let j = ((p[1] - self.y0) / self.res) as usize;
            if i < self.w && j < self.h && label[j * self.w + i] != usize::MAX {
                points[label[j * self.w + i]].push(*p);
            }
        }
        comps
            .iter()
            .enumerate()
            .map(|(n, c)| {
                let mut z_max = f32::MIN;
                let mut sum = (0.0, 0.0, 0);
                for &k in c {
                    z_max = z_max.max(self.z[k]);
                    sum.0 += self.x0 + ((k % self.w) as f64 + 0.5) * self.res;
                    sum.1 += self.y0 + ((k / self.w) as f64 + 0.5) * self.res;
                    sum.2 += 1;
                }
                let mut rings = self.trace(c, &label, n);
                rings.retain(|r| ring_area(r) > 0.0 || -ring_area(r) >= min_area);
                rings.sort_by(|a, b| ring_area(b).total_cmp(&ring_area(a)));
                Building {
                    rings,
                    area: c.len() as f64 * cell,
                    z_max: z_max as f64,
                    points: std::mem::take(&mut points[n]),
                    sum,
                }
            })
            .collect()
    }

    /// Boundary rings of component `n` along cell edges, interior on the left.
    fn trace(&self, comp: &[usize], label: &[usize], n: usize) -> Vec<Vec<[f64; 2]>> {
        let inside = |i: isize, j: isize| {
            i >= 0 && j >= 0 && (i as usize) < self.w && (j as usize) < self.h && label[j as usize * self.w + i as usize] == n
        };
        // Directed edges keyed by start corner.
        let mut out: BTreeMap<(isize, isize), Vec<(isize, isize)>> = BTreeMap::new();
        for &k in comp {
            let (i, j) = ((k % self.w) as isize, (k / self.w) as isize);
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
        let mut rings = Vec::new();
        // Starting where only one edge leaves, the walk is back at its start exactly when the ring closes.
        while let Some((&start, _)) =
            out.iter().find(|(_, v)| v.len() == 1).or_else(|| out.iter().find(|(_, v)| !v.is_empty()))
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
                // At a corner shared diagonally by two cells, turning right keeps them apart, as 4-connectivity does.
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
            let ring: Vec<[f64; 2]> = corners
                .iter()
                .map(|&(i, j)| [self.x0 + i as f64 * self.res, self.y0 + j as f64 * self.res])
                .collect();
            if ring.len() >= 4 {
                rings.push(ring);
            }
        }
        rings
    }
}

/// Squared distance (cells) from each cell to the nearest set one (Felzenszwalb and Huttenlocher).
fn edt(m: &[bool], w: usize, h: usize) -> Vec<f64> {
    const FAR: f64 = 1e12;
    let mut d: Vec<f64> = m.iter().map(|&v| if v { 0.0 } else { FAR }).collect();
    let mut line = Vec::new();
    for i in 0..w {
        line.clear();
        line.extend((0..h).map(|j| d[j * w + i]));
        for (j, v) in edt_1d(&line).into_iter().enumerate() {
            d[j * w + i] = v;
        }
    }
    for j in 0..h {
        let row = edt_1d(&d[j * w..(j + 1) * w]);
        d[j * w..(j + 1) * w].copy_from_slice(&row);
    }
    d
}

/// Lower envelope of the parabolas rooted at `f`.
fn edt_1d(f: &[f64]) -> Vec<f64> {
    let n = f.len();
    let (mut v, mut z) = (vec![0usize; n], vec![0.0f64; n + 1]);
    let mut k = 0;
    z[0] = f64::NEG_INFINITY;
    z[1] = f64::INFINITY;
    let cross = |q: usize, p: usize| ((f[q] + (q * q) as f64) - (f[p] + (p * p) as f64)) / (2.0 * (q as f64 - p as f64));
    for q in 1..n {
        let mut s = cross(q, v[k]);
        while s <= z[k] {
            k -= 1;
            s = cross(q, v[k]);
        }
        k += 1;
        v[k] = q;
        z[k] = s;
        z[k + 1] = f64::INFINITY;
    }
    k = 0;
    (0..n)
        .map(|q| {
            while z[k + 1] < q as f64 {
                k += 1;
            }
            let d = q as f64 - v[k] as f64;
            d * d + f[v[k]]
        })
        .collect()
}

/// Signed area, positive for counter-clockwise.
pub fn ring_area(r: &[[f64; 2]]) -> f64 {
    let n = r.len();
    (0..n).map(|i| r[i][0] * r[(i + 1) % n][1] - r[(i + 1) % n][0] * r[i][1]).sum::<f64>() / 2.0
}
