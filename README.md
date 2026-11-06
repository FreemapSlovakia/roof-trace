# roof-trace

Traces building outlines from the building points (class 6) of the Slovak national lidar survey
of ÚGKK SR, for adding and improving buildings in OpenStreetMap. Outlines follow the roof edge,
not the walls, and are regularized: right-angled where the building is, with free edges where it
is not.

Attached houses are split where their roofs differ (height, slope, ridge direction, a valley
between them). Splitting errs towards too many parts: merging two ways in an editor is quicker than
finding where one should be cut.

Each building gets `building=yes`, `height` (ridge above the lowest ground around it) and, where
the roof is one of the simple shapes, `roof:shape`, `roof:height` and `roof:orientation`.

## Build

```sh
cargo build --release
```

`osm-export` needs PROJ (`libproj-dev`).

## Input

Flight-line LAZ files in S-JTSK [JTSK03] (EPSG:8353), and an index of their extents: an SQLite
table `laz_index (min_x, max_x, min_y, max_y, file)` as made by
[laz2geotiff](https://github.com/FreemapSlovakia/laz2geotiff) `index`.

## Usage

Small areas can be read straight from the flight lines:

```sh
roof-trace --index index.sqlite --bbox=-296000,-1228800,-295600,-1228400 --out buildings.geojson
```

The flight lines are long strips of several hundred MB each, without spatial indexes, so for large
areas the points are first sorted into a tile store, which is then traced tile by tile:

```sh
roof-trace --index index.sqlite --bin area.sqlite --bbox=…
roof-trace --store area.sqlite --bbox=… --out buildings.geojson
```

`squares.sh` does both for each square of a grid over an area, from the flight lines of every LOT
around it, and deletes each square's store when done; a stopped run resumes:

```sh
INDEX=index.sqlite ./squares.sh -340000 -1200000 -320000 -1180000   # → work/sq/<x>_<y>.geojson, work/sq/<x>_<y>/*.osm
```

A building belongs to the area its centroid falls in, so the outputs of neighbouring areas do not
overlap.

### Inspecting

- `--dsm dsm.asc` writes the highest building point per 0.25 m cell.
- `ROOF_DEBUG=x,y[,radius]` prints how the buildings near a point were split and traced.
- `--check tests/<area>.cases` checks pairs of points that must end up in the same building or in
  different ones (`same|split x1 y1 x2 y2`).

## OSM export

```sh
osmium tags-filter -o buildings.osm.pbf slovakia.osm.pbf w/building r/building
osmium export -f geojsonseq -x print_record_separator=false --geometry-types=polygon \
  -o buildings.geojsonseq buildings.osm.pbf
osm-export --osm buildings.geojsonseq --out osm/ work/sq/*.geojson
```

This writes `osm/new.osm`, the buildings OSM lacks, and `osm/overlap.osm`, those overlapping an
existing OSM building, for comparing by hand. Corners of neighbouring buildings closer than 0.4 m
become one node, and a corner on a neighbour's wall becomes a node of that wall, so attached
buildings share their walls. Only the OSM tags above are written, with `source=LLS: ÚGKK SR`.

The files are marked `upload="never"`: they are for review in JOSM, not for uploading unchecked.

`publish.sh` copies the finished squares of `squares.sh`, an overview `done.geojson` and
`publish/README.md` (in Slovak, for the users of the files) to download.freemap.sk.

## Licence

GPL-3.0-or-later. Outlines derived from the lidar data must credit "ÚGKK SR".
