#!/bin/bash
# Traces the squares of a grid over an area, each from the flight lines of every LOT around it:
#   squares.sh xmin ymin xmax ymax [size]
# Squares are aligned to multiples of size (m, default 20000). Each gives work/sq/<x>_<y>.geojson
# and, for review in JOSM, work/sq/<x>_<y>/{new,overlap}.osm. Squares without flight lines and
# finished squares are skipped, so a stopped run resumes. Each square's store is deleted once traced.
# INDEX names the flight-line index (default ~/fm/dm-maker/index.sqlite), OSM the existing OSM
# buildings (default work/osm/buildings.geojsonseq, see README), STORE the directory for the
# squares' temporary stores, up to two at a time (default work/sq).
set -euo pipefail
cd "$(dirname "$0")"
read -r x0 y0 x1 y1 <<<"$1 $2 $3 $4"
size=${5:-20000}
index=${INDEX:-$HOME/fm/dm-maker/index.sqlite}
osm=${OSM:-work/osm/buildings.geojsonseq}
stores=${STORE:-work/sq}
[ -e "$osm" ] || { echo "$osm missing" >&2; exit 1; }
mkdir -p work/sq "$stores"
cargo build --release --quiet
# run threads memory command…
run() { systemd-run --user --wait --quiet --pipe -p MemoryMax="$2" -p MemorySwapMax=0 -E RAYON_NUM_THREADS="$1" \
  --working-directory="$PWD" "${@:3}"; }
fl() { awk -v v="$1" -v s="$size" 'BEGIN { f = int(v / s); if (f * s > v) f--; print f * s }'; }

# Traces and exports a binned square, while the next one is binned.
finish() {
  local name=$1 bbox=$2 binned=$3 store=$4
  local start=$(date +%s)
  if [ ! -e "$name.geojson" ]; then
    run 10 10G ./target/release/roof-trace --store "$store" --bbox="$bbox" --out "$name.tmp.geojson" 2>"$name.trace.log"
    mv "$name.tmp.geojson" "$name.geojson"
  fi
  run 1 4G ./target/release/osm-export --osm "$osm" --out "$name.tmp" "$name.geojson" 2>"$name.export.log"
  rm -rf "$name"; mv "$name.tmp" "$name"
  echo "$bbox: $binned trace $(($(date +%s) - start)) s, $(tail -1 "$name.trace.log"); export $(tr '\n' ' ' <"$name.export.log")"
  rm -f "$store"*
}

prev=
for ((x = $(fl "$x0"); x < x1; x += size)); do
  for ((y = $(fl "$y0"); y < y1; y += size)); do
    name=work/sq/${x}_${y}
    [ -e "$name/overlap.osm" ] && continue
    bbox=$x,$y,$((x + size)),$((y + size))
    lines=$(sqlite3 "$index" "select count(*) from laz_index
      where max_x > $x and min_x < $((x + size)) and max_y > $y and min_y < $((y + size))")
    [ "$lines" = 0 ] && continue
    # Two stores at once must fit next to everything else on the disk.
    if [ -n "$prev" ] && [ "$(df --output=avail -B1G "$stores" | tail -1)" -lt 25 ]; then wait "$prev"; prev=; fi
    store=$stores/${x}_${y}.sqlite
    binned="$lines flight lines,"
    if [ ! -e "$name.geojson" ]; then
      start=$(date +%s)
      run 12 16G ./target/release/roof-trace --index "$index" --bin "$store" --bbox="$bbox" 2>"$name.bin.log"
      binned="$binned $(du -m "$store" | cut -f1) MB, bin $(($(date +%s) - start)) s,"
    fi
    [ -n "$prev" ] && wait "$prev"
    finish "$name" "$bbox" "$binned" "$store" &
    prev=$!
  done
done
[ -n "$prev" ] && wait "$prev"
