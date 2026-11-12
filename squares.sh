#!/bin/bash
# Traces the squares of a grid over an area, each from the flight lines of every LOT around it:
#   squares.sh xmin ymin xmax ymax [size]
# Squares are aligned to multiples of size (m, default 20000) and done column by column from the
# east. Each gives work/sq/<x>_<y>.geojson and, for review in JOSM, work/sq/<x>_<y>.osm.
# Squares without flight lines and finished squares are skipped, so a stopped run resumes.
# INDEX names the flight-line index (default ~/fm/dm-maker/index.sqlite), STORE the directory for the
# squares' stores while they are worked on, up to two at a time (default work/sq), and KEEP, if set,
# the directory finished stores are moved to and taken from when a square is redone; without it
# they are deleted. BIN_ONLY=1 only bins the squares not in KEEP yet, into KEEP.
set -euo pipefail
cd "$(dirname "$0")"
read -r x0 y0 x1 y1 <<<"$1 $2 $3 $4"
size=${5:-20000}
index=${INDEX:-$HOME/fm/dm-maker/index.sqlite}
stores=${STORE:-work/sq}
keep=${KEEP:-}
bin_only=${BIN_ONLY:-}
[ -z "$bin_only" ] || [ -n "$keep" ] || { echo "BIN_ONLY needs KEEP" >&2; exit 1; }
mkdir -p work/sq "$stores" ${keep:+"$keep"}
cargo build --release --quiet
# run threads memory command…
run() { systemd-run --user --wait --quiet --pipe -p MemoryMax="$2" -p MemorySwapMax=0 -E RAYON_NUM_THREADS="$1" \
  --working-directory="$PWD" "${@:3}"; }
fl() { awk -v v="$1" -v s="$size" 'BEGIN { f = int(v / s); if (f * s > v) f--; print f * s }'; }

# Moves a finished store to KEEP, or deletes it.
put_away() {
  local store=$1 kept=$keep/${1##*/}
  if [ -z "$keep" ]; then
    rm -f "$store"*
  elif [ "$store" != "$kept" ]; then
    sqlite3 "$store" "pragma wal_checkpoint(truncate)" >/dev/null
    { cp "$store" "$kept.tmp" && mv "$kept.tmp" "$kept" && rm -f "$store"*; } || echo "keeping $store failed"
  fi
}

# Traces and exports a binned square, while the next one is binned.
finish() {
  local name=$1 bbox=$2 binned=$3 store=$4
  local start=$(date +%s)
  if [ ! -e "$name.geojson" ]; then
    run 10 10G ./target/release/roof-trace --store "$store" --bbox="$bbox" --out "$name.tmp.geojson" 2>"$name.trace.log"
    mv "$name.tmp.geojson" "$name.geojson"
  fi
  run 1 4G ./target/release/osm-export --out "$name.tmp.osm" "$name.geojson" 2>"$name.export.log"
  mv "$name.tmp.osm" "$name.osm"
  echo "$bbox: $binned trace $(($(date +%s) - start)) s, $(tail -1 "$name.trace.log"); export $(tr '\n' ' ' <"$name.export.log")"
  put_away "$store"
  ./publish.sh || echo "publishing failed"
}

prev=
for ((x = $(fl "$x1"); x >= $(fl "$x0"); x -= size)); do
  [ "$x" -lt "$x1" ] || continue
  for ((y = $(fl "$y0"); y < y1; y += size)); do
    name=work/sq/${x}_${y}
    if [ -n "$bin_only" ]; then
      [ -e "$keep/${x}_${y}.sqlite" ] && continue
    else
      [ -e "$name.osm" ] && continue
    fi
    bbox=$x,$y,$((x + size)),$((y + size))
    lines=$(sqlite3 "$index" "select count(*) from laz_index
      where max_x > $x and min_x < $((x + size)) and max_y > $y and min_y < $((y + size))")
    [ "$lines" = 0 ] && continue
    # Two stores at once must fit next to everything else on the disk.
    if [ -n "$prev" ] && [ "$(df --output=avail -B1G "$stores" | tail -1)" -lt 25 ]; then wait "$prev"; prev=; fi
    store=$stores/${x}_${y}.sqlite
    [ -n "$keep" ] && [ -e "$keep/${x}_${y}.sqlite" ] && store=$keep/${x}_${y}.sqlite
    binned="$lines flight lines,"
    if [ -n "$bin_only" ] || [ ! -e "$name.geojson" ]; then
      start=$(date +%s)
      # Binning alone can have the whole machine.
      if [ -n "$bin_only" ]; then threads=20 mem=24G; else threads=12 mem=16G; fi
      run $threads $mem ./target/release/roof-trace --index "$index" --bin "$store" --bbox="$bbox" 2>"$name.bin.log"
      binned="$binned $(du -m "$store" | cut -f1) MB, bin $(($(date +%s) - start)) s,"
    fi
    [ -n "$prev" ] && wait "$prev"
    if [ -n "$bin_only" ]; then
      echo "$bbox: ${binned%,}"
      put_away "$store" &
    else
      finish "$name" "$bbox" "$binned" "$store" &
    fi
    prev=$!
  done
done
[ -n "$prev" ] && wait "$prev"
