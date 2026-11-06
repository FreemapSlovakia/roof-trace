#!/bin/bash
# Publishes the finished squares of work/sq with publish/README.md and an overview of them,
# done.geojson (WGS 84), to DEST (default fm6:/fm/storage2/download/lidar-budovy/).
set -euo pipefail
cd "$(dirname "$0")"
dest=${DEST:-fm6:/fm/storage2/download/lidar-budovy/}
size=20000
cd work/sq
for d in */; do
  d=${d%/}
  [ -e "$d/overlap.osm" ] || continue
  IFS=_ read -r x y <<<"$d"
  # 5 points per side: the square's sides are curved in WGS 84.
  coords=$(for i in 0 1 2 3; do echo "$((x + i * size / 4)) $y"; done
    for i in 0 1 2 3; do echo "$((x + size)) $((y + i * size / 4))"; done
    for i in 0 1 2 3; do echo "$((x + size - i * size / 4)) $((y + size))"; done
    for i in 0 1 2 3 4; do echo "$x $((y + size - i * size / 4))"; done)
  coords=$(cs2cs -f %.7f EPSG:8353 EPSG:4326 <<<"$coords" | awk '{ printf "%s[%s,%s]", (NR > 1 ? "," : ""), $2, $1 }')
  echo "{\"type\":\"Feature\",\"properties\":{\"dir\":\"$d\"},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[$coords]]}}"
done | paste -sd, | sed 's/^/{"type":"FeatureCollection","features":[/; s/$/]}/' >done.geojson
files=(./done.geojson)
for d in */; do
  d=${d%/}
  [ -e "$d/overlap.osm" ] && files+=("./$d/new.osm" "./$d/overlap.osm" "./$d.geojson")
done
rsync -a --relative "${files[@]}" "$dest"
rsync -a ../../publish/README.md "$dest"
