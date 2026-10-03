#!/bin/bash
# The four report figures of one saved state, one after another, each dump
# run under the GPU lock: the state maps, the equatorial section, the
# equatorial panels and the mixed-layer deck, as
# <outdir>/<tag>_<figure>_dayNNNN.png with the dump's JSON and both
# steps' output beside it in .json and .log. Prints each figure's path and
# wall time and the deck's area and box statistics lines; a figure that
# fails prints its error and the rest still run, and the exit status is
# then 1.
#   scripts/figures/snapshot.sh <state.bin> <outdir> ["title"]
# Environment: OCEAN (JSON options for the GPU dumps' ocean,
# '{"everySteps":8}' as the run), GPULOCK (the lock wrapper, called as
# "$GPULOCK shared <command>"; none when the file does not exist), PYTHON.
STATE=$1; OUT=$2; TITLE=$3
if [ -z "$STATE" ] || [ -z "$OUT" ]; then echo 'usage: scripts/figures/snapshot.sh <state.bin> <outdir> ["title"]' >&2; exit 2; fi
HERE=$(cd "$(dirname "$0")" && pwd)
PYTHON=${PYTHON:-python3}
export OCEAN=${OCEAN:-'{"everySteps":8}'}
GPULOCK=${GPULOCK:-/private/tmp/claude-501/-Users-jlhawn-git-repos-jlhawn-geodesic/4ab28226-4f31-4147-a831-f542e02bf0fc/scratchpad/gpulock.sh}
LOCK=(); [ -x "$GPULOCK" ] && LOCK=("$GPULOCK" shared)
mkdir -p "$OUT" || exit 1
NAME=$(basename "$STATE"); NAME=${NAME%.gz}; NAME=${NAME%.bin}; NAME=${NAME%.json}
TAG=${NAME%%_day*}
DAY=$(node -e 'const b=require("fs").readFileSync(process.argv[1]);if(b.toString("latin1",0,4)!=="GCMS")process.exit(1);console.log(String(JSON.parse(b.toString("utf8",8,8+b.readUInt32LE(4))).day).padStart(4,"0"))' "$STATE" 2>/dev/null)
[ -z "$DAY" ] && DAY=$(echo "$NAME" | sed -nE 's/.*_day([0-9]+).*/\1/p')
[ -z "$DAY" ] && DAY=0000
failed=0
for FIGURE in stateMaps eqsection eqpanels mlmdeck; do
  BASE="$OUT/${TAG}_${FIGURE}_day$DAY"
  start=$(date +%s)
  if ! "${LOCK[@]}" node "$HERE/$FIGURE.mjs" "$STATE" "$BASE.json" > "$BASE.log" 2>&1; then
    echo "$FIGURE: the dump failed ($BASE.log):"; tail -n 5 "$BASE.log" | sed 's/^/  /'; failed=1; continue
  fi
  if ! "$PYTHON" "$HERE/$FIGURE.py" "$BASE.json" "$BASE.png" ${TITLE:+"$TITLE"} >> "$BASE.log" 2>&1; then
    echo "$FIGURE: the plot failed ($BASE.log):"; tail -n 5 "$BASE.log" | sed 's/^/  /'; failed=1; continue
  fi
  echo "$BASE.png ($(( $(date +%s) - start )) s)"
  [ "$FIGURE" = mlmdeck ] && grep -E '^(deck:|deck boxes:|night LWP)' "$BASE.log"
done
exit $failed
