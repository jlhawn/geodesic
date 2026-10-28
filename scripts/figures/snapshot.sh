#!/bin/bash
# The five report figures of one saved state, one after another: the state
# maps, the ocean panels, the equatorial section, the equatorial panels and
# the mixed-layer deck, as <outdir>/<tag>_<figure>_dayNNNN.png with the
# dump's JSON and both steps' output beside it in .json and .log. The
# equatorial panels' and the deck's dumps step the state on the GPU: with
# GPU_QUEUE each waits for it at level 1 (scripts/gpuQueue.mjs share) and
# holds it for that dump alone, and is rerun once when the queue fails it
# with status 75. The other dumps and every plot run without the GPU.
# Prints each figure's path and wall time, the deck's night-side LWP check
# and its area and box statistics lines; a figure that fails prints its
# error and the rest still run, and the exit status is then 1. A figure's
# earlier outputs are removed before it runs, and an <outdir> in the
# state's own directory is refused.
#   scripts/figures/snapshot.sh <state.bin> <outdir> ["title"]
# Environment: OCEAN (JSON options for the GPU dumps' ocean,
# '{"everySteps":8}' as the run), GPU_QUEUE (scripts/gpuQueue.mjs when it
# exists; empty runs the GPU dumps directly), GPU_JOB (the start of their
# lines in the queue's status), NODE, PYTHON.
STATE=$1; OUT=$2; TITLE=$3
if [ -z "$STATE" ] || [ -z "$OUT" ]; then echo 'usage: scripts/figures/snapshot.sh <state.bin> <outdir> ["title"]' >&2; exit 2; fi
HERE=$(cd "$(dirname "$0")" && pwd)
NODE=${NODE:-node}; PYTHON=${PYTHON:-python3}
export OCEAN=${OCEAN:-'{"everySteps":8}'}
GPU_QUEUE=${GPU_QUEUE-$([ -f "$HERE/../gpuQueue.mjs" ] && echo "$(cd "$HERE/.." && pwd)/gpuQueue.mjs")}
physical() { local d=$1 rest=; case "$d" in /*) ;; *) d=$PWD/$d;; esac
  while [ ! -d "$d" ]; do rest=/$(basename "$d")$rest; d=$(dirname "$d"); done
  echo "$(cd "$d" && pwd -P)$rest"; }
STATEDIR=$(physical "$(dirname "$STATE")")
case "$(physical "$OUT")/" in "$STATEDIR"/*) echo "snapshot.sh: <outdir> $OUT lies in the state's own directory $STATEDIR; give one elsewhere" >&2; exit 2;; esac
mkdir -p "$OUT" || exit 1
NAME=$(basename "$STATE"); NAME=${NAME%.gz}; NAME=${NAME%.bin}; NAME=${NAME%.json}
TAG=${NAME%%_day*}
DAY=$("$NODE" -e 'const b=require("fs").readFileSync(process.argv[1]);if(b.toString("latin1",0,4)!=="GCMS")process.exit(1);console.log(String(JSON.parse(b.toString("utf8",8,8+b.readUInt32LE(4))).day).padStart(4,"0"))' "$STATE" 2>/dev/null)
[ -z "$DAY" ] && DAY=$(echo "$NAME" | sed -nE 's/.*_day([0-9]+).*/\1/p')
[ -z "$DAY" ] && DAY=0000
failed=0
for FIGURE in stateMaps oceanPanels eqsection eqpanels mlmdeck; do
  BASE="$OUT/${TAG}_${FIGURE}_day$DAY"
  rm -f "$BASE.json" "$BASE.log" "$BASE.png"
  start=$(date +%s)
  if [ -n "$GPU_QUEUE" ] && { [ "$FIGURE" = eqpanels ] || [ "$FIGURE" = mlmdeck ]; }; then
    for try in 1 2; do
      GPU_JOB="${GPU_JOB:-$TAG figures}: the $FIGURE dump of day $DAY" "$NODE" "$GPU_QUEUE" share "$BASE.log" "$NODE" "$HERE/$FIGURE.mjs" "$STATE" "$BASE.json" 2>> "$BASE.log"
      rc=$?
      [ $rc -eq 75 ] && [ $try = 1 ] || break
      echo "$FIGURE: GPU command-buffer aborts during the dump (status 75); rerunning it"
    done
  else
    "$NODE" "$HERE/$FIGURE.mjs" "$STATE" "$BASE.json" > "$BASE.log" 2>&1; rc=$?
  fi
  if [ $rc -ne 0 ]; then
    echo "$FIGURE: the dump failed ($BASE.log):"; tail -n 5 "$BASE.log" | sed 's/^/  /'; failed=1; continue
  fi
  if ! "$PYTHON" "$HERE/$FIGURE.py" "$BASE.json" "$BASE.png" ${TITLE:+"$TITLE"} >> "$BASE.log" 2>&1; then
    echo "$FIGURE: the plot failed ($BASE.log):"; tail -n 5 "$BASE.log" | sed 's/^/  /'; failed=1; continue
  fi
  echo "$BASE.png ($(( $(date +%s) - start )) s)"
  [ "$FIGURE" = mlmdeck ] && grep -E '^(deck:|deck boxes:|night LWP)' "$BASE.log"
  [ "$FIGURE" = oceanPanels ] && grep -E '^ocean:' "$BASE.log"
done
exit $failed
