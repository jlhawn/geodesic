#!/bin/bash
# The GPU model's pace, kernel profile and memory on this machine, in one
# report for comparing machines. For each N of NS it takes the page's
# default state runs/eleven<N>_day1825.bin (fetched from SITE's parts into
# OUT/start when the checkout has none), runs one model day from it with
# scripts/spinup.mjs so that every field a continuing run carries is there,
# and from that day-1826 state runs scripts/profileGpu.mjs by pass and by
# kernel (SPLIT=1) and scripts/paceGpu.mjs with the steps awaited one at a
# time and queued BATCH to a submission. For each N of FRESH, a resolution
# with no saved run, it starts scripts/spinup.mjs from the atlas, stops it
# after FRESH_STEPS steps and times FRESH_DAYS of a model day from that
# checkpoint the same way. Where nvidia-smi exists the GPU's used memory
# and utilisation are sampled once a second through each pace run and
# their maxima reported. A failed stage is reported and the rest still run.
# The report goes to OUT/report.txt, every stage's output to OUT/<stage>.txt.
#   OUT=runs/bench scripts/gpuBench.sh
# Environment: OUT (runs/bench), NS ("64 128"), FRESH ("256"), FRESH_STEPS
# (32), FRESH_DAYS (0.25), PACE_DAYS (1; 2 at N=64), BATCH (8), SITE
# (https://gcm.echorelay.net), and whatever scripts/spinup.mjs and
# scripts/figures/figureState.mjs read (OCEAN, RADIATION, MOIST, ...).
cd "$(dirname "$0")/.."
OUT=${OUT:-runs/bench} NS=${NS-"64 128"} FRESH=${FRESH-"256"} FRESH_STEPS=${FRESH_STEPS:-32} FRESH_DAYS=${FRESH_DAYS:-0.25}
BATCH=${BATCH:-8} SITE=${SITE:-https://gcm.echorelay.net}
mkdir -p "$OUT/start"
OUT=$(cd "$OUT" && pwd)
REPORT=$OUT/report.txt
: > "$REPORT"
say() { echo "$*" | tee -a "$REPORT"; }

stage() {
  local name=$1 start=$SECONDS; shift
  if "$@" > "$OUT/$name.txt" 2>&1; then return 0; fi
  say "  $name FAILED after $((SECONDS - start)) s: $(tail -n 2 "$OUT/$name.txt" | tr '\n' ' ' | cut -c1-300)"
  return 1
}

sampled() {
  local name=$1 sampler=; shift
  if command -v nvidia-smi > /dev/null; then
    nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv,noheader,nounits -l 1 > "$OUT/$name.smi" 2> /dev/null &
    sampler=$!
  fi
  stage "$name" "$@"; local status=$?
  if [ -n "$sampler" ]; then
    kill "$sampler" 2> /dev/null; wait "$sampler" 2> /dev/null
    say "  $name: GPU memory in use at most $(awk -F', ' '$1 > m { m = $1 } END { print m + 0 }' "$OUT/$name.smi") MiB, utilisation at most $(awk -F', ' '$2 > m { m = $2 } END { print m + 0 }' "$OUT/$name.smi") %"
  fi
  return $status
}

fetch() {
  local name=$1 parts
  parts=$(curl -fsS "$SITE/runs/$name.parts.json" | node -e 'let s = ""; process.stdin.on("data", (d) => (s += d)).on("end", () => console.log(JSON.parse(s).parts.map((p) => p.file).join(" ")))') || return 1
  for part in $parts; do curl -fsS "$SITE/runs/$part" || return 1; done | gunzip > "$OUT/start/$name.bin"
}

timed() {
  local label=$1 state=$2 days=$3
  stage "profile$label" env WARM=8 STEPS=32 node scripts/profileGpu.mjs "$state" && { head -n 1 "$OUT/profile$label.txt" | cut -c1-400; sed -n '2,13p' "$OUT/profile$label.txt"; } | tee -a "$REPORT"
  stage "split$label" env WARM=8 STEPS=32 SPLIT=1 ROWS=400 node scripts/profileGpu.mjs "$state" && say "  by kernel: $OUT/split$label.txt"
  sampled "pace$label" env DAYS="$days" node scripts/paceGpu.mjs "$state" && say "  $(cat "$OUT/pace$label.txt")"
  sampled "pace${label}b$BATCH" env DAYS="$days" BATCH="$BATCH" node scripts/paceGpu.mjs "$state" && say "  $(cat "$OUT/pace${label}b$BATCH.txt")"
}

say "GPU benchmark at $(git rev-parse --short HEAD 2> /dev/null || echo '?') on $(hostname), $(date -u +%Y-%m-%dT%H:%M:%SZ), node $(node --version), $(nproc 2> /dev/null || sysctl -n hw.ncpu) cores"
command -v nvidia-smi > /dev/null && say "$(nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader)"

for N in $NS; do
  say "N=$N from eleven${N}_day1825"
  name=eleven${N}_day1825
  if [ -s "runs/$name.bin" ]; then cp "runs/$name.bin" "$OUT/start/$name.bin"; else stage "fetch$N" fetch "$name" || continue; fi
  rm -f "$OUT/start/eleven${N}_day1826.bin" "$OUT/start/eleven$N.log"
  start=$SECONDS
  stage "day$N" env N="$N" TAG="eleven$N" DAYS=1 STRATOSPHERE=1 OUT="$OUT/start" node scripts/spinup.mjs || continue
  [ -s "$OUT/start/eleven${N}_day1826.bin" ] || { say "  no day-1826 state"; continue; }
  say "  one spin-up day with its load and save: $((SECONDS - start)) s; $(grep -c '^day ' "$OUT/start/eleven$N.log") day line(s), $(grep '^day ' "$OUT/start/eleven$N.log" | tail -n 1 | cut -c1-160)"
  timed "$N" "$OUT/start/eleven${N}_day1826.bin" "${PACE_DAYS:-$([ "$N" -le 64 ] && echo 2 || echo 1)}"
done

for N in $FRESH; do
  say "N=$N from the atlas, $FRESH_STEPS steps in"
  rm -f "$OUT/start/fresh${N}"_day*.bin "$OUT/start/fresh$N.log"
  start=$SECONDS
  stage "fresh$N" env N="$N" TAG="fresh$N" DAYS=1 STOP_AFTER_STEPS="$FRESH_STEPS" MINUTES=60 STRATOSPHERE=1 OUT="$OUT/start" node scripts/spinup.mjs || continue
  state=$(ls "$OUT/start/fresh${N}"_day*_step*.bin 2> /dev/null | tail -n 1)
  [ -n "$state" ] || { say "  no checkpoint"; continue; }
  say "  fresh start to the checkpoint: $((SECONDS - start)) s ($(grep -o 'after [0-9]* s of setup' "$OUT/start/fresh$N.log" | head -n 1)), state $(du -m "$state" | cut -f1) MiB"
  timed "$N" "$state" "$FRESH_DAYS"
done
say "done in $SECONDS s"
