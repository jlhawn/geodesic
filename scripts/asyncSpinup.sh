#!/bin/bash
# Asynchronous coupling at one resolution: CYCLES cycles, each a coupled
# phase of COUPLED_YEARS model years (scripts/spinup.mjs segments ending on
# the PER_YEAR snapshot days of the model year, as in pairedSpinup.sh) that
# records the ocean's forcing over its last RECORD_YEARS years, then an
# ocean-only phase of OCEAN_YEARS years looping that record
# (scripts/oceanSpinup.mjs, from the coupled phase's last snapshot). The
# next coupled phase continues the coupled run from that snapshot with the
# ocean, sea ice and sea surface of the ocean-only phase's last year
# (OCEAN_FROM): the coupled calendar never counts the ocean-only years. A
# summary line per cycle goes to <OUT>/<PREFIX><N>_cycles.log.
#
# Everything is recovered from the files in OUT, so the script can be
# stopped and started again at any point, and one driver at a time holds
# <OUT>/<PREFIX><N>_cycles.lock (another exits 0). It stops cleanly at
# <OUT>/STOP_<PREFIX> or STOP_<PREFIX><N> before its next segment, and on
# SIGTERM or SIGINT, which it passes to the running segment, once that has
# saved. A segment that fails is retried from its last snapshot after
# RETRY_WAIT seconds, up to five failures in a row; NaN (exit 2) stops it.
#
#   N=64 scripts/asyncSpinup.sh
#   N=128 CYCLES=3 OCEAN_ONLY='{"everySteps":8}' scripts/asyncSpinup.sh
#
# Environment: N (64), PREFIX (async), OUT (runs/), CYCLES (9),
# COUPLED_YEARS (10), RECORD_YEARS (1), OCEAN_YEARS (100), PER_YEAR (4),
# KEEP (4, coupled snapshots), OCEAN_KEEP (2), SNAPSHOT_DAYS (30),
# RESTORE (30 W/m²/K), OCEAN and RADIATION (JSON, both phases), OCEAN_ONLY
# (JSON ocean options for the ocean-only phase, OCEAN by default),
# CONVERGED (off; a drift in K: after an ocean-only phase whose global SST
# moved less than that over its last 10 years, no further cycle starts),
# FINAL_YEARS (0: coupled years to run after the last ocean-only phase),
# KEEP_FORCING (0: a cycle's record is deleted once its ocean-only phase is
# done), RESTORE_CMD (a shell command run once before anything else, e.g.
# to fetch OUT from a bucket), SYNC_CMD, BATCH and LAND_FROM (passed on),
# YEAR_DAYS (365; shorter years only for tests), RETRY_WAIT (30).
cd "$(dirname "$0")/.."
N=${N:-64} PREFIX=${PREFIX:-async} CYCLES=${CYCLES:-9} COUPLED_YEARS=${COUPLED_YEARS:-10} RECORD_YEARS=${RECORD_YEARS:-1} OCEAN_YEARS=${OCEAN_YEARS:-100}
PER_YEAR=${PER_YEAR:-4} KEEP=${KEEP:-4} OCEAN_KEEP=${OCEAN_KEEP:-2} SNAPSHOT_DAYS=${SNAPSHOT_DAYS:-30} RESTORE=${RESTORE:-30}
FINAL_YEARS=${FINAL_YEARS:-0} KEEP_FORCING=${KEEP_FORCING:-0} YEAR_DAYS=${YEAR_DAYS:-365} RETRY_WAIT=${RETRY_WAIT:-30}
OCEAN=${OCEAN:-'{}'} RADIATION=${RADIATION:-'{}'} OCEAN_ONLY=${OCEAN_ONLY:-$OCEAN}
export OUT=${OUT:-$PWD/runs}
unset RECORD OCEAN_FROM
TAG=$PREFIX$N
LOG=$OUT/${TAG}_cycles.log LEDGER=$OUT/${TAG}_cycles.txt
child= stopping= failures=0

note() { echo "$(date '+%Y-%m-%d %H:%M') $*" >> "$LOG"; }
sync_file() { [ -z "$SYNC_CMD" ] || /bin/sh -c "$SYNC_CMD" sync "$1" || note "SYNC_CMD failed on $1"; }
ledger() { [ -f "$LEDGER" ] && sed -n "s/^$1=//p" "$LEDGER" | tail -1; }
record() { echo "$1=$2" >> "$LEDGER"; sync_file "$LEDGER"; }
coupled_day() { ls "$OUT" | grep -E "^${TAG}_day[0-9]+\.bin$" | sed -E 's/.*_day0*([0-9]+)\.bin$/\1/' | sort -n | tail -1; }
cycle_tag() { printf '%s_c%02d' "$TAG" "$1"; }
ocean_final() { printf '%s/%s_year%04d.bin' "$OUT" "$(cycle_tag "$1")" "$OCEAN_YEARS"; }
end_of() { echo $((start + $1 * COUPLED_YEARS * YEAR_DAYS)); }
stopped() { [ -n "$stopping" ] || [ -f "$OUT/STOP_$PREFIX" ] || [ -f "$OUT/STOP_$TAG" ]; }

trap 'stopping=1; [ -n "$child" ] && kill -TERM "$child" 2>/dev/null' TERM INT
wait_child() {
  local status
  while :; do
    wait "$child"; status=$?
    kill -0 "$child" 2>/dev/null || { child=; return $status; }
  done
}
# supervise OUTPUT ENV... runs node with ENV (NAME=value words) and the
# script last, retrying failures; returns 0 when it ran to its end or was
# stopped, 2 on NaN and 1 after five failures in a row.
supervise() {
  local output=$1 status
  shift
  while :; do
    stopped && return 0
    env "$@" > "$output" 2>&1 &
    child=$!
    [ -n "$stopping" ] && kill -TERM "$child" 2>/dev/null
    wait_child; status=$?
    if [ $status -eq 0 ]; then failures=0; return 0; fi
    [ -n "$stopping" ] && return 0
    [ $status -eq 2 ] && { note "NaN: $(tail -2 "$output" | tr '\n' ' ')"; return 2; }
    failures=$((failures + 1))
    note "failed (exit $status), attempt $failures: $(tail -3 "$output" | tr '\n' ' ')"
    [ $failures -ge 5 ] && { note "five failures in a row; giving up"; return 1; }
    sleep "$RETRY_WAIT" & child=$!
    wait_child
  done
}
# coupled_segment END [RECORD_FROM HANDOFF_DAY OCEAN_FILE] runs spinup.mjs
# from the newest coupled snapshot to the next snapshot day before END,
# recording from day RECORD_FROM on and taking the ocean of OCEAN_FILE when
# it starts at HANDOFF_DAY.
coupled_segment() {
  local end=$1 from=$2 handoff=$3 ocean=$4 day target
  day=$(coupled_day); day=${day:-$start}
  target=$(awk -v d="$day" -v q="$PER_YEAR" -v y="$YEAR_DAYS" 'BEGIN { for (k = 1; ; k++) { t = int(k * y / q + 0.5); if (t > d) { print t; exit } } }')
  [ -n "$from" ] && [ "$from" -gt "$day" ] && [ "$from" -lt "$target" ] && target=$from
  [ "$end" -lt "$target" ] && target=$end
  set -- N="$N" TAG="$TAG" DAYS="$target" MINUTES=100000 KEEP="$KEEP" OCEAN="$OCEAN" RADIATION="$RADIATION"
  [ -n "$from" ] && [ "$day" -ge "$from" ] && set -- "$@" RECORD="$OUT/$(cycle_tag "$cycle")_forcing"
  if [ -n "$ocean" ] && [ "$day" -eq "$handoff" ]; then
    [ -f "$ocean" ] || { note "cycle $cycle: $ocean is missing, so the ocean cannot be handed to day $day"; return 1; }
    set -- "$@" OCEAN_FROM="$ocean"
  fi
  supervise "$OUT/$TAG.segment.out" "$@" node scripts/spinup.mjs
}
summarize() {
  local c=$1 log line ocean column sst ice
  log=$OUT/$(cycle_tag "$c").log
  ocean=$(grep "^ocean after year $OCEAN_YEARS:" "$log" 2>/dev/null | tail -1 | sed 's/^ocean after year [0-9]*: //')
  line=$(grep "^year $OCEAN_YEARS:" "$log" 2>/dev/null | tail -1 | sed 's/^year [0-9]*: //')
  column=$(echo "$line" | tr ';' '\n' | grep 'Southern Ocean' | sed 's/^ *//')
  sst=$(echo "$line" | tr ';' '\n' | grep 'global SST' | sed 's/^ *//')
  ice=$(echo "$line" | tr ';' '\n' | grep 'ice extent' | sed 's/^ *//')
  note "cycle $c: $((c * COUPLED_YEARS)) coupled years to day $(end_of "$c"), $((c * OCEAN_YEARS)) ocean-only years; ${ocean:-no ocean line}; ${column:-no column}; ${sst:-no SST}; ${ice:-no ice extent}"
  if [ -n "$CONVERGED" ]; then
    local recent
    recent=$(echo "$sst" | sed -nE 's/.*, (-?[0-9.]+) K over the last [0-9]+ years.*/\1/p')
    if [ -n "$recent" ] && awk -v d="$recent" -v t="$CONVERGED" 'BEGIN { exit !((d < 0 ? -d : d) < t) }'; then
      note "cycle $c: the SST moved $recent K over the ocean-only phase's last years, under CONVERGED=$CONVERGED; no further cycle"
      record converged "$c"
    fi
  fi
}
finish_cycle() {
  local c=$1 tag
  tag=$(cycle_tag "$c")
  summarize "$c"
  [ "$KEEP_FORCING" = 1 ] || rm -rf "$OUT/${tag}_forcing"
  for file in $(ls "$OUT" | grep -E "^${tag}_year[0-9]+(_day[0-9]+(_step[0-9]+)?)?\.bin$"); do
    [ "$OUT/$file" = "$(ocean_final "$c")" ] || rm -f "$OUT/$file"
  done
  record done "$c"
  sync_file "$LOG"
}

mkdir -p "$OUT"
LOCK=$OUT/${TAG}_cycles.lock
if ! mkdir "$LOCK" 2>/dev/null; then
  holder=$(cat "$LOCK/pid" 2>/dev/null)
  if [ -n "$holder" ] && ps -p "$holder" -o args= 2>/dev/null | grep -q asyncSpinup; then note "another driver (pid $holder) holds $LOCK; not starting"; exit 0; fi
  note "taking over $LOCK from pid ${holder:-unknown}, which is gone"
fi
trap 'rm -rf "$LOCK"' EXIT
if [ -n "$RESTORE_CMD" ]; then
  /bin/sh -c "$RESTORE_CMD" || { echo "RESTORE_CMD failed; not starting" >&2; exit 1; }
fi
echo $$ > "$LOCK/pid"
start=$(ledger start)
if [ -z "$start" ]; then
  if ls "$OUT" | grep -qE "^${TAG}_c[0-9]+_"; then note "cycle files of $TAG in $OUT but no $LEDGER; not starting"; exit 1; fi
  start=$(coupled_day); start=${start:-0}
  record start "$start"
fi
last=$CYCLES converged=$(ledger converged)
[ -n "$converged" ] && [ "$converged" -lt "$last" ] && last=$converged
note "asynchronous spin-up of $TAG from day $start: $last cycles of $COUPLED_YEARS coupled years (the last $RECORD_YEARS recorded) and $OCEAN_YEARS ocean-only years, commit $(git rev-parse --short HEAD 2>/dev/null)"

status=0 cycle=1
while [ $cycle -le "$last" ] && [ $status -eq 0 ] && ! stopped; do
  finished=$(ledger done)
  if [ -n "$finished" ] && [ "$finished" -ge "$cycle" ]; then cycle=$((cycle + 1)); continue; fi
  end=$(end_of "$cycle") day=$(coupled_day); day=${day:-$start}
  if [ "$day" -lt "$end" ]; then
    previous=
    [ "$cycle" -gt 1 ] && previous=$(ocean_final $((cycle - 1)))
    coupled_segment "$end" $((end - RECORD_YEARS * YEAR_DAYS)) "$(end_of $((cycle - 1)))" "$previous"; status=$?
    continue
  fi
  tag=$(cycle_tag "$cycle")
  if [ ! -f "$(ocean_final "$cycle")" ]; then
    if [ "$day" -gt "$end" ]; then note "cycle $cycle: the coupled run is at day $day, past $end, without $(ocean_final "$cycle")"; status=1; break; fi
    forcing=$OUT/${tag}_forcing
    count=$(ls "$forcing" 2>/dev/null | grep -cE '^forcing-[0-9]+\.bin$')
    if [ "$count" -lt $((RECORD_YEARS * YEAR_DAYS)) ]; then note "cycle $cycle: $forcing holds $count of $((RECORD_YEARS * YEAR_DAYS)) days"; status=1; break; fi
    supervise "$OUT/$tag.out" N="$N" TAG="$tag" STATE="$OUT/$(printf '%s_day%04d.bin' "$TAG" "$end")" FORCING="$forcing" YEARS="$OCEAN_YEARS" DAYS_PER_YEAR="$YEAR_DAYS" LOOP_DAYS=$((RECORD_YEARS * YEAR_DAYS)) KEEP="$OCEAN_KEEP" SNAPSHOT_DAYS="$SNAPSHOT_DAYS" RESTORE="$RESTORE" OCEAN="$OCEAN_ONLY" RADIATION="$RADIATION" node scripts/oceanSpinup.mjs; status=$?
    continue
  fi
  finish_cycle "$cycle"
  converged=$(ledger converged)
  [ -n "$converged" ] && [ "$converged" -lt "$last" ] && last=$converged
  cycle=$((cycle + 1))
done

if [ $status -eq 0 ] && ! stopped && [ "$FINAL_YEARS" -gt 0 ] && [ $cycle -gt "$last" ]; then
  handoff=$(end_of "$last") final=$(( $(end_of "$last") + FINAL_YEARS * YEAR_DAYS ))
  while [ $status -eq 0 ] && ! stopped; do
    day=$(coupled_day); [ "$day" -ge "$final" ] && break
    coupled_segment "$final" "" "$handoff" "$(ocean_final "$last")"; status=$?
  done
fi

if [ -n "$stopping" ]; then note "stopped by a signal at coupled day $(coupled_day)"; sync_file "$LOG"; exit 0; fi
if stopped; then note "stopped at STOP_$PREFIX or STOP_$TAG at coupled day $(coupled_day)"; sync_file "$LOG"; exit 0; fi
if [ $status -eq 0 ]; then note "done: $last cycles, coupled day $(coupled_day)"; fi
sync_file "$LOG"
exit $status
