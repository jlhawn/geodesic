#!/bin/bash
# Spins several resolutions up in step, one at a time so that each has the
# GPU to itself: every round, each resolution runs scripts/spinup.mjs (from
# its newest runs/<PREFIX><N>_dayNNNN.bin, or fresh) to the next multiple
# of EVERY simulated days, keeping its two newest snapshots, and then
# scripts/compareStates.mjs appends the round's states side by side to
# <OUT>/<PREFIX>_compare.md. Stops at <OUT>/STOP_<PREFIX>. A resolution
# that hits NaN drops out and the others go on; any other failure is
# retried three times from the last snapshot.
#   NS="64 128" EVERY=90 PREFIX=twin scripts/pairedSpinup.sh
cd "$(dirname "$0")/.."
NS=${NS:-"64 128"} EVERY=${EVERY:-90} PREFIX=${PREFIX:-twin}
export OUT=${OUT:-$PWD/runs}
LOG=$OUT/$PREFIX.log
alive=$NS

newest() { ls "$OUT" | grep -E "^$PREFIX$1_day[0-9]+\.bin$" | sort | tail -1; }
dayOf() { local file; file=$(newest "$1"); [ -n "$file" ] && echo "$file" | sed -E 's/.*_day0*([0-9]+)\.bin$/\1/' || echo 0; }
note() { echo "$(date '+%Y-%m-%d %H:%M') $*" >> "$LOG"; }

note "paired spin-up of N=$NS every $EVERY days, commit $(git rev-parse --short HEAD)"
while [ ! -f "$OUT/STOP_$PREFIX" ] && [ -n "$alive" ]; do
  lowest=
  for n in $alive; do day=$(dayOf "$n"); [ -z "$lowest" ] || [ "$day" -lt "$lowest" ] && lowest=$day; done
  target=$(( (lowest / EVERY + 1) * EVERY ))
  for n in $alive; do
    [ -f "$OUT/STOP_$PREFIX" ] && break
    failures=0
    while [ "$(dayOf "$n")" -lt "$target" ]; do
      N=$n TAG=$PREFIX$n DAYS=$target MINUTES=100000 node scripts/spinup.mjs > "$OUT/$PREFIX$n.segment.out" 2>&1
      status=$?
      [ $status -eq 0 ] && continue
      if [ $status -eq 2 ]; then note "N=$n stopped on NaN before day $target"; alive=$(echo " $alive " | sed "s/ $n / /;s/^ *//;s/ *$//"); break; fi
      failures=$((failures + 1))
      note "N=$n failed (exit $status), attempt $failures: $(tail -3 "$OUT/$PREFIX$n.segment.out" | tr '\n' ' ')"
      if [ $failures -ge 3 ]; then note "N=$n dropped after three failures"; alive=$(echo " $alive " | sed "s/ $n / /;s/^ *//;s/ *$//"); break; fi
      sleep 30
    done
  done
  states=
  for n in $NS; do file=$OUT/$PREFIX${n}_day$(printf %04d "$target").bin; [ -f "$file" ] && states="$states $file"; done
  if [ -n "$states" ]; then
    node scripts/compareStates.mjs $states --window "$EVERY" >> "$OUT/${PREFIX}_compare.md" 2>> "$LOG"
    note "day $target done:$(for f in $states; do printf ' %s' "$(basename "$f")"; done)"
  fi
done
note "stopped ($( [ -f "$OUT/STOP_$PREFIX" ] && echo "STOP_$PREFIX" || echo 'no resolution left'))"
