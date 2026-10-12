#!/bin/bash
# The asynchronous schedule of scripts/asyncSpinup.sh for several
# resolutions, one after the other so that each has the GPU to itself:
# for every N in NS in turn, asyncSpinup.sh runs that resolution's whole
# schedule (its cycles and FINAL_YEARS) from the files in OUT, and the
# next N starts once it has ended. It takes the place of
# scripts/pairedSpinup.sh for scripts/verdaEleven.sh (DRIVER=async) and
# logs to the same <OUT>/<PREFIX>.log: "N=<n> ended: <why>" when a
# resolution's schedule ends (done, NaN, failed, or its own
# STOP_<PREFIX><N>), and "stopped (<why>)" once every resolution has ended
# or at <OUT>/STOP_<PREFIX>. A resolution with an "ended" line since the
# log's last "stopped (" line is not run again, so a start after an
# interruption carries on with the resolution in progress and then the
# rest; once the schedule has stopped, a new start gives every resolution
# its schedule again, which for a finished one runs nothing more. When
# asyncSpinup.sh was stopped by a signal without STOP_<PREFIX> in OUT, or
# another driver holds its lock, this exits without the stopped line (143
# or 1), as it does on SIGTERM or SIGINT, which it passes on to
# asyncSpinup.sh; STOP_<PREFIX> and a SIGTERM to asyncSpinup.sh stop the
# schedule at once, mid-phase, with the stopped line. It exits 0 when
# every resolution's schedule is done or STOP_<PREFIX> stopped it, 1
# otherwise.
#   NS="64 128" PREFIX=twelve CYCLES=1 COUPLED_YEARS=3 OCEAN_YEARS=20 FINAL_YEARS=2 scripts/asyncPair.sh
# Environment: NS ("64 128"), PREFIX (async), OUT (runs/), ASYNC
# (scripts/asyncSpinup.sh); everything else reaches asyncSpinup.sh, with N
# set to each resolution in turn.
cd "$(dirname "$0")/.."
NS=${NS:-"64 128"} PREFIX=${PREFIX:-async} ASYNC=${ASYNC:-scripts/asyncSpinup.sh}
export PREFIX OUT=${OUT:-$PWD/runs}
LOG=$OUT/$PREFIX.log
STAMP='^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2} '
child= stopping=

note() { echo "$(date '+%Y-%m-%d %H:%M') $*" >> "$LOG"; }
since_stop() { local at; at=$(grep -n 'stopped (' "$LOG" 2>/dev/null | tail -1 | cut -d: -f1); tail -n +"$((${at:-0} + 1))" "$LOG" 2>/dev/null; }
ended() { since_stop | grep -E "${STAMP}N=$1 ended: " | tail -1 | sed 's/.* ended: //'; }
coupled_day() { ls "$OUT" | grep -E "^$PREFIX$1_day[0-9]+\.bin$" | sed -E 's/.*_day0*([0-9]+)\.bin$/\1/' | sort -n | tail -1; }
summary() {
  local n why day all=
  for n in $NS; do
    why=$(ended "$n") day=$(coupled_day "$n")
    [ -n "$why" ] || why=$([ -n "$day" ] && echo "at coupled day $day" || echo 'not started')
    all="$all${all:+; }N=$n $why"
  done
  echo "$all"
}
interrupted() { note "N=$1 interrupted by a signal at coupled day $(coupled_day "$1" | grep . || echo 0); a new start carries on from there"; exit 143; }

trap 'stopping=1; [ -n "$child" ] && kill -TERM "$child" 2>/dev/null' TERM INT

mkdir -p "$OUT"
note "asynchronous schedule of N=$NS one after the other through $ASYNC, commit $(git rev-parse --short HEAD 2>/dev/null)"
for n in $NS; do
  [ -n "$stopping" ] && interrupted "$n"
  [ -f "$OUT/STOP_$PREFIX" ] && break
  why=$(ended "$n")
  if [ -n "$why" ]; then note "N=$n already ended ($why); next"; continue; fi
  cycles=$OUT/$PREFIX${n}_cycles.log
  before=$(cat "$cycles" 2>/dev/null | wc -l | tr -d ' ')
  note "N=$n: $ASYNC from coupled day $(coupled_day "$n" | grep . || echo 0)"
  N=$n bash "$ASYNC" &
  child=$!
  [ -n "$stopping" ] && kill -TERM "$child" 2>/dev/null
  while :; do wait "$child"; code=$?; kill -0 "$child" 2>/dev/null || break; done
  child=
  new=$(tail -n +"$((${before:-0} + 1))" "$cycles" 2>/dev/null)
  day=$(coupled_day "$n"); day=${day:-0}
  if echo "$new" | grep -qE "${STAMP}another driver \(pid"; then
    note "N=$n not run: $(echo "$new" | grep -E "${STAMP}another driver \(pid" | tail -1 | cut -c18-)"; exit 1
  fi
  if [ $code -eq 0 ] && echo "$new" | grep -qE "${STAMP}done: [0-9]+ cycles, coupled day"; then note "N=$n ended: done at coupled day $day"; continue; fi
  [ -f "$OUT/STOP_$PREFIX" ] && break
  if [ -n "$stopping" ] || echo "$new" | grep -qE "${STAMP}stopped by a signal at coupled day"; then interrupted "$n"; fi
  if [ $code -eq 0 ] && [ -f "$OUT/STOP_$PREFIX$n" ]; then why="STOP_$PREFIX$n at coupled day $day"
  elif [ $code -eq 2 ]; then why="NaN at coupled day $day"
  else why="failed (exit $code) at coupled day $day, see $PREFIX${n}_cycles.log"
  fi
  note "N=$n ended: $why"
done
if [ -f "$OUT/STOP_$PREFIX" ]; then note "stopped (STOP_$PREFIX: $(summary))"; exit 0; fi
all=$(summary)
note "stopped (schedule ended: $all)"
echo "$all" | tr ';' '\n' | grep -qvE 'done at coupled day' && exit 1
exit 0
