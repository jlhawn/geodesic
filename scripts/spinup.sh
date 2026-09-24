#!/bin/bash
# Runs scripts/spinup.mjs segments back to back (15 wall minutes each by default, each continuing from the
# newest snapshot) until runs/STOP_<TAG> exists. A NaN stops the loop; any other failure is retried up to
# three times in a row from the last snapshot. Takes the same environment as spinup.mjs, e.g.
#   N=128 TAG=spin128 OCEAN='{"closureHours":3}' scripts/spinup.sh
cd "$(dirname "$0")/.."
export N=${N:-128} MINUTES=${MINUTES:-15}
export TAG=${TAG:-spin$N}
LOG=runs/$TAG.log
segment=0 failures=0
while [ ! -f "runs/STOP_$TAG" ]; do
  segment=$((segment + 1))
  node scripts/spinup.mjs > "runs/$TAG.segment.out" 2>&1
  status=$?
  if [ $status -eq 0 ]; then failures=0; continue; fi
  if [ $status -eq 2 ]; then echo "loop stopped on NaN in segment $segment at $(date '+%H:%M')" >> "$LOG"; exit 2; fi
  failures=$((failures + 1))
  { echo "segment $segment failed (exit $status) at $(date '+%H:%M'), attempt $failures:"; tail -5 "runs/$TAG.segment.out"; } >> "$LOG"
  [ $failures -ge 3 ] && { echo "loop stopped after three failures in a row" >> "$LOG"; exit 1; }
  sleep 30
done
echo "loop stopped by runs/STOP_$TAG at $(date '+%H:%M %Y-%m-%d') after $segment segments" >> "$LOG"
