#!/bin/bash
# The Verda startup script for the instances of the spin-up PREFIX (register
# it with `verda --agent startup-script add --name gcm-eleven
# --file scripts/verdaElevenStartup.sh -o json` and pass its id to
# scripts/verdaRelaunch.sh; for another run register a copy whose PREFIX
# line names it, e.g. `sed 's/^PREFIX=eleven$/PREFIX=twelve/'
# scripts/verdaElevenStartup.sh > gcm-twelve.sh`, since Verda runs the
# script's text as it was registered, without an environment). Verda runs
# it as root when it creates an instance; on one recreated from the kept OS
# volume, scripts/verdaEleven.sh resume continues the run in RUNS/<PREFIX>
# if it had been started there and has neither ended nor been stopped,
# with the settings the run wrote there. On a fresh image there is no
# checkout yet and it does nothing.
PREFIX=eleven
export HOME=${HOME:-/root} PATH=/usr/local/bin:/usr/bin:/bin:/usr/local/sbin:/usr/sbin:/sbin:$PATH
REPO=${REPO:-/root/geodesic} RUNS=${RUNS:-/root/runs}
OUT=$RUNS/$PREFIX
[ -x "$REPO/scripts/verdaEleven.sh" ] || exit 0
mkdir -p "$OUT"
echo "$(date '+%Y-%m-%d %H:%M:%S') startup script on $(hostname)" >> "$OUT/$PREFIX.boot.log"
cd "$REPO" && OUT=$OUT PREFIX=$PREFIX scripts/verdaEleven.sh resume --detach >> "$OUT/$PREFIX.boot.log" 2>&1
