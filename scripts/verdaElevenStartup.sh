#!/bin/bash
# The Verda startup script for the instances of the paired spin-up eleven
# (register it with `verda --agent startup-script add --name gcm-eleven
# --file scripts/verdaElevenStartup.sh -o json` and pass its id to
# scripts/verdaRelaunch.sh). Verda runs it as root when it creates an
# instance; on one recreated from the kept OS volume, scripts/verdaEleven.sh
# resume continues the run if it had been started there and has neither
# ended nor been stopped. On a fresh image there is no checkout yet and it
# does nothing.
export HOME=${HOME:-/root} PATH=/usr/local/bin:/usr/bin:/bin:/usr/local/sbin:/usr/sbin:/sbin:$PATH
REPO=/root/geodesic
[ -x "$REPO/scripts/verdaEleven.sh" ] || exit 0
mkdir -p /root/runs/eleven
echo "$(date '+%Y-%m-%d %H:%M:%S') startup script on $(hostname)" >> /root/runs/eleven/eleven.boot.log
cd "$REPO" && OUT=/root/runs/eleven scripts/verdaEleven.sh resume --detach >> /root/runs/eleven/eleven.boot.log 2>&1
