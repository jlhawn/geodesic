#!/bin/bash
# Keeps one Verda spot instance alive for scripts/asyncSpinup.sh. Every POLL
# seconds it lists the account's instances (verda --agent vm list -o json)
# and looks at the one called NAME:
#   running                    left alone; its OS volume is remembered in
#                              STATE_FILE
#   new, ordered, provisioning, validating or deleting
#                              waited for
#   offline                    started (verda vm start)
#   absent, discontinued, failed or anything else
#                              created again as a spot instance whose OS
#                              volume is kept detached when the spot
#                              instance is discontinued: on the remembered
#                              OS volume (or OS_VOLUME) once it is detached,
#                              on the image OS while none is known yet.
# The OS volume holds the repository, node_modules and the run's files, so
# the startup script (STARTUP_SCRIPT) only has to cd there and start
# asyncSpinup.sh, which continues from the newest files. Every action and
# every change of the instance's status is logged with a timestamp to LOG
# and stdout. It stops at STOP_FILE.
#
#   NAME=gcm64 INSTANCE_TYPE=1A100.22V OS=ubuntu-24.04-cuda-12.8-open-docker \
#     SSH_KEY=<id> STARTUP_SCRIPT=<id> scripts/verdaRelaunch.sh
#
# Environment: NAME (hostname), INSTANCE_TYPE, LOCATION (FIN-01), OS (the
# image of the first instance), OS_VOLUME (a detached OS volume to start
# from, in LOCATION; otherwise the one remembered), OS_VOLUME_SIZE (100 GiB),
# SSH_KEY and STARTUP_SCRIPT (IDs from verda ssh-key list and verda
# startup-script list), POLL (120 s), STATE_FILE ($HOME/.verda-relaunch-NAME),
# LOG ($HOME/verda-relaunch-NAME.log), STOP_FILE ($HOME/STOP_verda-relaunch-NAME),
# VERDA (the CLI, verda), ROUNDS (for tests: stop after that many polls).
cd "$(dirname "$0")/.."
: "${NAME:?NAME must name the instance}" "${INSTANCE_TYPE:?INSTANCE_TYPE must be set}" "${SSH_KEY:?SSH_KEY must be set}" "${STARTUP_SCRIPT:?STARTUP_SCRIPT must be set}"
LOCATION=${LOCATION:-FIN-01} OS_VOLUME_SIZE=${OS_VOLUME_SIZE:-100} POLL=${POLL:-120} VERDA=${VERDA:-verda}
STATE_FILE=${STATE_FILE:-$HOME/.verda-relaunch-$NAME} LOG=${LOG:-$HOME/verda-relaunch-$NAME.log} STOP_FILE=${STOP_FILE:-$HOME/STOP_verda-relaunch-$NAME}
HELPER=$PWD/scripts/verdaInstances.mjs

log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a "$LOG"; }
verda_json() { local out; out=$($VERDA --agent "$@" -o json 2>&1) || { log "verda $* failed: $(echo "$out" | tr '\n' ' ' | cut -c1-300)"; return 1; }; echo "$out"; }
remembered() { [ -n "$OS_VOLUME" ] && echo "$OS_VOLUME" || { [ -f "$STATE_FILE" ] && cat "$STATE_FILE"; }; }

create() {
  local os=$OS volume volumes status out
  volume=$(remembered)
  if [ -n "$volume" ]; then
    volumes=$(verda_json volume list) || return
    status=$(echo "$volumes" | node "$HELPER" volume "$volume")
    case "$status" in
      detached) os=$volume ;;
      "") log "OS volume $volume is not listed; not creating $NAME on another (set OS_VOLUME, or remove $STATE_FILE to start from $OS)"; return ;;
      *) log "waiting for OS volume $volume to detach (it is $status)"; return ;;
    esac
  fi
  [ -n "$os" ] || { log "no OS volume known and no OS image set; not creating $NAME"; return; }
  log "creating $NAME: $INSTANCE_TYPE spot in $LOCATION on $os"
  if out=$($VERDA --agent vm create --kind gpu --instance-type "$INSTANCE_TYPE" --location "$LOCATION" --is-spot \
      --os "$os" --os-volume-size "$OS_VOLUME_SIZE" --os-volume-on-spot-discontinue keep_detached \
      --ssh-key "$SSH_KEY" --hostname "$NAME" --startup-script "$STARTUP_SCRIPT" --wait -o json 2>&1); then
    log "created $NAME: $(echo "$out" | tr '\n' ' ' | cut -c1-300)"
  else
    log "creating $NAME failed: $(echo "$out" | tr '\n' ' ' | cut -c1-300)"
  fi
}

log "keeping $NAME alive every $POLL s (state $STATE_FILE, stop at $STOP_FILE)"
seen= round=0
while [ ! -f "$STOP_FILE" ]; do
  round=$((round + 1))
  if list=$(verda_json vm list); then
    read -r id status volume <<< "$(echo "$list" | node "$HELPER" instance "$NAME")"
    [ "$status" != "$seen" ] && log "$NAME is ${status:-absent}${id:+ ($id)}"
    seen=$status
    case "$status" in
      running)
        if [ "$volume" = - ] && described=$(verda_json vm describe "$id"); then volume=$(echo "$described" | node "$HELPER" osvolume); fi
        if [ -n "$volume" ] && [ "$volume" != - ] && [ "$(cat "$STATE_FILE" 2>/dev/null)" != "$volume" ]; then
          echo "$volume" > "$STATE_FILE"
          log "remembering OS volume $volume of $NAME"
        fi ;;
      new|ordered|provisioning|validating|deleting) ;;
      offline)
        log "starting $NAME ($id)"
        verda_json vm start "$id" > /dev/null && log "started $NAME" ;;
      *) create; seen= ;;
    esac
  fi
  [ -n "$ROUNDS" ] && [ "$round" -ge "$ROUNDS" ] && break
  sleep "$POLL"
done
if [ -f "$STOP_FILE" ]; then log "stopped at $STOP_FILE"; fi
