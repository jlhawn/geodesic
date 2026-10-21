#!/bin/bash
# The instance side of a spin-up on a Verda GPU instance, paired (PREFIX
# eleven: NS "64 128" on bl36 from the atlas, PER_YEAR 36, three years)
# or asynchronous (DRIVER=async, PREFIX twelve: each N in NS in turn
# through the cycles of scripts/asyncSpinup.sh), run from the checkout
# scripts/verdaPush.sh puts on the instance's OS volume. Every mode can be
# run again: after an eviction the kept OS volume still holds the
# installation, the reports and the run. Modes:
#   setup      installs what is missing (node 22 or later from NodeSource,
#              rsync, the Vulkan loader and vulkan-tools; node_modules by
#              npm ci, with build-essential and python3 only if that
#              fails), then checks that node's webgpu reaches the GPU:
#              nvidia-smi, vulkaninfo's device and the adapter's vendor
#              (EXPECT_VENDOR, nvidia on Linux), into OUT/<PREFIX>_gpu.txt;
#              and autostart
#   autostart  where systemd runs (Linux, or SYSTEMD_DIR set), installs and
#              enables gcm-<PREFIX>-resume.service in SYSTEMD_DIR
#              (/etc/systemd/system) through SYSTEMCTL (systemctl): resume
#              at every boot of the OS volume, a restart of the same
#              instance included, on which Verda's startup script may not
#              run again
#   suite      scripts/suiteReport.sh: every test file as a concurrent
#              runner, failures with their output in OUT/<PREFIX>_suite.txt
#   bench      scripts/verdaBenchmark.mjs: fresh starts of CASES timed, the
#              run's hours and cost at PRICE in OUT/<PREFIX>_benchmark.txt
#   bootstrap  setup, suite and bench
#   run        the run through PAIRED: with DRIVER=paired the paired run
#              of scripts/pairedSpinup.sh to day UNTIL; with DRIVER=async
#              scripts/asyncPair.sh, which gives each N in NS in turn the
#              schedule of scripts/asyncSpinup.sh (CYCLES cycles of
#              COUPLED_YEARS coupled years recording the last
#              RECORD_YEARS, each followed by OCEAN_YEARS ocean-only
#              years, then FINAL_YEARS coupled years), UNTIL being then
#              the coupled day the schedule ends on. It starts once the
#              disk has room for every file to the end (for async: KEEP
#              snapshots, CYCLES + OCEAN_KEEP ocean-only year files and the
#              forcing of a cycle, about 2.3 MB a day at N=64 and 9 at
#              N=128, or of every cycle with KEEP_FORCING=1); every snapshot
#              stays (KEEP) for the Mac to pull, each saved file
#              flushed to the disk (SYNC_CMD 'sync "$1"') so that a power
#              cut cannot leave a state named but empty. It marks OUT with
#              STARTED_<PREFIX>, and with ENDED_<PREFIX> only when the
#              driver logs its own stop ("stopped (…)" in <PREFIX>.log: day
#              UNTIL reached, every resolution's schedule ended,
#              STOP_<PREFIX>, or no resolution left); killed by a signal,
#              as at an eviction's shutdown, it leaves no mark
#   resume     what the startup script and the autostart unit run at boot
#              (a second one exits on the lock): if run was
#              started here and has neither ended nor been stopped, waits
#              up to GPU_WAIT seconds for the GPU and runs on from the
#              newest files (for async: the resolution in progress, then
#              the ones after it); otherwise nothing
#   status     the newest state of each N (for async with its cycle and
#              phase and the last line of its cycles log), the run's last
#              log lines, the markers and the free disk
# `--detach` puts the mode in the background, its output appended to
# OUT/<PREFIX>.<mode>.out. The GPU modes take OUT/<PREFIX>.lock, so a
# second one exits while the first is alive (a lock left by an earlier
# boot is taken over). Every stage's outcome goes to OUT/<PREFIX>.stages.log.
#   scripts/verdaEleven.sh bootstrap --detach
#   PRICE=1.88 scripts/verdaEleven.sh bench
# Environment: PREFIX (eleven), OUT ($HOME/runs/<PREFIX>), DRIVER (paired
# or async), NS ("64 128"), PER_YEAR (36), YEARS (3), UNTIL (365·YEARS;
# `run` writes the value it uses
# to OUT/UNTIL_<PREFIX>, which a later run or resume without UNTIL takes, so
# a run extended by hand stays extended after an eviction), KEEP (1000), LEVELS (bl60),
# OCEAN ('{"everySteps":8}'), STRATOSPHERE (1), SEGMENT_MINUTES (none:
# pairedSpinup.sh's segment length in minutes of wall time, each segment
# ending with a snapshot), RADIATION and MOIST (none: spinup.mjs's option
# objects, as JSON), for DRIVER=async CYCLES
# (9), COUPLED_YEARS (10), RECORD_YEARS (1), OCEAN_YEARS (100),
# FINAL_YEARS (0), OCEAN_ONLY (OCEAN), RESTORE (30), OCEAN_KEEP (2),
# SNAPSHOT_DAYS (30) and KEEP_FORCING (0) (asyncSpinup.sh's defaults),
# PRICE (1.85 $/h), CASES
# ("128:2 64:5"), JOBS (the core count), GPU_WAIT (600), EXPECT_VENDOR,
# MARGIN_GB (5, free space kept beyond the states), FORCE_DISK (1 runs
# without the room), PAIRED (scripts/pairedSpinup.sh, or
# scripts/asyncPair.sh for DRIVER=async), ASYNC (asyncPair.sh's) and
# SKIP_GPU_CHECK=1 (for tests). `run` writes the run's settings (DRIVER,
# NS, PER_YEAR, KEEP, LEVELS, OCEAN, STRATOSPHERE, SEGMENT_MINUTES,
# RADIATION, MOIST and for async the schedule's) to OUT/SETTINGS_<PREFIX>,
# one NAME=value a line, and every mode takes from there those left
# unset, so that the resume at boot, given OUT and PREFIX alone, runs on
# the run that was started.
SELF=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
cd "$(dirname "$0")/.."
PREFIX=${PREFIX:-eleven}
OUT=${OUT:-$HOME/runs/$PREFIX}
SETTINGS=$OUT/SETTINGS_$PREFIX
COMMON="DRIVER NS PER_YEAR KEEP LEVELS OCEAN STRATOSPHERE SEGMENT_MINUTES RADIATION MOIST"
SCHEDULE="CYCLES COUPLED_YEARS RECORD_YEARS OCEAN_YEARS FINAL_YEARS OCEAN_ONLY RESTORE OCEAN_KEEP SNAPSHOT_DAYS KEEP_FORCING"
for name in $COMMON $SCHEDULE; do
  [ -z "${!name}" ] && grep -q "^$name=" "$SETTINGS" 2>/dev/null && printf -v "$name" '%s' "$(sed -n "s/^$name=//p" "$SETTINGS" | tail -1)"
done
DRIVER=${DRIVER:-paired}
NS=${NS:-"64 128"} PER_YEAR=${PER_YEAR:-36} YEARS=${YEARS:-3} KEEP=${KEEP:-1000}
UNTIL=${UNTIL:-$(cat "$OUT/UNTIL_$PREFIX" 2>/dev/null || echo $((365 * YEARS)))}
LEVELS=${LEVELS:-bl60} OCEAN=${OCEAN:-'{"everySteps":8}'} STRATOSPHERE=${STRATOSPHERE:-1} SEGMENT_MINUTES=${SEGMENT_MINUTES:-}
for name in RADIATION MOIST; do if [ -n "${!name}" ]; then export "$name"; else unset "$name"; fi; done
CYCLES=${CYCLES:-9} COUPLED_YEARS=${COUPLED_YEARS:-10} RECORD_YEARS=${RECORD_YEARS:-1} OCEAN_YEARS=${OCEAN_YEARS:-100} FINAL_YEARS=${FINAL_YEARS:-0}
OCEAN_ONLY=${OCEAN_ONLY:-$OCEAN} RESTORE=${RESTORE:-30} OCEAN_KEEP=${OCEAN_KEEP:-2} SNAPSHOT_DAYS=${SNAPSHOT_DAYS:-30} KEEP_FORCING=${KEEP_FORCING:-0}
PRICE=${PRICE:-1.85} CASES=${CASES:-"128:2 64:5"} GPU_WAIT=${GPU_WAIT:-600} MARGIN_GB=${MARGIN_GB:-5}
JOBS=${JOBS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)}
[ "$(uname)" = Linux ] && EXPECT_VENDOR=${EXPECT_VENDOR-nvidia}
case "$DRIVER" in
  paired) PAIRED=${PAIRED:-scripts/pairedSpinup.sh} ;;
  async) PAIRED=${PAIRED:-scripts/asyncPair.sh} ;;
  *) echo "DRIVER is paired or async, not $DRIVER" >&2; exit 2 ;;
esac
LOCK=$OUT/$PREFIX.lock
SUDO=; [ "$(id -u)" -ne 0 ] && SUDO=sudo

mode= detach=
for arg in "$@"; do case "$arg" in --detach) detach=1 ;; *) mode=$arg ;; esac; done
case "$mode" in setup|autostart|suite|bench|bootstrap|run|resume|status) ;; *) echo "usage: $0 setup|autostart|suite|bench|bootstrap|run|resume|status [--detach]" >&2; exit 2 ;; esac
mkdir -p "$OUT"

say() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a "$OUT/$PREFIX.stages.log"; }

if [ -n "$detach" ]; then
  out=$OUT/$PREFIX.$mode.out
  if command -v setsid > /dev/null; then setsid nohup bash "$SELF" "$mode" >> "$out" 2>&1 < /dev/null &
  else nohup bash "$SELF" "$mode" >> "$out" 2>&1 < /dev/null &
  fi
  echo "$mode running in the background (pid $!), output in $out"
  exit 0
fi

boot() { cat /proc/sys/kernel/random/boot_id 2>/dev/null || sysctl -n kern.boottime 2>/dev/null | tr -dc '0-9' | cut -c1-12; }
lock() {
  local pid owner_boot owner_mode
  if ! mkdir "$LOCK" 2>/dev/null; then
    read -r pid owner_boot owner_mode < "$LOCK/owner" 2>/dev/null || { sleep 2; read -r pid owner_boot owner_mode < "$LOCK/owner" 2>/dev/null; }
    if [ -n "$pid" ] && [ "$owner_boot" = "$(boot)" ] && kill -0 "$pid" 2>/dev/null; then
      say "$mode not started: $owner_mode is running (pid $pid)"; exit 1
    fi
    say "taking over the lock of ${owner_mode:-an unknown mode} (pid ${pid:-?}, gone)"
    rm -rf "$LOCK"; mkdir "$LOCK" || exit 1
  fi
  echo "$$ $(boot) $mode" > "$LOCK/owner"
  trap 'rm -rf "$LOCK"' EXIT
}

apt_install() { $SUDO env DEBIAN_FRONTEND=noninteractive NEEDRESTART_MODE=a apt-get -o DPkg::Lock::Timeout=600 install -y -qq "$@" > /dev/null; }
apt_updated=
apt_ready() {
  command -v apt-get > /dev/null || { say "setup: $1 is missing and there is no apt-get"; return 1; }
  [ -n "$apt_updated" ] || { $SUDO apt-get -o DPkg::Lock::Timeout=600 update -qq > /dev/null && apt_updated=1; }
}
node_major() { node -e 'console.log(process.versions.node.split(".")[0])' 2>/dev/null || echo 0; }

install() {
  if [ "$(node_major)" -lt 22 ]; then
    apt_ready "node 22 or later" || return 1
    command -v curl > /dev/null || apt_install curl ca-certificates || return 1
    curl -fsSL https://deb.nodesource.com/setup_22.x | ${SUDO:+$SUDO -E} bash - > /dev/null && apt_install nodejs || { say "setup: installing node failed"; return 1; }
  fi
  command -v rsync > /dev/null || { apt_ready rsync && apt_install rsync; } || return 1
  if [ "$(uname)" = Linux ] && ! command -v vulkaninfo > /dev/null; then apt_ready vulkaninfo && apt_install libvulkan1 vulkan-tools || return 1; fi
  if [ ! -f node_modules/webgpu/package.json ]; then
    npm ci --no-audit --no-fund > /dev/null 2>&1 || { apt_ready "a compiler for npm ci" && apt_install build-essential python3 && npm ci --no-audit --no-fund > /dev/null; } || { say "setup: npm ci failed"; return 1; }
  fi
  say "setup: node $(node --version), $(rsync --version 2>/dev/null | head -1 | cut -c1-40), webgpu $(node -p 'require("./node_modules/webgpu/package.json").version' 2>/dev/null)"
}

autostart() {
  local dir=${SYSTEMD_DIR:-/etc/systemd/system} ctl=${SYSTEMCTL:-systemctl} unit=gcm-$PREFIX-resume.service as=
  if [ "$(uname)" != Linux ] && [ -z "$SYSTEMD_DIR" ]; then say "autostart: no systemd here; nothing installed"; return 0; fi
  command -v "$ctl" > /dev/null || { say "autostart: $ctl is missing; the startup script alone resumes the run"; return 0; }
  [ -w "$dir" ] || as=$SUDO
  printf '%s\n' '[Unit]' "Description=Resume the spin-up $PREFIX at boot" 'Wants=network-online.target' 'After=network-online.target' '' \
    '[Service]' 'Type=simple' "Environment=HOME=$HOME OUT=$OUT PREFIX=$PREFIX" "WorkingDirectory=$PWD" "ExecStart=/bin/bash $SELF resume" \
    "StandardOutput=append:$OUT/$PREFIX.resume.out" "StandardError=append:$OUT/$PREFIX.resume.out" 'TimeoutStopSec=300' '' \
    '[Install]' 'WantedBy=multi-user.target' | $as tee "$dir/$unit" > /dev/null || { say "autostart: cannot write $dir/$unit"; return 1; }
  $as "$ctl" daemon-reload && $as "$ctl" enable "$unit" > /dev/null 2>&1 || { say "autostart: $ctl could not enable $unit"; return 1; }
  say "autostart: $unit enabled ($dir), resume runs at every boot"
}

gpu_check() {
  [ "$SKIP_GPU_CHECK" = 1 ] && return 0
  local file=$OUT/${PREFIX}_gpu.txt ok=0
  {
    echo "== $(date -u '+%Y-%m-%dT%H:%M:%SZ') on $(hostname)"
    if command -v nvidia-smi > /dev/null; then nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader || ok=1; fi
    if command -v vulkaninfo > /dev/null; then
      vulkan=$(vulkaninfo --summary 2>/dev/null | grep -E 'deviceName|deviceType|driverName|apiVersion')
      echo "$vulkan"
      [ "$EXPECT_VENDOR" = nvidia ] && ! echo "$vulkan" | grep -qi nvidia && { echo "vulkaninfo shows no NVIDIA device"; ok=1; }
    fi
    EXPECT_VENDOR=$EXPECT_VENDOR node --input-type=module -e "
      const { create } = await import('webgpu');
      const adapter = await create([]).requestAdapter();
      if (!adapter) { console.log('no WebGPU adapter'); process.exit(1); }
      const i = adapter.info ?? {};
      console.log('adapter:', [i.vendor, i.architecture, i.device, i.description].filter(Boolean).join(' | '));
      const want = process.env.EXPECT_VENDOR;
      if (want && !String(i.vendor).toLowerCase().includes(want)) { console.log('the adapter is not ' + want); process.exit(1); }
      await adapter.requestDevice();
      console.log('device ok');
      process.exit(0);" || ok=1
  } > "$file" 2>&1
  [ $ok -eq 0 ] && say "gpu: $(grep -E '^adapter:' "$file" | head -1)" || say "gpu: check failed, see $file: $(tail -2 "$file" | tr '\n' ' ')"
  return $ok
}

suite() {
  local report=$OUT/${PREFIX}_suite.txt
  JOBS=$JOBS scripts/suiteReport.sh "$report" > /dev/null
  say "suite: $(head -1 "$report")"
}

bench() {
  local report=$OUT/${PREFIX}_benchmark.txt
  if OUT=$OUT/${PREFIX}_bench REPORT=$report CASES=$CASES LEVELS=$LEVELS OCEAN=$OCEAN STRATOSPHERE=$STRATOSPHERE PRICE=$PRICE PER_YEAR=$PER_YEAR UNTIL=$UNTIL node scripts/verdaBenchmark.mjs; then
    say "bench: $(grep '^RESULT' "$report")"
  else
    say "bench: failed: $(tail -1 "$report" 2>/dev/null)"; return 1
  fi
}

newest() { ls "$OUT" | grep -E "^$PREFIX$1_day[0-9]+\.bin$" | sort | tail -1; }
dayOf() { local file; file=$(newest "$1"); [ -n "$file" ] && echo "$file" | sed -E 's/.*_day0*([0-9]+)\.bin$/\1/' || echo 0; }
remaining() { awk -v d="$1" -v q="$PER_YEAR" -v u="$UNTIL" 'BEGIN { n = 0; for (k = 1; ; k++) { t = int(k * 365 / q + 0.5); if (t > d) n++; if (t >= u) break } print n }'; }
stateBytes() {
  local file; file=$(newest "$1")
  if [ -n "$file" ]; then wc -c < "$OUT/$file" | tr -d ' '; return; fi
  node -e 'const j = JSON.parse(require("fs").readFileSync(process.argv[1], "utf8")); const m = j.measured.find((c) => c.N === Number(process.argv[2])); if (!m) process.exit(1); console.log(m.bytes)' "$OUT/${PREFIX}_benchmark.json" "$1" 2>/dev/null \
    || awk -v n="$1" 'BEGIN { printf "%d\n", 425e6 * (n / 128) ^ 2 }'
}

ledgerOf() { sed -n "s/^$2=//p" "$OUT/$PREFIX$1_cycles.txt" 2>/dev/null | tail -1; }
asyncStart() { local start; start=$(ledgerOf "$1" start); echo "${start:-$(dayOf "$1")}"; }
asyncLast() { local last=$CYCLES converged; converged=$(ledgerOf "$1" converged); [ -n "$converged" ] && [ "$converged" -lt "$last" ] && last=$converged; echo "$last"; }
asyncEnd() { echo $(($(asyncStart "$1") + 365 * ($(asyncLast "$1") * COUPLED_YEARS + FINAL_YEARS))); }
asyncUntil() { local n end until=0; for n in $NS; do end=$(asyncEnd "$n"); [ "$end" -gt "$until" ] && until=$end; done; echo "$until"; }
snapshotsBetween() { awk -v a="$1" -v b="$2" -v q="$PER_YEAR" 'BEGIN { if (a >= b) { print 0; exit } n = 1; for (k = 1; ; k++) { t = int(k * 365 / q + 0.5); if (t >= b) break; if (t > a) n++ } print n }'; }
forcingBytes() {
  local file; file=$(ls "$OUT/$PREFIX$1"_c*_forcing/forcing-*.bin 2>/dev/null | head -1)
  if [ -n "$file" ]; then wc -c < "$file" | tr -d ' '; else awk -v n="$1" 'BEGIN { printf "%d\n", 2.3e6 * (n / 64) ^ 2 }'; fi
}

room() {
  local total=0 need=0 line= n day left bytes free_kb start end snaps years forcing all present at
  for n in $NS; do
    bytes=$(stateBytes "$n")
    if [ "$DRIVER" = async ]; then
      start=$(asyncStart "$n") end=$(asyncEnd "$n")
      snaps=$(snapshotsBetween "$start" "$end"); [ "$snaps" -gt "$KEEP" ] && snaps=$KEEP
      years=$((CYCLES + OCEAN_KEEP))
      forcing=$(awk -v d="$(forcingBytes "$n")" -v r="$RECORD_YEARS" -v c="$([ "$KEEP_FORCING" = 1 ] && echo "$CYCLES" || echo 1)" 'BEGIN { print d * 365 * r * c }')
      all=$(awk -v s="$snaps" -v y="$years" -v b="$bytes" -v f="$forcing" 'BEGIN { print (s + y) * b + f }')
      present=$(du -skc "$OUT/$PREFIX${n}_"* 2>/dev/null | tail -1 | cut -f1)
      total=$(awk -v t="$total" -v a="$all" 'BEGIN { print t + a }')
      need=$(awk -v t="$need" -v a="$all" -v p="${present:-0}" 'BEGIN { r = a - p * 1024; print t + (r > 0 ? r : 0) }')
      line="$line N=$n from day $start to $end, $snaps snapshots and $years ocean-only year files × $(awk -v b="$bytes" 'BEGIN { printf "%.0f", b / 1e6 }') MB and $(awk -v f="$forcing" 'BEGIN { printf "%.1f", f / 1e9 }') GB of forcing;"
    else
      day=$(dayOf "$n") left=$(remaining "$(dayOf "$n")")
      total=$(awk -v t="$total" -v a="$(remaining 0)" -v b="$bytes" 'BEGIN { print t + a * b }')
      need=$(awk -v t="$need" -v a="$left" -v b="$bytes" 'BEGIN { print t + a * b }')
      line="$line N=$n day $day, $left segments left × $(awk -v b="$bytes" 'BEGIN { printf "%.0f", b / 1e6 }') MB;"
    fi
  done
  at=$([ "$DRIVER" = async ] && echo "the files take" || echo "the states take")
  free_kb=$(df -Pk "$OUT" | awk 'NR == 2 { print $4 }')
  say "disk:$line $at $(awk -v t="$total" 'BEGIN { printf "%.1f", t / 1e9 }') GB at day $UNTIL, $(awk -v t="$need" 'BEGIN { printf "%.1f", t / 1e9 }') GB of it still to come, $(awk -v f="$free_kb" 'BEGIN { printf "%.1f", f * 1024 / 1e9 }') GB free"
  awk -v f="$free_kb" -v t="$need" -v m="$MARGIN_GB" 'BEGIN { exit !(f * 1024 >= t + m * 1e9) }'
}

run() {
  if [ -f "$OUT/STOP_$PREFIX" ]; then say "run: STOP_$PREFIX is in $OUT; remove it to run"; return 1; fi
  if ! room && [ "$FORCE_DISK" != 1 ]; then say "run: not started, the disk lacks room for the states plus $MARGIN_GB GB (FORCE_DISK=1 runs anyway)"; return 1; fi
  [ -f "$OUT/STARTED_$PREFIX" ] || echo "$(date -u '+%Y-%m-%dT%H:%M:%SZ') commit $(git rev-parse --short HEAD 2>/dev/null)" > "$OUT/STARTED_$PREFIX"
  echo "$UNTIL" > "$OUT/UNTIL_$PREFIX"
  local name schedule=()
  [ "$DRIVER" = async ] && { schedule=(DRIVER=async); for name in $SCHEDULE; do schedule+=("$name=${!name}"); done; }
  { for name in $COMMON; do echo "$name=${!name}"; done; [ "$DRIVER" = async ] && for name in $SCHEDULE; do echo "$name=${!name}"; done; } > "$SETTINGS"
  rm -f "$OUT/ENDED_$PREFIX"
  local stops
  stops=$(grep -c 'stopped (' "$OUT/$PREFIX.log" 2>/dev/null); stops=${stops:-0}
  say "run: NS=\"$NS\" PREFIX=$PREFIX LEVELS=$LEVELS PER_YEAR=$PER_YEAR KEEP=$KEEP OCEAN=$OCEAN STRATOSPHERE=$STRATOSPHERE SEGMENT_MINUTES=$SEGMENT_MINUTES${RADIATION:+ RADIATION=$RADIATION}${MOIST:+ MOIST=$MOIST} UNTIL=$UNTIL OUT=$OUT${schedule[*]:+ ${schedule[*]}} $PAIRED"
  env -u FROM -u OCEAN_FROM -u LAND_FROM -u ICE_FROM -u RECORD -u STOP_AFTER_STEPS -u N SYNC_CMD='sync "$1"' \
    NS="$NS" PREFIX="$PREFIX" LEVELS="$LEVELS" PER_YEAR="$PER_YEAR" KEEP="$KEEP" OCEAN="$OCEAN" STRATOSPHERE="$STRATOSPHERE" SEGMENT_MINUTES="$SEGMENT_MINUTES" UNTIL="$UNTIL" OUT="$OUT" "${schedule[@]}" bash "$PAIRED"
  local code=$? reason
  reason=$(grep 'stopped (' "$OUT/$PREFIX.log" 2>/dev/null | tail -n +"$((stops + 1))" | tail -1)
  if [ -z "$reason" ]; then say "run: interrupted (exit $code) before $PAIRED logged a stop; not marked ended, so resume continues it"; return 1; fi
  echo "$(date -u '+%Y-%m-%dT%H:%M:%SZ') exit $code: $reason" > "$OUT/ENDED_$PREFIX"
  say "run: ended (exit $code): $reason"
}

resume() {
  if [ ! -f "$OUT/STARTED_$PREFIX" ]; then say "resume: the run was never started in $OUT; nothing to do"; return 0; fi
  if [ -f "$OUT/ENDED_$PREFIX" ]; then say "resume: the run ended ($(cat "$OUT/ENDED_$PREFIX")); nothing to do"; return 0; fi
  if [ -f "$OUT/STOP_$PREFIX" ]; then say "resume: STOP_$PREFIX is in $OUT; nothing to do"; return 0; fi
  local waited=0
  until gpu_check; do
    [ "$waited" -ge "$GPU_WAIT" ] && { say "resume: no usable GPU after $GPU_WAIT s"; return 1; }
    sleep 15; waited=$((waited + 15))
  done
  run
}

asyncPhase() {
  local n=$1 day start last finished c end from tag year
  day=$(dayOf "$n") start=$(asyncStart "$n") last=$(asyncLast "$n") finished=$(ledgerOf "$n" done); finished=${finished:-0}
  if [ ! -f "$OUT/$PREFIX${n}_cycles.txt" ]; then echo "not started"; return; fi
  if [ "$finished" -lt "$last" ]; then
    c=$((finished + 1)) end=$((start + 365 * (finished + 1) * COUPLED_YEARS)); from=$((end - 365 * RECORD_YEARS))
    if [ "$day" -lt "$end" ]; then
      echo "cycle $c of $last: coupled, day $day of $end, recording the forcing $([ "$day" -ge "$from" ] && echo since || echo from) day $from"
    else
      tag=$(printf '%s%s_c%02d' "$PREFIX" "$n" "$c")
      year=$(ls "$OUT" | grep -E "^${tag}_year[0-9]+(_day[0-9]+(_step[0-9]+)?)?\.bin$" | sort | tail -1)
      if [ -n "$year" ]; then echo "cycle $c of $last: ocean-only, $(echo "$year" | sed -E 's/.*_year0*([0-9]+).*/\1/') of $OCEAN_YEARS years done ($year)"
      else echo "cycle $c of $last: ocean-only, starting from day $end"
      fi
    fi
  elif [ "$day" -lt "$(asyncEnd "$n")" ]; then echo "after the last cycle: coupled, day $day of $(asyncEnd "$n")"
  else echo "done at coupled day $day"
  fi
}

status() {
  local n file
  for n in $NS; do
    file=$(newest "$n")
    if [ "$DRIVER" = async ]; then
      echo "N=$n: ${file:-no state yet}; $(asyncPhase "$n"); the schedule ends at coupled day $(asyncEnd "$n")"
      [ -s "$OUT/$PREFIX${n}_cycles.log" ] && echo "  $(tail -1 "$OUT/$PREFIX${n}_cycles.log")"
    else
      echo "N=$n: ${file:-no state yet} (day $(dayOf "$n") of $UNTIL, $(remaining "$(dayOf "$n")") segments left)"
    fi
  done
  for marker in STARTED ENDED STOP; do [ -f "$OUT/${marker}_$PREFIX" ] && echo "${marker}_$PREFIX: $(cat "$OUT/${marker}_$PREFIX")"; done
  [ -f "$LOCK/owner" ] && echo "lock: $(cat "$LOCK/owner")"
  [ -f "$OUT/$PREFIX.log" ] && { echo "-- $PREFIX.log"; tail -n 3 "$OUT/$PREFIX.log"; }
  df -h "$OUT" | tail -1
}

[ "$DRIVER" = async ] && UNTIL=$(asyncUntil)
case "$mode" in
  setup) lock; install && { autostart; gpu_check; } ;;
  autostart) autostart ;;
  suite) lock; suite ;;
  bench) lock; bench ;;
  bootstrap) lock; install && { autostart; gpu_check; } || exit 1; suite; bench ;;
  run) lock; gpu_check && run ;;
  resume) [ -f "$OUT/STARTED_$PREFIX" ] && [ ! -f "$OUT/ENDED_$PREFIX" ] && [ ! -f "$OUT/STOP_$PREFIX" ] && lock; resume ;;
  status) status ;;
esac
