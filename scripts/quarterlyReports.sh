#!/bin/bash
# Quarterly reports of a paired spin-up whose resolutions may run on
# different machines. For each snapshot day round(k·365/PER_YEAR) up to
# UNTIL with a state <PREFIX><N>_dayDDDD.bin for every N in NS, writes
# <REPORTS>/dayDDDD/<PREFIX>_dayDDDD.md: scripts/compareStates.mjs side by
# side over the quarter, scripts/woaSstBias.mjs merged by box,
# scripts/cloudRegimes.mjs per resolution, scripts/rounds.py and a
# top-of-the-atmosphere table by quarter over the logs up to that day, and
# the five figures of scripts/figures/snapshot.sh per resolution, which sit
# beside it; the tools' raw outputs stay under <REPORTS>/dayDDDD/parts, and a
# part that failed is marked <name>.failed there and not retried (delete the
# marker to retry). The figure dumps use the GPU, so they wait while another
# node script or a page server on port 8000 has it; a report is rewritten
# whenever a part or figure lands. Loops every INTERVAL seconds until every
# report up to UNTIL is complete; ONCE=1 makes a single pass.
#
#   OUT=runs/thirteen DIR_128=runs/verda-thirteen UNTIL=1825 scripts/quarterlyReports.sh
#
# Environment: OUT (the run directory; REPORTS defaults to OUT/reports),
# DIR_<N> (the directory holding resolution N's states and log, OUT by
# default), PREFIX (thirteen), NS ("64 128"), PER_YEAR (4), UNTIL (1825),
# INTERVAL (300), NOTE (a line for each report's header), NODE, PYTHON.
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/.." && pwd)
OUT=${OUT:?OUT: the run directory}; PREFIX=${PREFIX:-thirteen}; NS=${NS:-"64 128"}; PER_YEAR=${PER_YEAR:-4}
UNTIL=${UNTIL:-1825}; INTERVAL=${INTERVAL:-300}; NOTE=${NOTE:-}; NODE=${NODE:-node}; PYTHON=${PYTHON:-python3}
REPORTS=${REPORTS:-$OUT/reports}; LOG=$REPORTS/reports.log
SNAPSHOT=${SNAPSHOT:-$HERE/figures/snapshot.sh}; FIGURES="stateMaps oceanPanels eqsection eqpanels mlmdeck"
mkdir -p "$REPORTS" || exit 1
cd "$REPO" || exit 1

log() { echo "$(date '+%F %T') $*" >> "$LOG"; }
dir_of() { local v=DIR_$1; echo "${!v:-$OUT}"; }
state_of() { printf '%s/%s%s_day%04d.bin' "$(dir_of "$1")" "$PREFIX" "$1" "$2"; }
tag_of() { printf '%s_day%04d' "$PREFIX" "$1"; }
day_dir() { printf '%s/day%04d' "$REPORTS" "$1"; }
parts_of() { echo "$(day_dir "$1")/parts"; }
qday() { echo $(( ($1 * 365 * 2 + PER_YEAR) / (2 * PER_YEAR) )); }
calendar() { "$PYTHON" -c "from datetime import date, timedelta; print((date(2001, 3, 20) + timedelta(days=$1)).strftime('%b %-d'))"; }
gpu_user() {
  local p; p=$(ps -axo command | grep -E '^([^ ]*/)?node (--test|([^ ]*/)?scripts/|[^ ]*/scratchpad/)' | grep -vE '^([^ ]*/)?node ([^ ]*/)?scripts/(compareStates|woaSstBias|cloudRegimes)\.mjs' | head -n 1 | cut -c1-80)
  [ -n "$p" ] && { echo "$p"; return; }
  lsof -nP -iTCP:8000 -sTCP:LISTEN > /dev/null 2>&1 && echo "a page server on port 8000"
}

run_part() { # <day> <name> <command...>: the command's stdout becomes <name> in the day's parts
  local tag=$(tag_of "$1") parts=$(parts_of "$1") name=$2; shift 2
  [ -s "$parts/$name" ] || [ -e "$parts/$name.failed" ] && return 0
  log "$tag: $name"
  if "$@" > "$parts/$name.tmp" 2> "$parts/$name.err"; then mv "$parts/$name.tmp" "$parts/$name"
  else rm -f "$parts/$name.tmp"; touch "$parts/$name.failed"; log "$tag: $name failed: $(tail -n 1 "$parts/$name.err" | cut -c1-200)"; fi
}

log_to_day() { # <N> <day> <dir>: a copy of resolution N's log up to that day
  local dst=$3/$PREFIX$1.log
  awk -v d="$2" '/^day [0-9]+ / { if ($2 + 0 > d) exit } { print }' "$(dir_of "$1")/$PREFIX$1.log" > "$dst"; echo "$dst"
}

toa_table() { # <day> <log>...: quarter means of the daily lines up to that day, beside CERES
  "$PYTHON" - "$PER_YEAR" "$@" <<'PY'
import re, sys
per_year, D, logs = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3:]
pat = re.compile(r'^day (\d+) .*?Ts (-?[\d.]+) °C, ASR ([\d.]+)(?: \(atmosphere ([\d.]+)\))? OLR ([\d.]+).*?albedo ([\d.]+)(?:, SWCRE (-?[\d.]+) LWCRE (-?[\d.]+))?(?:, clear-sky reflectance ([\d.]+))?')
qend = lambda k: int(k * 365 / per_year + 0.5)
f = lambda x, d=1: '' if x is None else f'{x:.{d}f}'
cols = ['Ts (°C)', 'ASR', 'OLR', 'ASR − OLR', 'albedo', 'SWCRE', 'LWCRE', 'clear-sky OLR', 'absorbed in the air', 'clear-sky reflectance']
print('| Run | Window | ' + ' | '.join(cols) + ' |'); print('|---|---|' + '---|' * len(cols))
for path in logs:
    name, rows = path.split('/')[-1].rsplit('.', 1)[0], {}
    for line in open(path):
        m = pat.match(line)
        if m: rows[int(m.group(1))] = [None if x is None else float(x) for x in m.groups()[1:]]
    k = 1
    while qend(k) <= D:
        a, b = qend(k - 1) + 1, qend(k); win = [rows[d] for d in range(a, b + 1) if d in rows]
        if len(win) >= (b - a) - 2:
            mean = lambda i: (lambda v: sum(v) / len(v) if v else None)([r[i] for r in win if r[i] is not None])
            ts, asr, atm, olr, alb, sw, lw, csr = [mean(i) for i in range(8)]
            net = None if asr is None or olr is None else asr - olr; clear = None if olr is None or lw is None else olr + lw
            print(f'| {name} | Q{k} days {a}–{b} | {f(ts, 2)} | {f(asr)} | {f(olr)} | {f(net)} | {f(alb, 3)} | {f(sw)} | {f(lw)} | {f(clear)} | {f(atm)} | {f(csr, 3)} |')
        k += 1
print('| CERES EBAF, annual mean | | 14.85 | 241 | 240 | +0.8 | 0.29 | −46 | +27 | 267 | about 79 | about 0.15 |')
print()
print("W/m² except as marked; the log's daily values are means over the day's steps. The model's orbit is circular, 340.5 W/m² incoming in every quarter, so its absorbed sunlight compares with the annual 241 in every quarter; the observed outgoing longwave runs about 241, 243, 239 and 237 over the four quarters from the March equinox.")
PY
}

make_parts() { # <day> <k>
  local D=$1 k=$2 parts=$(parts_of "$1") N states=() logs=()
  local window=$((D - $(qday $((k - 1)))))
  mkdir -p "$parts"
  for N in $NS; do states+=("$(state_of "$N" "$D")"); done
  run_part "$D" compare.md "$NODE" scripts/compareStates.mjs "${states[@]}" --window "$window"
  for N in $NS; do
    run_part "$D" "woa$N.md" "$NODE" scripts/woaSstBias.mjs "$(state_of "$N" "$D")" --json "$parts/woa$N.json"
    run_part "$D" "cloud$N.txt" "$NODE" scripts/cloudRegimes.mjs "$(state_of "$N" "$D")"
  done
  if [ ! -s "$parts/rounds.txt" ] || [ ! -s "$parts/toa.md" ]; then
    for N in $NS; do logs+=("$(log_to_day "$N" "$D" "$parts")"); done
    PER_YEAR=$PER_YEAR run_part "$D" rounds.txt "$PYTHON" scripts/rounds.py "${logs[@]}"
    run_part "$D" toa.md toa_table "$D" "${logs[@]}"
  fi
}

figure_png() { printf '%s/%s%s_%s_day%04d.png' "$(day_dir "$3")" "$PREFIX" "$1" "$2" "$3"; }
figures_missing() { # <N> <day>: 0 when a figure is still to be dumped
  local F; for F in $FIGURES; do [ -s "$(figure_png "$1" "$F" "$2")" ] || return 0; done; return 1
}

make_figures() { # <day>: 1 when a dump waits for the GPU
  local D=$1 parts=$(parts_of "$1") N F busy tmp start png aborts
  for N in $NS; do
    figures_missing "$N" "$D" || continue
    [ -e "$parts/figures$N.failed" ] && continue
    busy=$(gpu_user); [ -n "$busy" ] && { log "$(tag_of "$D"): figures N=$N wait for the GPU ($busy)"; return 1; }
    tmp=$(mktemp -d "${TMPDIR:-/tmp}/quarterly-figures.XXXXXX"); start=$(date +%s)
    log "$(tag_of "$D"): figures N=$N"
    bash "$SNAPSHOT" "$(state_of "$N" "$D")" "$tmp" "$PREFIX N=$N, day $D" > "$parts/figures$N.log" 2>&1
    for F in $FIGURES; do
      png=$(ls "$tmp"/*_"$F"_day*.png 2> /dev/null | head -n 1)
      [ -n "$png" ] && mv "$png" "$(figure_png "$N" "$F" "$D")"
      cat "$tmp"/*_"$F"_day*.log >> "$parts/figures$N.log" 2> /dev/null
    done
    rm -rf "$tmp"
    aborts=$(/usr/bin/log show --last "$(( ($(date +%s) - start) / 60 + 2 ))m" --predicate 'eventMessage CONTAINS "aborted"' --style compact 2> /dev/null | grep 'aborted' | grep -vc 'com.apple.log')
    log "$(tag_of "$D"): figures N=$N done in $(( $(date +%s) - start )) s; 'aborted' lines in the system log over the dump: $aborts"
    figures_missing "$N" "$D" && { touch "$parts/figures$N.failed"; log "$(tag_of "$D"): figures N=$N incomplete: $(grep -E 'failed' "$parts/figures$N.log" | head -n 2 | tr '\n' ' ' | cut -c1-200)"; }
  done
  return 0
}

figures_done() { local N; for N in $NS; do figures_missing "$N" "$1" && [ ! -e "$(parts_of "$1")/figures$N.failed" ] && return 1; done; return 0; }

woa_table() { # <parts>
  local args=() N
  for N in $NS; do [ -s "$1/woa$N.json" ] && args+=("$N" "$1/woa$N.json"); done
  [ ${#args[@]} -eq 0 ] && { echo "_missing_"; return; }
  "$NODE" -e '
    const fs = require("node:fs"), runs = [];
    for (let a = 1; a + 1 < process.argv.length; a += 2) runs.push({ n: process.argv[a], rows: JSON.parse(fs.readFileSync(process.argv[a + 1], "utf8")).rows });
    const f = (x, signed) => (x == null || Number.isNaN(x) ? "" : (signed && x >= 0 ? "+" : "") + x.toFixed(2));
    const pair = (at, key, signed) => at.map((x) => f(x[key], signed)).join(" / ");
    const head = ["Box", "WOA23 (°C)", "Model (°C)", "Bias (K)", "RMS (K)", "Iced"];
    if (runs.length > 1) head.push(`Bias N=${runs[runs.length - 1].n} − N=${runs[0].n} (K)`);
    console.log(`Each cell ${runs.map((r) => `N=${r.n}`).join(" / ")}.`); console.log();
    console.log(`| ${head.join(" | ")} |`); console.log(`|${head.map(() => "---").join("|")}|`);
    for (const row of runs[0].rows) {
      const at = runs.map((r) => r.rows.find((x) => x.box === row.box) ?? {});
      const cells = [row.box, pair(at, "woa"), pair(at, "model"), pair(at, "bias", true), pair(at, "rms"), row.cells == null ? "" : at.map((x) => x.iced == null ? "" : `${Math.round(x.iced * 100)} %`).join(" / ")];
      if (runs.length > 1) cells.push(at[at.length - 1].bias == null || at[0].bias == null ? "" : f(at[at.length - 1].bias - at[0].bias, true));
      console.log(`| ${cells.join(" | ")} |`);
    }' "${args[@]}"
}

woa_dates() { # <parts>: what each resolution's table compares, over which days, and the notes on it
  local N any=0
  for N in $NS; do [ -s "$1/woa$N.md" ] && { sed -n "1s/^- [^:]*: /- N=$N: /p" "$1/woa$N.md"; any=1; }; done
  [ "$any" = 1 ] || return 0
  echo; for N in $NS; do [ -s "$1/woa$N.md" ] && sed -n '3s/^- //p' "$1/woa$N.md"; done | awk '!seen[$0]++'; echo
}

section() { if [ -s "$1" ]; then cat "$1"; else echo "_$(basename "$1") is missing; see reports.log._"; fi; echo; }
code() { if [ -s "$1" ]; then echo '```text'; cat "$1"; echo '```'; else echo "_$(basename "$1") is missing; see reports.log._"; fi; echo; }

assemble() { # <day> <k>
  local D=$1 k=$2 parts=$(parts_of "$1") md=$(day_dir "$1")/$(tag_of "$1").md N F png
  local year=$(( (k - 1) / PER_YEAR + 1 )) quarter=$(( (k - 1) % PER_YEAR + 1 )) prev=$(qday $((k - 1)))
  {
    echo "# $PREFIX, day $D: year $year, quarter $quarter ($(calendar "$D"))"; echo
    [ -n "$NOTE" ] && { echo "$NOTE"; echo; }
    echo "States:$(for N in $NS; do printf ' %s' "$(state_of "$N" "$D")"; done). The log means cover days $((prev + 1))–$D. Written $(date '+%F %H:%M')."; echo
    echo "## Side by side"; echo; section "$parts/compare.md"
    echo "## Sea surface temperature against WOA23"; echo
    woa_dates "$parts"
    woa_table "$parts"; echo
    echo "## Cloud by regime"; echo
    for N in $NS; do echo "### N=$N"; echo; code "$parts/cloud$N.txt"; done
    echo "## Top of the atmosphere by quarter"; echo; section "$parts/toa.md"
    echo "## Rounds"; echo; code "$parts/rounds.txt"
    echo "## Figures"; echo
    for N in $NS; do
      echo "### N=$N"; echo
      for F in $FIGURES; do
        png=$(figure_png "$N" "$F" "$D")
        if [ -s "$png" ]; then echo "![$F, N=$N, day $D]($(basename "$png"))"; else echo "- $F: not dumped yet"; fi; echo
      done
    done
  } > "$md.tmp" && mv "$md.tmp" "$md"
}

pass() { # 0 once every report up to UNTIL is complete
  local k D N all complete=1 waiting=0
  for ((k = 1; ; k++)); do
    D=$(qday "$k"); [ "$D" -gt "$UNTIL" ] && break
    all=1; for N in $NS; do [ -s "$(state_of "$N" "$D")" ] || all=0; done
    [ "$all" = 1 ] || { complete=0; continue; }
    make_parts "$D" "$k"
    [ "$waiting" = 1 ] || make_figures "$D" || waiting=1
    assemble "$D" "$k"
    figures_done "$D" || complete=0
  done
  [ "$complete" = 1 ]
}

log "reports of $PREFIX N=$NS under $REPORTS (states from$(for N in $NS; do printf ' %s' "$(dir_of "$N")"; done)), $PER_YEAR a year up to day $UNTIL"
while :; do
  pass && { log "every report up to day $UNTIL is complete"; exit 0; }
  [ -n "$ONCE" ] && exit 0
  sleep "$INTERVAL"
done
