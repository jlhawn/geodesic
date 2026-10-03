#!/bin/bash
# Runs the test files as concurrent runners, one node --test process per
# file and JOBS (the core count) at a time, each under TIMEOUT_MIN (60)
# minutes where `timeout` exists, its TAP output kept in DIR/<file>.tap,
# and writes REPORT: the files passed and failed, the tests passed and
# failed, every failing test's name with its TAP block (error, expected
# and actual values, stack), the last lines of a file that failed without
# naming a test (a crash or a timeout), the file's last 40 lines of console
# output, and every file's seconds, slowest first. Exits 1 when anything
# failed.
#   scripts/suiteReport.sh runs/suite.txt
#   FILES="test/grid.test.mjs test/mesh.test.mjs" JOBS=4 scripts/suiteReport.sh /tmp/r.txt
# Environment: FILES (test/*.test.mjs), JOBS, TIMEOUT_MIN, DIR (REPORT less
# its .txt).
cd "$(dirname "$0")/.."
REPORT=${1:-runs/suite.txt}
DIR=${DIR:-${REPORT%.txt}}
JOBS=${JOBS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)}
FILES=${FILES:-$(ls test/*.test.mjs)}
TIMEOUT_MIN=${TIMEOUT_MIN:-60}
LIMIT=; command -v timeout > /dev/null && LIMIT="timeout ${TIMEOUT_MIN}m"
mkdir -p "$(dirname "$REPORT")" && rm -rf "$DIR" && mkdir -p "$DIR"
unset NODE_TEST_CONTEXT

start=$(date +%s)
printf '%s\n' $FILES | xargs -P "$JOBS" -I{} bash -c '
  f=$1; name=$(basename "$f" .test.mjs); s=$(date +%s)
  $3 node --test --test-reporter=tap "$f" > "$2/$name.tap" 2>&1
  echo "$? $(( $(date +%s) - s )) $f" > "$2/$name.status"' _ {} "$DIR" "$LIMIT"
wall=$(( $(date +%s) - start ))

failed= files=0 bad=0 pass=0 fail=0
for status in "$DIR"/*.status; do
  read -r code seconds file < "$status"
  files=$((files + 1))
  tap=${status%.status}.tap
  p=$(grep -E '^# pass [0-9]+' "$tap" | tail -1 | awk '{print $3}'); f=$(grep -E '^# fail [0-9]+' "$tap" | tail -1 | awk '{print $3}')
  pass=$((pass + ${p:-0})) fail=$((fail + ${f:-0}))
  [ "$code" -ne 0 ] && { bad=$((bad + 1)); failed="$failed $status"; }
done

{
  echo "test suite at $(git rev-parse --short HEAD 2>/dev/null || echo unknown) on $(hostname), $(date -u '+%Y-%m-%dT%H:%M:%SZ'): $files files, $((files - bad)) passed, $bad failed; tests: $pass passed, $fail failed; ${wall} s with $JOBS runners"
  echo "node $(node --version); full output of each file in $DIR/<file>.tap"
  if [ -n "$failed" ]; then
    echo; echo "failing tests:"
    for status in $failed; do
      read -r code seconds file < "$status"
      tap=${status%.status}.tap
      echo; echo "== $file (exit $code$([ "$code" -eq 124 ] && echo ', timed out'), $seconds s)"
      blocks=$(awk '
        block { print; if ($0 == endmark) block = 0; next }
        /^ *not ok / { block = 1; indent = match($0, /[^ ]/) - 1; endmark = ""; for (i = 0; i < indent + 2; i++) endmark = endmark " "; endmark = endmark "..."; print }
      ' "$tap")
      if [ -n "$blocks" ]; then echo "$blocks"; else echo "(no failing test named; the file's last lines)"; tail -n 30 "$tap"; fi
      output=$(grep -E '^ *# ' "$tap" | grep -vE '^ *# Subtest: |^# (tests|suites|pass|fail|cancelled|skipped|todo|duration_ms) ' | tail -n 40)
      [ -n "$output" ] && { echo "-- the file's console output (last 40 lines):"; echo "$output"; }
    done
  fi
  echo; echo "seconds per file, slowest first:"
  cat "$DIR"/*.status | sort -k2,2nr | awk '{ printf "%6d  %s%s\n", $2, $3, ($1 != 0 ? "  (failed)" : "") }'
} > "$REPORT"
head -1 "$REPORT"
[ "$bad" -eq 0 ]
