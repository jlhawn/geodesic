#!/bin/bash
# Pulls a run from the Verda instance called NAME to this machine every
# INTERVAL (600) seconds until STOP_FILE: each round asks verda for the
# instance's address (it changes when the relauncher recreates it), lists
# the run's files on the instance with scripts/verdaFiles.mjs (the
# whole-day states <PREFIX><N>_dayDDDD.bin, logs, the comparison, the
# suite and benchmark reports and markers; no .partial files, in-day
# checkpoints or benchmark states), copies them by rsync into DEST and
# compares sizes: a state here of the size listed there counts as pulled,
# a missing one is fetched next round, and one of another size is logged
# and left as it is, since a round never replaces a state already here.
# Nothing is ever deleted, here or there. A round without a running
# instance (an eviction gap), or one that fails, is logged and the next
# one tries again.
# --verify pulls once more and then compares the SHA-256 of every listed
# file on both sides, pulls whatever differs again by checksum, keeping
# the copy it replaces here as <file>.replaced-<UTC time>, and
# compares once more, writing DEST/verify.txt: each file's size and sum,
# the problems left (missing here, different size or sum) and whether the
# run had ended (ENDED_<PREFIX> on the instance); it exits 0 only when
# every file there is here and identical.
#   nohup scripts/verdaPull.sh > /dev/null 2>&1 &
#   scripts/verdaPull.sh --verify
# Environment: NAME (gcm-eleven), PREFIX (eleven), REMOTE_OUT
# (/root/runs/<PREFIX>), REMOTE_REPO (/root/geodesic, whose
# scripts/verdaFiles.mjs lists the files there), DEST (runs/verda-<PREFIX>
# of the main checkout, made if missing), INTERVAL, STREAMS (4 rsyncs of
# the states at a time, each with its share of the list), ROUNDS (stop after
# that many rounds), STOP_FILE (DEST/STOP_pull), LOG (DEST/pull.log), and
# HOST, SSH_USER, VERDA and SSH of scripts/verdaHost.sh.
cd "$(dirname "$0")/.."
REPO=$PWD
NAME=${NAME:-gcm-eleven} PREFIX=${PREFIX:-eleven}
REMOTE_OUT=${REMOTE_OUT:-/root/runs/$PREFIX} REMOTE_REPO=${REMOTE_REPO:-/root/geodesic}
MAIN=$(dirname "$(git rev-parse --path-format=absolute --git-common-dir 2>/dev/null || echo "$REPO/.git")")
DEST=${DEST:-$MAIN/runs/verda-$PREFIX}
INTERVAL=${INTERVAL:-600} STREAMS=${STREAMS:-4} STOP_FILE=${STOP_FILE:-$DEST/STOP_pull} LOG=${LOG:-$DEST/pull.log}
verify=; [ "$1" = --verify ] && verify=1
. "$REPO/scripts/verdaHost.sh"
mkdir -p "$DEST" || exit 1
work=$(mktemp -d "${TMPDIR:-/tmp}/verda-pull.XXXXXX"); trap 'rm -rf "$work"' EXIT

log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a "$LOG"; }
listing() { remote "cd '$REMOTE_REPO' && node scripts/verdaFiles.mjs '$REMOTE_OUT' '$PREFIX' $1"; }
listed_here() { node "$REPO/scripts/verdaFiles.mjs" "$DEST" "$PREFIX" $1; }

last=
round() {
  resolve || { log "$why; nothing pulled"; return 1; }
  [ "$target" != "$last" ] && log "instance at $target"
  last=$target
  listing > "$work/there" 2> "$work/error" || { log "listing $REMOTE_OUT on $target failed: $(tr '\n' ' ' < "$work/error" | cut -c1-200)"; return 1; }
  [ -s "$work/there" ] || { log "nothing in $REMOTE_OUT yet"; return 0; }
  cut -f3 "$work/there" | grep -E '\.bin$' > "$work/states"
  cut -f3 "$work/there" | grep -vE '\.bin$' > "$work/files"
  if [ -s "$work/states" ]; then
    rm -f "$work/states."*
    awk -v n="$STREAMS" -v out="$work/states." '{ print > (out NR % n) }' "$work/states"
    for part in "$work/states."*; do
      [ -s "$part" ] || continue
      rsync -a --ignore-existing --files-from="$part" -e "$rsh" "$(at "$REMOTE_OUT")/" "$DEST/" 2> "$part.error" || log "rsync: $(tr '\n' ' ' < "$part.error" | cut -c1-200)" &
    done
    wait
  fi
  if [ -s "$work/files" ]; then
    rsync -a --files-from="$work/files" -e "$rsh" "$(at "$REMOTE_OUT")/" "$DEST/" 2> "$work/error" || log "rsync: $(tr '\n' ' ' < "$work/error" | cut -c1-200)"
  fi
  listed_here > "$work/here"
  local summary
  summary=$(awk -F'\t' '
    NR == FNR { there[$3] = $1; next } { here[$3] = $1 }
    END {
      for (p in there) {
        if (p !~ /\.bin$/) { other++; continue }
        states++; bytes += there[p]
        if (!(p in here)) { missing++; names = names " " p " (not here)" }
        else if (here[p] != there[p]) { differ++; names = names " " p " (" here[p] " of " there[p] " bytes)" }
        else pulled++
      }
      printf "%d states there (%.1f GB): %d here of the same size, %d missing, %d of another size; %d other files%s", states, bytes / 1e9, pulled, missing, differ, other, (names ? ";" names : "")
    }' "$work/there" "$work/here")
  log "$summary; newest:$(awk -F'\t' -v p="$PREFIX" '$3 ~ "^" p "[0-9]+_day[0-9]+\\.bin$" { n = $3; sub(/_day.*/, "", n); if ($3 > last[n]) last[n] = $3 } END { for (n in last) printf " %s", last[n] }' "$work/there")"
  echo "$summary" | grep -qE ' 0 missing, 0 of another size'
}

if [ -n "$verify" ]; then
  round
  resolve || { log "cannot verify: $why"; exit 1; }
  listing --sha256 > "$work/there" 2> "$work/error" || { log "cannot verify: listing on $target failed: $(tr '\n' ' ' < "$work/error" | cut -c1-200)"; exit 1; }
  ended=$(awk -F'\t' -v m="ENDED_$PREFIX" '$3 == m { found = 1 } END { exit !found }' "$work/there" && remote "cat '$REMOTE_OUT/ENDED_$PREFIX'")
  check() {
    listed_here --sha256 > "$work/here"
    awk -F'\t' -v ended="$ended" -v prefix="$PREFIX" -v where="$target:$REMOTE_OUT" -v dest="$DEST" -v when="$(date -u '+%Y-%m-%dT%H:%M:%SZ')" -v repulled="$1" -v bad_list="$work/bad" '
      NR == FNR { size[$3] = $1; sum[$3] = $2; order[++n] = $3; next } { hsize[$3] = $1; hsum[$3] = $2 }
      END {
        printf "" > bad_list
        for (i = 1; i <= n; i++) {
          p = order[i]; bytes += size[p]
          if (!(p in hsize)) problem[++bad] = p ": not here"
          else if (hsize[p] != size[p]) problem[++bad] = p ": " hsize[p] " bytes here, " size[p] " there"
          else if (hsum[p] != sum[p]) problem[++bad] = p ": SHA-256 differs"
          else { ok++; continue }
          print p > bad_list
        }
        printf "verified %s: %d files on %s (%.2f GB), %d identical in %s, %d problems\n", when, n, where, bytes / 1e9, ok, dest, bad
        print (ended ? "the run had ended: " ended : "the run had NOT ended (no ENDED_" prefix " on the instance): files still being written may differ")
        if (repulled) print repulled " files that differed were pulled again by checksum before this comparison"
        if (bad) { print ""; print "problems:"; for (i = 1; i <= bad; i++) print "  " problem[i] }
        print ""; print "files (bytes, SHA-256, path), identical on both sides unless listed above:"
        for (i = 1; i <= n; i++) printf "%s\t%s\t%s\n", size[order[i]], sum[order[i]], order[i]
        exit (bad > 0)
      }' "$work/there" "$work/here" > "$DEST/verify.txt"
  }
  check
  code=$?
  if [ $code -ne 0 ] && [ -s "$work/bad" ]; then
    rsync -a --checksum --backup --suffix=".replaced-$(date -u '+%Y%m%dT%H%M%SZ')" --files-from="$work/bad" -e "$rsh" "$(at "$REMOTE_OUT")/" "$DEST/" 2> "$work/error" || log "rsync: $(tr '\n' ' ' < "$work/error" | cut -c1-200)"
    check "$(wc -l < "$work/bad" | tr -d ' ')"
    code=$?
  fi
  log "$(head -1 "$DEST/verify.txt"); report in $DEST/verify.txt"
  exit $code
fi

log "pulling $NAME:$REMOTE_OUT into $DEST every $INTERVAL s (stop at $STOP_FILE)"
r=0
while [ ! -f "$STOP_FILE" ]; do
  r=$((r + 1))
  round
  [ -n "$ROUNDS" ] && [ "$r" -ge "$ROUNDS" ] && break
  sleep "$INTERVAL"
done
[ -f "$STOP_FILE" ] && log "stopped at $STOP_FILE"
exit 0
