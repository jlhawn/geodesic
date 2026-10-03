#!/bin/bash
# Puts this checkout's commit on the Verda instance called NAME for
# scripts/verdaEleven.sh. A partial, sparse clone of HEAD is staged in a
# temporary directory: every tracked file but runs/ (the page's default
# states), with a .git that holds the commit and the blobs checked out, so
# that `git rev-parse HEAD` there names the commit the run used. It goes
# to REMOTE_REPO (/root/geodesic) by rsync, which deletes nothing there;
# rsync itself is installed on the instance first if it is missing, and
# node_modules is left to verdaEleven.sh setup. The data the model reads
# is tracked and travels with it (data/topography_0p25.bin,
# woa_annual_1deg.bin, subgrid_N*.bin, radiationBenchmark.json); after
# the copy the instance's HEAD and the size of every file in data/ are
# checked against the stage. A checkout with uncommitted changes to
# tracked files is refused unless ALLOW_DIRTY=1 (the stage is the
# commit either way).
#   scripts/verdaPush.sh
# Environment: NAME (gcm-eleven), REMOTE_REPO (/root/geodesic), ALLOW_DIRTY,
# STAGE (the directory to stage in, made and removed when unset), and
# HOST, SSH_USER, VERDA and SSH of scripts/verdaHost.sh.
cd "$(dirname "$0")/.."
REPO=$PWD
NAME=${NAME:-gcm-eleven} REMOTE_REPO=${REMOTE_REPO:-/root/geodesic}
. "$REPO/scripts/verdaHost.sh"

if [ -n "$(git status --porcelain --untracked-files=no)" ] && [ "$ALLOW_DIRTY" != 1 ]; then
  echo "uncommitted changes to tracked files; commit them or set ALLOW_DIRTY=1 to push HEAD without them" >&2; exit 1
fi
commit=$(git rev-parse HEAD)
if [ -z "$STAGE" ]; then STAGE=$(mktemp -d "${TMPDIR:-/tmp}/verda-push.XXXXXX"); trap 'rm -rf "$STAGE"' EXIT; fi
rm -rf "$STAGE/repo"
git clone -q --no-checkout --depth 1 --filter=blob:none --upload-pack='git -c uploadpack.allowFilter=true upload-pack' "file://$REPO" "$STAGE/repo" 2>/dev/null \
  && git -C "$STAGE/repo" sparse-checkout set --no-cone '/*' '!/runs/' 2>/dev/null \
  && git -C "$STAGE/repo" checkout -q "$commit" 2>/dev/null || { echo "staging $commit failed" >&2; exit 1; }
[ "$(git -C "$STAGE/repo" rev-parse HEAD)" = "$commit" ] || { echo "the stage is not at $commit" >&2; exit 1; }
for f in data/topography_0p25.bin data/woa_annual_1deg.bin data/subgrid_N64.bin data/subgrid_N128.bin; do
  [ -s "$STAGE/repo/$f" ] || { echo "$f is not in $commit" >&2; exit 1; }
done
kb() { du -sk "$1" | awk '{ print $1 }'; }
files=$(cd "$STAGE/repo" && find . -path ./.git -prune -o -type f -print | wc -l | tr -d ' ')
all=$(kb "$STAGE/repo") git_kb=$(kb "$STAGE/repo/.git") data_kb=$(kb "$STAGE/repo/data")
echo "staged $(git log -1 --format='%h %s' "$commit" | cut -c1-100): $files files, $((all / 1024)) MB ($((data_kb / 1024)) MB of data/, $((git_kb / 1024)) MB of .git; runs/ left out)"

resolve || { echo "$why" >&2; exit 1; }
if [ "$target" != local ]; then
  remote "command -v rsync > /dev/null || { apt-get -o DPkg::Lock::Timeout=600 update -qq && DEBIAN_FRONTEND=noninteractive apt-get -o DPkg::Lock::Timeout=600 install -y -qq rsync; } > /dev/null" || { echo "rsync is missing on $target and could not be installed" >&2; exit 1; }
fi
remote "mkdir -p '$REMOTE_REPO'" || { echo "cannot reach $target" >&2; exit 1; }
rsync -a --exclude node_modules -e "$rsh" "$STAGE/repo/" "$(at "$REMOTE_REPO")/" || { echo "rsync to $target failed" >&2; exit 1; }

sizes() { (cd "$1" && git rev-parse HEAD && for f in data/*; do echo "$(wc -c < "$f" | tr -d ' ') $f"; done); }
here=$(sizes "$STAGE/repo")
there=$(remote "cd '$REMOTE_REPO' && git rev-parse HEAD && for f in data/*; do echo \"\$(wc -c < \"\$f\" | tr -d ' ') \$f\"; done")
if [ "$here" != "$there" ]; then
  echo "the instance's copy differs from the stage:" >&2; diff <(echo "$here") <(echo "$there") >&2; exit 1
fi
echo "pushed to $(at "$REMOTE_REPO"): HEAD $(git rev-parse --short "$commit") and $(($(echo "$here" | wc -l) - 1)) data files of the same size there"
