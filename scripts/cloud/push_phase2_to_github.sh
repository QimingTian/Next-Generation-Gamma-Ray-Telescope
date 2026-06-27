#!/usr/bin/env bash
# Run on MaxCloudON: pack Phase 2 CSVs, split tarball, push to GitHub.
set -euo pipefail
ENERGY="${1:-30}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DATA="$ROOT/data"
BD="$ROOT/build/data"
PART_PREFIX="$DATA/phase2_E${ENERGY}_events.tar.gz.part-"

cd "$ROOT"
bash /tmp/pack_phase2.sh "$ENERGY" 2>/dev/null || bash scripts/cloud/_pack_phase2.sh "$ENERGY"

TAR="$DATA/phase2_events.tar.gz"
DEST="$DATA/phase2_E${ENERGY}_events.tar.gz"
mv -f "$TAR" "$DEST"
rm -f "${PART_PREFIX}"*
split -b 95m "$DEST" "$PART_PREFIX"
ls -lh "$DEST" "${PART_PREFIX}"*

git fetch -q origin
git reset -q --hard origin/main
git pull -q --rebase origin "$(git rev-parse --abbrev-ref HEAD)" 2>/dev/null || git pull -q
git add "$DEST" "${PART_PREFIX}"*
git add scripts/cloud/push_phase2_to_github.sh scripts/cloud/fetch_phase2_from_github.sh 2>/dev/null || true
if git diff --cached --quiet; then
  echo "Nothing new to push"
  exit 0
fi
git commit -m "Add Phase 2 E${ENERGY} GeV direction_scan tarball (split for GitHub)."
git push origin HEAD
echo "Pushed Phase 2 E${ENERGY} GeV to GitHub"
