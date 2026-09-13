#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run.sh — run the validation study, tier by tier.
#
#   ./run.sh smoke                       # ~10 min, sizes the main tier
#   ./run.sh main                        # ~2-3 h
#   ./run.sh ext                         # optional n = 200 arm
#   EMPHASIS_LIB=/path/to/rlib ./run.sh smoke
#   ./run.sh main --dry-run              # job table only, nothing executed
#
# Order is mandatory: tier 0 (references) gates everything; 01 and 02 are
# deterministic and build-independent; 03 is resumable, so re-running it after
# an interruption continues where it stopped.
#
# EMPHASIS_LIB selects the library emphasis is loaded from.  Leave it unset to
# use the normal library.  Numbers from a pre-wave-1 build are not final.
# ---------------------------------------------------------------------------
set -euo pipefail

TIER="${1:-smoke}"; shift || true
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R="$HERE/R"
WORKERS="${WORKERS:-10}"
LIBARG=()
if [[ -n "${EMPHASIS_LIB:-}" ]]; then LIBARG=(--lib "$EMPHASIS_LIB"); fi

case "$TIER" in
  smoke|main|ext) ;;
  *) echo "usage: $0 {smoke|main|ext} [extra args passed to 03-fit.R]" >&2; exit 2 ;;
esac

echo "=== tier $TIER  workers=$WORKERS  lib=${EMPHASIS_LIB:-<default>} ==="

# --- tier 0: reference identities and build sentinels (hard gate) ----------
if [[ ! -f "$HERE/results/GATE.ok" ]]; then
  echo "--- 00-selfcheck (tier 0) ---"
  Rscript "$R/00-selfcheck.R" "${LIBARG[@]}"
else
  echo "--- 00-selfcheck: results/GATE.ok present, skipping ---"
fi

echo "--- 01-simulate ---"
Rscript "$R/01-simulate.R" --tier "$TIER" "${LIBARG[@]}"

echo "--- 02-reference ---"
Rscript "$R/02-reference.R" --tier "$TIER" --workers "$WORKERS" "${LIBARG[@]}"

echo "--- 03-fit ---"
Rscript "$R/03-fit.R" --tier "$TIER" --workers "$WORKERS" "${LIBARG[@]}" "$@"

echo "--- 04-analyse ---"
Rscript "$R/04-analyse.R" --tier "$TIER" "${LIBARG[@]}"

echo "=== done: $HERE/results/$TIER/report.md ==="
