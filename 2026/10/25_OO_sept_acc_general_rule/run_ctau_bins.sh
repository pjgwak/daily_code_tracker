#!/usr/bin/env bash
# Prefits and final 2D fit (FINAL_S3) into roots/ and figs/.
# Usage: bash run_ctau_bins.sh [--mc] [--pr] [--np] [--mass] [--bkg] [--2d] [--redraw] [--publish]
#                              [--rapidity=both|fwd|mid] [--pt=LOW:HIGH] [--max-jobs=N]
# No stage option = all stages in order (mc, pr, np, mass, bkg, 2d).
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
stages=(); extra=()
for arg in "$@"; do
  case "$arg" in
    --mc) stages+=(mc) ;; --pr) stages+=(pr) ;; --np) stages+=(np) ;; --mass) stages+=(mass) ;; --bkg) stages+=(bkg) ;; --2d) stages+=(2d) ;;
    --redraw) extra+=(--redraw) ;; --publish) extra+=(--publish) ;;
    --rapidity=*) extra+=(--rapidity "${arg#*=}") ;; --pt=*) extra+=(--pt "${arg#*=}") ;; --max-jobs=*) extra+=(--jobs "${arg#*=}") ;;
    *) echo "Unknown option: $arg" >&2; sed -n 3,5p "$0" >&2; exit 1 ;;
  esac
done
[ ${#stages[@]} -gt 0 ] || stages=(mc pr np mass bkg 2d)
exec python3 study/production.py "${stages[@]}" "${extra[@]}"
