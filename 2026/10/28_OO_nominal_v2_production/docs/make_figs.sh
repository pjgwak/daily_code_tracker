#!/usr/bin/env bash
# PNG copies of the final-fit and summary figures for docs/ (REPORT.md, oo_fit_rules.html); no fitting.
# Run after run_publish.sh and drawing/refresh_all.sh. Prefit figures (*_mc_mass, *_mass_prefit, *_ctau_pr, *_ctau_np,
# *_ctau_bkg) are not touched: the prefits are those of the 25 folder (bit-identical rerun, docs/FINAL_V2.md).
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
out=docs/figs
png() { # png SRC.pdf DEST.png WIDTH
  pdftoppm -png -singlefile -scale-to-x "$3" -scale-to-y -1 "$1" "${2%.png}"
  convert "$2" -trim +repage "$2"
}
for y in y1.60_2.40 y0.00_1.60; do
  r=fwd; [ $y = y0.00_1.60 ] && r=mid
  for f in figs_publish/$y/fit2d/mass_fit_${y}_pt*.pdf; do
    pt=${f##*_pt}; pt=${pt%.pdf}; lo=${pt%%.*}; hi=${pt#*_}; hi=${hi%%.*}; b=${r}_${lo}_${hi}
    png "$f" "$out/${b}_mass.png" 842
    png "figs_publish/$y/fit2d/lifetime_fit_${y}_pt${pt}.pdf" "$out/${b}_lifetime.png" 842
    png "drawing/projection_diagnostics/fit2d_mass_${y}_pt${pt}.pdf" "$out/${b}_diag_mass.png" 1200
    png "drawing/projection_diagnostics/fit2d_ctau_${y}_pt${pt}.pdf" "$out/${b}_diag_ctau.png" 1200
  done
done
for b in 1_2 2_3; do
  lo=${b%_*}; hi=${b#*_}
  png "drawing/ctau_zoom/lifetime_fit_y1.60_2.40_pt${lo}.00_${hi}.00.pdf" "$out/fwd_${b}_lifetime_zoom.png" 842
done
for r in fwd mid; do
  png "drawing/b_frac/bFrac_OO_${r}.pdf" "$out/sum_bfrac_${r}.png" 900
  png "drawing/raw_yields/raw_yield_OO_${r}_ratio.pdf" "$out/sum_yield_ratio_${r}.png" 900
  png "drawing/b_frac_comparison/bfraction_nominal_with_${r}_syst.pdf" "$out/sum_nominal_syst_${r}.png" 900
done
png drawing/b_frac_comparison/bfraction_production_vs_nominal.pdf "$out/sum_vs_nominal_fwd.png" 900
png drawing/b_frac_comparison/bfraction_mid_production_vs_nominal.pdf "$out/sum_vs_nominal_mid.png" 900
for s in mc_mass ctau_pr ctau_np mass ctau_bkg; do png "drawing/model_count/model_count_${s}.pdf" "$out/sum_model_count_${s}.png" 900; done
png drawing/projection_diagnostics/binning_overview.pdf "$out/sum_binning_overview.png" 1200
echo "updated $out"
