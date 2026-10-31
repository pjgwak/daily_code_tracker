#!/usr/bin/env bash
# PNG copies of every prefit, final-fit and summary figure for the web page (docs/web/figs/); no fitting.
# Run after run_publish.sh and drawing/refresh_all.sh.
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
out=docs/web/figs; mkdir -p "$out"
png() { # png SRC.pdf DEST.png WIDTH
  pdftoppm -png -singlefile -scale-to-x "$3" -scale-to-y -1 "$1" "${2%.png}"
  convert "$2" -trim +repage "$2"
}
for y in y1.60_2.40 y0.00_1.60; do
  r=fwd; [ $y = y0.00_1.60 ] && r=mid
  for f in figs_publish/$y/fit2d/mass_fit_${y}_pt*.pdf; do
    pt=${f##*_pt}; pt=${pt%.pdf}; lo=${pt%%.*}; hi=${pt#*_}; hi=${hi%%.*}; b=${r}_${lo}_${hi}
    png "$(ls figs_publish/$y/mc_mass/mass_fit_${y}_pt${pt}*.pdf | head -1)" "$out/${b}_mc_mass.png" 842 &  # MC mass figure may carry a _bkgOff suffix
    png "figs_publish/$y/ctau_pr/lifetime_fit_${y}_pt${pt}.pdf" "$out/${b}_ctau_pr.png" 842 &
    png "figs_publish/$y/ctau_np/lifetime_fit_${y}_pt${pt}.pdf" "$out/${b}_ctau_np.png" 842 &
    png "figs_publish/$y/mass/mass_fit_${y}_pt${pt}.pdf" "$out/${b}_mass_prefit.png" 842 &
    png "figs_publish/$y/ctau_bkg/lifetime_fit_${y}_pt${pt}.pdf" "$out/${b}_ctau_bkg.png" 842 &
    png "$f" "$out/${b}_mass.png" 842 &
    png "figs_publish/$y/fit2d/lifetime_fit_${y}_pt${pt}.pdf" "$out/${b}_lifetime.png" 842 &
    png "drawing/projection_diagnostics/fit2d_mass_${y}_pt${pt}.pdf" "$out/${b}_diag_mass.png" 1200 &
    png "drawing/projection_diagnostics/fit2d_ctau_${y}_pt${pt}.pdf" "$out/${b}_diag_ctau.png" 1200 &
    wait
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
for s in bfrac sys_np sys_pr; do cp docs/figs/final_${s}.png "$out/sum_v3_${s}.png"; cp docs/figs/baseline_${s}.png "$out/sum_v2_${s}.png"; done
echo "updated $out ($(ls $out | wc -l) files)"
