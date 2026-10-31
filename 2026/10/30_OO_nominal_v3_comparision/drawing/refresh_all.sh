#!/usr/bin/env bash
# Re-draw the saved production fits (FINAL_V4) and all summaries; no minimization.
set -euo pipefail
DRAWING_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
BASE_DIR="$(cd "${DRAWING_DIR}/.." && pwd -P)"
cd "${BASE_DIR}"
LOG_DIR="${DRAWING_DIR}/logs/refresh_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${LOG_DIR}"
echo "Drawing saved ROOT results; logs: ${LOG_DIR}"
bash run_publish.sh --max-jobs="${MAX_JOBS:-4}" > "${LOG_DIR}/saved_fits.log" 2>&1
# fit2d replay with OO_DIAGNOSTICS_DIR (also refreshes figs_publish/*/fit2d); prefit macros have no diagnostics hook
python3 drawing/run_projection_diagnostics.py --jobs="${MAX_JOBS:-4}" > "${LOG_DIR}/projection_replay.log" 2>&1
python3 drawing/summarize_projection_diagnostics.py > "${LOG_DIR}/projection_diagnostics.log" 2>&1
for region in fwd mid; do
  root -n -l -b -q "drawing/b_frac/bFraction_OO_${region}_pp.C()" > "${LOG_DIR}/bfrac_${region}.log" 2>&1
  root -n -l -b -q "drawing/raw_yields/raw_yield_${region}.C()" > "${LOG_DIR}/yield_${region}.log" 2>&1
  root -n -l -b -q "drawing/b_frac_comparison/plot_bfraction_vs_nominal.C(\"\",\"\",2,\"${region}\")" > "${LOG_DIR}/comparison_${region}.log" 2>&1
  cp "drawing/b_frac/bFrac_OO_${region}.pdf" drawing/all_plots_in_one/
  cp "drawing/raw_yields/raw_yield_OO_${region}.pdf" drawing/all_plots_in_one/
done
root -n -l -b -q 'drawing/b_frac_comparison/plot_bfraction_nominal_with_fwd_syst.C()' > "${LOG_DIR}/fwd_syst.log" 2>&1
root -n -l -b -q 'drawing/b_frac_comparison/plot_bfraction_nominal_with_mid_syst.C()' > "${LOG_DIR}/mid_syst.log" 2>&1
root -n -l -b -q 'drawing/model_count/plot_model_counts.C()' > "${LOG_DIR}/model_counts.log" 2>&1
(cd drawing/fit_convergence && root -n -l -b -q 'check_roofit_results.C("../../roots",1e-3,false)') > "${LOG_DIR}/convergence.log" 2>&1
(cd drawing/ctau_zoom && OO_RULE_DIR="${BASE_DIR}/study/rules/FINAL_V4" root -n -l -b -q 'plot_ctau_zoom.C+()' && rm -f plot_ctau_zoom_C*) > "${LOG_DIR}/ctau_zoom.log" 2>&1
bash drawing/all_plots_in_one/run_all_panels.sh > "${LOG_DIR}/panels.log" 2>&1
echo "Updated B fractions, yields, nominal comparisons, model counts, convergence table, zooms, diagnostics, and all fit panels."
