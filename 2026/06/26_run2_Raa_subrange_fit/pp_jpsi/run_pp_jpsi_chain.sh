#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}"

dry_run=0
publish=false
run_mc_mass=true
run_mass=true
run_ctau_pr=true
run_ctau_np=true
run_subrange=true
ct_min=-8.0
ct_max=10.0
ct_step=2.0
max_slices=-1
draw_mass_slices=true
trim_sparse_fwd=true
apply_quantile_cut=true
timestamp="$(date +"%Y%m%d_%H%M%S")"
log_dir="logs/run_pp_jpsi_chain_${timestamp}"

usage() {
  cat <<USAGE
Usage: bash $(basename "$0") [options]

Runs the pp J/psi per-bin chain:
  mc_mass.C -> mass.C -> ctau_pr.C -> ctau_np.C -> subrange_method.C

Options:
  --publish           Draw publish figures for macros that support publish mode.
                      This uses saved fits and writes figures under figs_publish.
  --dry-run           Print ROOT calls without executing them.
  --only=STEP         Run one step only. STEP: mc_mass, mass, ctau_pr, ctau_np, subrange_method.
  --skip=STEP         Skip one step. Can be repeated.
  --ct-min=X          subrange ctau lower edge. Default: ${ct_min}
  --ct-max=X          subrange ctau upper edge. Default: ${ct_max}
  --ct-step=X         subrange ctau step. Default: ${ct_step}
  --max-slices=N      cap subrange slices for tests. Default: ${max_slices}
  --no-draw-slices    do not save per-slice subrange mass PDFs.
  --no-trim-sparse    disable low-pT forward ctau-range trimming in subrange.
  --no-quantile-cut   disable subrange ctau/ctauErr quantile cleanup.
  -h, --help          show this help.

Edit bin_groups inside this script for the production bin list.
USAGE
}

disable_all_steps() {
  run_mc_mass=false
  run_mass=false
  run_ctau_pr=false
  run_ctau_np=false
  run_subrange=false
}

set_step_flag() {
  local step="$1"
  local value="$2"
  case "$step" in
    mc_mass) run_mc_mass="$value" ;;
    mass) run_mass="$value" ;;
    ctau_pr) run_ctau_pr="$value" ;;
    ctau_np) run_ctau_np="$value" ;;
    subrange_method|subrange) run_subrange="$value" ;;
    *)
      echo "[ERROR] unknown step: ${step}" >&2
      usage >&2
      exit 1
      ;;
  esac
}

for arg in "$@"; do
  case "$arg" in
    --publish)
      publish=true
      ;;
    --dry-run)
      dry_run=1
      ;;
    --only=*)
      disable_all_steps
      set_step_flag "${arg#*=}" true
      ;;
    --skip=*)
      set_step_flag "${arg#*=}" false
      ;;
    --ct-min=*)
      ct_min="${arg#*=}"
      ;;
    --ct-max=*)
      ct_max="${arg#*=}"
      ;;
    --ct-step=*)
      ct_step="${arg#*=}"
      ;;
    --max-slices=*)
      max_slices="${arg#*=}"
      ;;
    --no-draw-slices)
      draw_mass_slices=false
      ;;
    --no-trim-sparse)
      trim_sparse_fwd=false
      ;;
    --no-quantile-cut)
      apply_quantile_cut=false
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[ERROR] unknown option: $arg" >&2
      usage >&2
      exit 1
      ;;
  esac
done

require_int() {
  local name="$1"
  local value="$2"
  if ! [[ "$value" =~ ^-?[0-9]+$ ]]; then
    echo "[ERROR] ${name} must be an integer: ${value}" >&2
    exit 1
  fi
}

require_int "--max-slices" "$max_slices"
mkdir -p "$log_dir"

echo "[INFO] logs: ${log_dir}"
echo "[INFO] publish=${publish} dry_run=${dry_run}"
echo "[INFO] steps: mc_mass=${run_mc_mass} mass=${run_mass} ctau_pr=${run_ctau_pr} ctau_np=${run_ctau_np} subrange_method=${run_subrange}"
echo "[INFO] subrange: ct=[${ct_min}, ${ct_max}], step=${ct_step}, max_slices=${max_slices}, draw=${draw_mass_slices}, trim_sparse_fwd=${trim_sparse_fwd}, quantile=${apply_quantile_cut}"

# Format:
#   "yLow yHigh|ptLow1 ptHigh1;ptLow2 ptHigh2;..."
bin_groups=(
  # "0.0 1.6|6.5 9.0;"
  # "0.0 1.6|6.5 9.0;9.0 12.0"
  # "1.6 2.4|3.0 6.5;6.5 12.0"
  "1.6 2.4|3.0 6.5;"
)

format_tag() {
  local value="$1"
  value="${value//./p}"
  value="${value//-/m}"
  echo "$value"
}

root_call() {
  local macro_file="$1"
  local call_expr="$2"
  local log_file="$3"
  local cmd=(root -l -b -q -e ".L ${macro_file}" -e "${call_expr}")

  echo "[RUN] ${call_expr} -> ${log_file}"
  printf '%q ' "${cmd[@]}" > "${log_file}.cmd"
  printf '\n' >> "${log_file}.cmd"
  if [[ "$dry_run" -eq 1 ]]; then
    printf '[DRY-RUN] ' | tee "$log_file"
    printf '%q ' "${cmd[@]}" | tee -a "$log_file"
    printf '\n' | tee -a "$log_file"
    return 0
  fi
  "${cmd[@]}" > "$log_file" 2>&1
}

run_bin() {
  local pt_low="$1"
  local pt_high="$2"
  local y_low="$3"
  local y_high="$4"

  local y_low_tag y_high_tag pt_low_tag pt_high_tag label
  y_low_tag="$(format_tag "$y_low")"
  y_high_tag="$(format_tag "$y_high")"
  pt_low_tag="$(format_tag "$pt_low")"
  pt_high_tag="$(format_tag "$pt_high")"
  label="y${y_low_tag}_${y_high_tag}_pt${pt_low_tag}_${pt_high_tag}"

  echo "============================================================"
  echo "pt: ${pt_low}-${pt_high}, |y|: ${y_low}-${y_high}"
  echo "============================================================"

  if [[ "$run_mc_mass" == true ]]; then
    root_call "mc_mass.C" \
      "mc_mass(${pt_low},${pt_high},${y_low},${y_high},0,1,false,false,false,false,2,1,false,${publish})" \
      "${log_dir}/${label}_mc_mass.log"
  fi

  if [[ "$run_mass" == true ]]; then
    root_call "mass.C" \
      "mass(${pt_low},${pt_high},${y_low},${y_high},false,${publish})" \
      "${log_dir}/${label}_mass.log"
  fi

  if [[ "$run_ctau_pr" == true ]]; then
    root_call "ctau_pr.C" \
      "ctau_pr(${pt_low},${pt_high},${y_low},${y_high},false,${publish})" \
      "${log_dir}/${label}_ctau_pr.log"
  fi

  if [[ "$run_ctau_np" == true ]]; then
    root_call "ctau_np.C" \
      "ctau_np(${pt_low},${pt_high},${y_low},${y_high},0.0,1,0.08,2,false,${publish})" \
      "${log_dir}/${label}_ctau_np.log"
  fi

  if [[ "$run_subrange" == true ]]; then
    root_call "subrange_method.C" \
      "subrange_method(${pt_low},${pt_high},${y_low},${y_high},${ct_min},${ct_max},${ct_step},${max_slices},${draw_mass_slices},${trim_sparse_fwd},${apply_quantile_cut})" \
      "${log_dir}/${label}_subrange_method.log"
  fi
}

bins_failed=0
bins_total=0
for group in "${bin_groups[@]}"; do
  IFS='|' read -r y_bin pt_bin_list <<< "$group"
  read -r y_low y_high <<< "$y_bin"

  IFS=';' read -r -a pt_bins <<< "$pt_bin_list"
  for pt_bin in "${pt_bins[@]}"; do
    [[ -z "${pt_bin// }" ]] && continue
    read -r pt_low pt_high <<< "$pt_bin"
    bins_total=$((bins_total + 1))
    if ! run_bin "$pt_low" "$pt_high" "$y_low" "$y_high"; then
      bins_failed=$((bins_failed + 1))
      echo "[ERROR] bin failed: y=${y_low}-${y_high}, pt=${pt_low}-${pt_high}" >&2
    fi
  done
done

summary_file="${log_dir}/summary.log"
{
  echo "[INFO] publish=${publish} dry_run=${dry_run}"
  echo "[INFO] steps mc_mass=${run_mc_mass} mass=${run_mass} ctau_pr=${run_ctau_pr} ctau_np=${run_ctau_np} subrange_method=${run_subrange}"
  echo "[INFO] ct_min=${ct_min} ct_max=${ct_max} ct_step=${ct_step} max_slices=${max_slices} draw_mass_slices=${draw_mass_slices} trim_sparse_fwd=${trim_sparse_fwd} apply_quantile_cut=${apply_quantile_cut}"
  echo "[INFO] total_bins=${bins_total}"
  echo "[INFO] failed_bins=${bins_failed}"
  echo "[INFO] generated_logs=${log_dir}"
} > "$summary_file"

echo "[INFO] summary: ${summary_file}"
if (( bins_failed > 0 )); then
  exit 1
fi
