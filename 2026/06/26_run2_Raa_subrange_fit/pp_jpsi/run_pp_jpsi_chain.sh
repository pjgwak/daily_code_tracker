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
run_subrange_mass=true
run_subrange_ctau=true
ct_min=-8.0
ct_max=10.0
ct_step=2.0
max_slices=-1
draw_mass_slices=true
trim_sparse_fwd=true
apply_quantile_cut=true
reuse_mass_slices=false
center_veto_half_width=0.0
first_slice=0
n_slices=-1
parallel_mass_slices=false
max_parallel_jobs=4
parallel_bins=true
max_bin_jobs=12
max_log_mb=100
log_check_interval=2
mass_slices_only=false
target_pt_low=""
target_pt_high=""
target_y_low=""
target_y_high=""
subrange_slice_list=""
timestamp="$(date +"%Y%m%d_%H%M%S")"
log_dir="logs/run_pp_jpsi_chain_${timestamp}"

# Command cookbook:
#   Full production chain for all bins, bin-parallel by default with 12 jobs:
#     bash run_pp_jpsi_chain.sh
#   Override bin-parallel job count:
#     bash run_pp_jpsi_chain.sh --parallel-bins --bin-jobs=4
#   Full production chain for one bin:
#     bash run_pp_jpsi_chain.sh --bin=6.5,9.0,0.0,1.6
#   Individual pre-subrange steps for one bin:
#     bash run_pp_jpsi_chain.sh --only=mc_mass --bin=6.5,9.0,0.0,1.6
#     bash run_pp_jpsi_chain.sh --only=mass --bin=6.5,9.0,0.0,1.6
#     bash run_pp_jpsi_chain.sh --only=ctau_pr --bin=6.5,9.0,0.0,1.6
#     bash run_pp_jpsi_chain.sh --only=ctau_np --bin=6.5,9.0,0.0,1.6
#   Subrange mass + ctau only, after mc_mass/mass/ctau_pr/ctau_np already exist:
#     bash run_pp_jpsi_chain.sh --only=subrange --bin=6.5,9.0,0.0,1.6
#   Subrange mass-slice fits only; writes/refreshes sliceSummary for subrange_ctau:
#     bash run_pp_jpsi_chain.sh --only=subrange_mass --bin=6.5,9.0,0.0,1.6
#   Subrange ctau fit only; reads the saved sliceSummary:
#     bash run_pp_jpsi_chain.sh --only=subrange_ctau --bin=6.5,9.0,0.0,1.6
#   Refit selected mass slices, refresh sliceSummary, then run ctau fit:
#     bash run_pp_jpsi_chain.sh --subrange-slices=100 --bin=6.5,9.0,0.0,1.6
#     bash run_pp_jpsi_chain.sh --only=subrange_ctau --bin=6.5,9.0,0.0,1.6
#   Refit a slice range:
#     bash run_pp_jpsi_chain.sh --subrange-slices=98,100-102 --bin=6.5,9.0,0.0,1.6
#   Refit final ctau only after changing center veto/manual ctau controls:
#     bash run_pp_jpsi_chain.sh --reuse-mass-slices --bin=6.5,9.0,0.0,1.6
#   Check commands without running ROOT fits:
#     bash run_pp_jpsi_chain.sh --dry-run --only=subrange --bin=6.5,9.0,0.0,1.6
#     bash run_pp_jpsi_chain.sh --dry-run --parallel-bins --bin-jobs=4

usage() {
  cat <<USAGE
Usage: bash $(basename "$0") [options]

Runs the pp J/psi per-bin chain:
  mc_mass.C -> mass.C -> ctau_pr.C -> ctau_np.C -> subrange_mass.C -> subrange_ctau.C

Examples:
  bash $(basename "$0")
  bash $(basename "$0") --parallel-bins --bin-jobs=4
  bash $(basename "$0") --bin=6.5,9.0,0.0,1.6
  bash $(basename "$0") --only=mc_mass --bin=6.5,9.0,0.0,1.6
  bash $(basename "$0") --only=mass --bin=6.5,9.0,0.0,1.6
  bash $(basename "$0") --only=ctau_pr --bin=6.5,9.0,0.0,1.6
  bash $(basename "$0") --only=ctau_np --bin=6.5,9.0,0.0,1.6
  bash $(basename "$0") --only=subrange --bin=6.5,9.0,0.0,1.6
  bash $(basename "$0") --only=subrange_mass --bin=6.5,9.0,0.0,1.6
  bash $(basename "$0") --only=subrange_ctau --bin=6.5,9.0,0.0,1.6
  bash $(basename "$0") --subrange-slices=100 --bin=6.5,9.0,0.0,1.6
  bash $(basename "$0") --reuse-mass-slices --bin=6.5,9.0,0.0,1.6
  bash $(basename "$0") --dry-run --only=subrange --bin=6.5,9.0,0.0,1.6

Options:
  --publish           Draw publish figures for macros that support publish mode.
                      This uses saved fits and writes figures under figs_publish.
  --dry-run           Print ROOT calls without executing them.
  --parallel-bins     Run production bins in parallel. This is enabled by default.
                      Each bin still runs its selected steps serially.
  --bin-jobs=N        Max parallel production-bin jobs. Default: ${max_bin_jobs}
  --max-log-mb=N      Kill a running ROOT/bin job when its log exceeds N MB. Default: ${max_log_mb}
  --log-check-interval=N
                      Log size monitor interval in seconds. Default: ${log_check_interval}
  --only=STEP         Run one step only. STEP: mc_mass, mass, ctau_pr, ctau_np,
                      subrange_mass, subrange_ctau. STEP=subrange runs both subrange steps.
  --skip=STEP         Skip one step. Can be repeated.
  --ct-min=X          subrange ctau lower edge. Default: ${ct_min}
  --ct-max=X          subrange ctau upper edge. Default: ${ct_max}
  --ct-step=X         subrange ctau step. Default: ${ct_step}
  --max-slices=N      cap subrange slices for tests. Default: ${max_slices}
  --no-draw-slices    do not save per-slice subrange mass PDFs.
  --no-trim-sparse    disable low-pT forward ctau-range trimming in subrange.
  --no-quantile-cut   disable subrange ctau/ctauErr quantile cleanup.
  --reuse-mass-slices skip subrange mass fits and run the final ctau fit from saved sliceSummary.
  --center-veto=X     exclude |ctau| < X from the final ctau fit. Default: ${center_veto_half_width}
  --slice-index=N     start subrange processing at adaptive ctau slice index N. Default: ${first_slice}
  --n-slices=N        number of adaptive ctau slices to process. Default: ${n_slices}
  --parallel-subrange run subrange mass-slice fits in parallel ROOT subprocesses.
  --subrange-jobs=N   max parallel subrange mass-slice jobs. Default: ${max_parallel_jobs}
  --mass-slices-only  fit selected subrange mass slices and skip final ctau fit.
  --bin=PTLOW,PTHIGH,YLOW,YHIGH
                      run only one production bin.
  --pt-low=X          target one production bin. Use with --pt-high/--y-low/--y-high.
  --pt-high=X         target one production bin. Use with --pt-low/--y-low/--y-high.
  --y-low=X           target one production bin. Use with --pt-low/--pt-high/--y-high.
  --y-high=X          target one production bin. Use with --pt-low/--pt-high/--y-low.
  --subrange-slices=L fit only listed subrange mass slices in the target bin.
                      L accepts comma-separated indices/ranges, e.g. 100 or 98,100-102.
                      This implies --only=subrange_mass --mass-slices-only and skips ctau fit.
  -h, --help          show this help.

Edit bin_groups inside this script for the production bin list.
USAGE
}

disable_all_steps() {
  run_mc_mass=false
  run_mass=false
  run_ctau_pr=false
  run_ctau_np=false
  run_subrange_mass=false
  run_subrange_ctau=false
}

set_step_flag() {
  local step="$1"
  local value="$2"
  case "$step" in
    mc_mass) run_mc_mass="$value" ;;
    mass) run_mass="$value" ;;
    ctau_pr) run_ctau_pr="$value" ;;
    ctau_np) run_ctau_np="$value" ;;
    subrange_mass) run_subrange_mass="$value" ;;
    subrange_ctau) run_subrange_ctau="$value" ;;
    subrange_method|subrange)
      run_subrange_mass="$value"
      run_subrange_ctau="$value"
      ;;
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
    --parallel-bins)
      parallel_bins=true
      ;;
    --bin-jobs=*)
      max_bin_jobs="${arg#*=}"
      ;;
    --max-log-mb=*)
      max_log_mb="${arg#*=}"
      ;;
    --log-check-interval=*)
      log_check_interval="${arg#*=}"
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
    --reuse-mass-slices)
      reuse_mass_slices=true
      ;;
    --center-veto=*)
      center_veto_half_width="${arg#*=}"
      ;;
    --slice-index=*)
      first_slice="${arg#*=}"
      ;;
    --n-slices=*)
      n_slices="${arg#*=}"
      ;;
    --parallel-subrange)
      parallel_mass_slices=true
      ;;
    --subrange-jobs=*)
      max_parallel_jobs="${arg#*=}"
      ;;
    --mass-slices-only)
      mass_slices_only=true
      ;;
    --bin=*)
      bin_spec="${arg#*=}"
      bin_spec="${bin_spec//,/ }"
      read -r target_pt_low target_pt_high target_y_low target_y_high <<< "$bin_spec"
      ;;
    --pt-low=*)
      target_pt_low="${arg#*=}"
      ;;
    --pt-high=*)
      target_pt_high="${arg#*=}"
      ;;
    --y-low=*)
      target_y_low="${arg#*=}"
      ;;
    --y-high=*)
      target_y_high="${arg#*=}"
      ;;
    --subrange-slices=*)
      subrange_slice_list="${arg#*=}"
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
require_int "--slice-index" "$first_slice"
require_int "--n-slices" "$n_slices"
require_int "--subrange-jobs" "$max_parallel_jobs"
require_int "--bin-jobs" "$max_bin_jobs"
require_int "--max-log-mb" "$max_log_mb"
require_int "--log-check-interval" "$log_check_interval"
if (( max_bin_jobs < 1 )); then
  echo "[ERROR] --bin-jobs must be >= 1" >&2
  exit 1
fi
if (( max_log_mb < 1 )); then
  echo "[ERROR] --max-log-mb must be >= 1" >&2
  exit 1
fi
if (( log_check_interval < 1 )); then
  echo "[ERROR] --log-check-interval must be >= 1" >&2
  exit 1
fi

has_target_bin=false
if [[ -n "$target_pt_low" || -n "$target_pt_high" || -n "$target_y_low" || -n "$target_y_high" ]]; then
  if [[ -z "$target_pt_low" || -z "$target_pt_high" || -z "$target_y_low" || -z "$target_y_high" ]]; then
    echo "[ERROR] target bin requires all of --pt-low/--pt-high/--y-low/--y-high, or --bin=PTLOW,PTHIGH,YLOW,YHIGH" >&2
    exit 1
  fi
  has_target_bin=true
fi

if [[ -n "$subrange_slice_list" ]]; then
  if [[ "$has_target_bin" != true ]]; then
    echo "[ERROR] --subrange-slices requires a target bin: use --bin=PTLOW,PTHIGH,YLOW,YHIGH" >&2
    exit 1
  fi
  if ! [[ "$subrange_slice_list" =~ ^[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*$ ]]; then
    echo "[ERROR] --subrange-slices must be comma-separated indices/ranges, e.g. 100 or 98,100-102" >&2
    exit 1
  fi
  run_mc_mass=false
  run_mass=false
  run_ctau_pr=false
  run_ctau_np=false
  run_subrange_mass=true
  run_subrange_ctau=false
  reuse_mass_slices=false
  mass_slices_only=true
  parallel_mass_slices=false
fi

if [[ "$reuse_mass_slices" == true ]]; then
  run_subrange_mass=false
fi
if [[ "$mass_slices_only" == true ]]; then
  run_subrange_ctau=false
fi
max_log_bytes=$((max_log_mb * 1024 * 1024))
mkdir -p "$log_dir"

echo "[INFO] logs: ${log_dir}"
echo "[INFO] publish=${publish} dry_run=${dry_run} parallel_bins=${parallel_bins} bin_jobs=${max_bin_jobs}"
echo "[INFO] steps: mc_mass=${run_mc_mass} mass=${run_mass} ctau_pr=${run_ctau_pr} ctau_np=${run_ctau_np} subrange_mass=${run_subrange_mass} subrange_ctau=${run_subrange_ctau}"
echo "[INFO] subrange: ct=[${ct_min}, ${ct_max}], step=${ct_step}, max_slices=${max_slices}, draw=${draw_mass_slices}, trim_sparse_fwd=${trim_sparse_fwd}, quantile=${apply_quantile_cut}"
echo "[INFO] subrange controls: reuse_mass_slices=${reuse_mass_slices}, center_veto=${center_veto_half_width}, first_slice=${first_slice}, n_slices=${n_slices}, parallel=${parallel_mass_slices}, jobs=${max_parallel_jobs}, mass_slices_only=${mass_slices_only}"
echo "[INFO] log monitor: max_log_mb=${max_log_mb}, check_interval=${log_check_interval}s"
if [[ "$has_target_bin" == true ]]; then
  echo "[INFO] target bin: pt=${target_pt_low}-${target_pt_high}, |y|=${target_y_low}-${target_y_high}"
fi
if [[ -n "$subrange_slice_list" ]]; then
  echo "[INFO] rerun subrange mass slice list: ${subrange_slice_list}"
fi

# Format:
#   "yLow yHigh|ptLow1 ptHigh1;ptLow2 ptHigh2;..."
bin_groups=(
  "1.6 2.4|3.5 6.5;6.5 9.0;9.0 12.0;12.0 40.0;3.5 40.0"
  "0.0 1.6|6.5 9.0;9.0 12.0;12.0 15.0;15.0 20.0;20.0 25.0;25.0 40.0;6.5 40.0"
  # "1.6 2.4|3.5 6.5;"
  # "1.6 2.4|6.5 9;"
)

format_tag() {
  local value="$1"
  value="${value//./p}"
  value="${value//-/m}"
  echo "$value"
}

same_number() {
  [[ "$(printf "%.6g" "$1")" == "$(printf "%.6g" "$2")" ]]
}

bin_matches_target() {
  local pt_low="$1"
  local pt_high="$2"
  local y_low="$3"
  local y_high="$4"
  if [[ "$has_target_bin" != true ]]; then
    return 0
  fi
  same_number "$pt_low" "$target_pt_low" &&
    same_number "$pt_high" "$target_pt_high" &&
    same_number "$y_low" "$target_y_low" &&
    same_number "$y_high" "$target_y_high"
}

expand_slice_list() {
  local list="$1"
  local item start stop i
  local expanded=()
  IFS=',' read -r -a items <<< "$list"
  for item in "${items[@]}"; do
    if [[ "$item" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      start="${BASH_REMATCH[1]}"
      stop="${BASH_REMATCH[2]}"
      if (( start > stop )); then
        echo "[ERROR] bad --subrange-slices range: ${item}" >&2
        return 1
      fi
      for ((i = start; i <= stop; ++i)); do
        expanded+=("$i")
      done
    else
      expanded+=("$item")
    fi
  done
  printf '%s\n' "${expanded[@]}"
}

monitor_log_size() {
  local job_pid="$1"
  local log_file="$2"
  local label="$3"
  local max_bytes="$4"
  local interval="$5"

  while kill -0 "$job_pid" 2>/dev/null; do
    if [[ -f "$log_file" ]]; then
      local size
      size="$(stat -c '%s' "$log_file" 2>/dev/null || echo 0)"
      if [[ "$size" =~ ^[0-9]+$ ]] && (( size > max_bytes )); then
        {
          echo
          echo "[MONITOR] log exceeded ${max_bytes} bytes for ${label}; stopping job pid=${job_pid}"
        } >> "$log_file"

        pkill -TERM -P "$job_pid" 2>/dev/null || true
        kill -TERM "$job_pid" 2>/dev/null || true
        sleep 3

        if kill -0 "$job_pid" 2>/dev/null; then
          pkill -KILL -P "$job_pid" 2>/dev/null || true
          kill -KILL "$job_pid" 2>/dev/null || true
        fi

        printf '[MONITOR] killed %s at %s bytes\n' "$label" "$size" > "${log_file}.killed"
        return 0
      fi
    fi
    sleep "$interval"
  done
}

root_call() {
  local macro_file="$1"
  local call_expr="$2"
  local log_file="$3"
  local cmd=(root -l -b -q -e ".L ${macro_file}" -e "${call_expr}")

  echo
  echo "[RUN] ${call_expr}"
  # Command sidecar output disabled to avoid creating ${log_file}.cmd files.
  # printf '%q ' "${cmd[@]}" > "${log_file}.cmd"
  # printf '\n' >> "${log_file}.cmd"
  if [[ "$dry_run" -eq 1 ]]; then
    printf '[DRY-RUN] '
    printf '%q ' "${cmd[@]}"
    printf '\n'
    return 0
  fi

  "${cmd[@]}" &
  local root_pid=$!
  monitor_log_size "$root_pid" "$log_file" "$call_expr" "$max_log_bytes" "$log_check_interval" &
  local monitor_pid=$!

  wait "$root_pid"
  local root_status=$?
  wait "$monitor_pid" 2>/dev/null || true
  return "$root_status"
}

run_bin() {
  local pt_low="$1"
  local pt_high="$2"
  local y_low="$3"
  local y_high="$4"

  local y_low_tag y_high_tag pt_low_tag pt_high_tag label bin_log_file
  y_low_tag="$(format_tag "$y_low")"
  y_high_tag="$(format_tag "$y_high")"
  pt_low_tag="$(format_tag "$pt_low")"
  pt_high_tag="$(format_tag "$pt_high")"
  label="y${y_low_tag}_${y_high_tag}_pt${pt_low_tag}_${pt_high_tag}"
  bin_log_file="${log_dir}/${label}_chain.log"

  : > "$bin_log_file"

  {
  echo "============================================================"
  echo "pt: ${pt_low}-${pt_high}, |y|: ${y_low}-${y_high}"
  echo "log: ${bin_log_file}"
  echo "============================================================"

  if [[ "$run_mc_mass" == true ]]; then
    root_call "mc_mass.C" \
      "mc_mass(${pt_low},${pt_high},${y_low},${y_high},0,1,false,false,false,false,2,1,false,${publish})" \
      "$bin_log_file"
  fi

  if [[ "$run_mass" == true ]]; then
    root_call "mass.C" \
      "mass(${pt_low},${pt_high},${y_low},${y_high},false,${publish})" \
      "$bin_log_file"
  fi

  if [[ "$run_ctau_pr" == true ]]; then
    root_call "ctau_pr.C" \
      "ctau_pr(${pt_low},${pt_high},${y_low},${y_high},false,${publish})" \
      "$bin_log_file"
  fi

  if [[ "$run_ctau_np" == true ]]; then
    root_call "ctau_np.C" \
      "ctau_np(${pt_low},${pt_high},${y_low},${y_high},0.0,1,0.08,2,false,${publish})" \
      "$bin_log_file"
  fi

  if [[ "$run_subrange_mass" == true ]]; then
    if [[ -n "$subrange_slice_list" ]]; then
      local slice
      while IFS= read -r slice; do
        [[ -z "$slice" ]] && continue
        root_call "subrange_mass.C" \
          "subrange_mass(${pt_low},${pt_high},${y_low},${y_high},${ct_min},${ct_max},${ct_step},${max_slices},${draw_mass_slices},${trim_sparse_fwd},${apply_quantile_cut},false,${center_veto_half_width},${slice},1,false,${max_parallel_jobs},true)" \
          "$bin_log_file"
      done < <(expand_slice_list "$subrange_slice_list")
    else
      root_call "subrange_mass.C" \
        "subrange_mass(${pt_low},${pt_high},${y_low},${y_high},${ct_min},${ct_max},${ct_step},${max_slices},${draw_mass_slices},${trim_sparse_fwd},${apply_quantile_cut},false,${center_veto_half_width},${first_slice},${n_slices},${parallel_mass_slices},${max_parallel_jobs},${mass_slices_only})" \
        "$bin_log_file"
    fi
  fi

  if [[ "$run_subrange_ctau" == true ]]; then
    root_call "subrange_ctau.C" \
      "subrange_ctau(${pt_low},${pt_high},${y_low},${y_high},${ct_min},${ct_max},${ct_step},${max_slices},${draw_mass_slices},${trim_sparse_fwd},${apply_quantile_cut},true,${center_veto_half_width},${first_slice},${n_slices},${parallel_mass_slices},${max_parallel_jobs},false)" \
      "$bin_log_file"
  fi

  } >> "$bin_log_file" 2>&1
}

start_bin_job() {
  local pt_low="$1"
  local pt_high="$2"
  local y_low="$3"
  local y_high="$4"
  local bin_log_file="$5"
  local label="$6"

  (
    run_bin "$pt_low" "$pt_high" "$y_low" "$y_high"
  ) &
  local bin_pid=$!

  monitor_log_size "$bin_pid" "$bin_log_file" "$label" "$max_log_bytes" "$log_check_interval" &
  local monitor_pid=$!

  wait "$bin_pid"
  local bin_status=$?
  wait "$monitor_pid" 2>/dev/null || true
  return "$bin_status"
}

wait_for_running_jobs() {
  local pid
  local status

  for pid in "${running_jobs[@]}"; do
    if wait "$pid"; then
      continue
    fi
    status=$?
    bins_failed=$((bins_failed + 1))
    echo "[WARN] bin job pid=${pid} exited with status=${status}" >&2
  done
  running_jobs=()
}

bins_failed=0
bins_total=0
running_jobs=()
for group in "${bin_groups[@]}"; do
  IFS='|' read -r y_bin pt_bin_list <<< "$group"
  read -r y_low y_high <<< "$y_bin"

  IFS=';' read -r -a pt_bins <<< "$pt_bin_list"
  for pt_bin in "${pt_bins[@]}"; do
    [[ -z "${pt_bin// }" ]] && continue
    read -r pt_low pt_high <<< "$pt_bin"
    if ! bin_matches_target "$pt_low" "$pt_high" "$y_low" "$y_high"; then
      continue
    fi
    bins_total=$((bins_total + 1))
    if [[ "$parallel_bins" == true ]]; then
      y_low_tag="$(format_tag "$y_low")"
      y_high_tag="$(format_tag "$y_high")"
      pt_low_tag="$(format_tag "$pt_low")"
      pt_high_tag="$(format_tag "$pt_high")"
      bin_label="y${y_low_tag}_${y_high_tag}_pt${pt_low_tag}_${pt_high_tag}"
      bin_log_file="${log_dir}/${bin_label}_chain.log"
      echo "[RUN] bin ${bin_label} -> ${bin_log_file}"
      start_bin_job "$pt_low" "$pt_high" "$y_low" "$y_high" "$bin_log_file" "$bin_label" &
      running_jobs+=("$!")
      if (( ${#running_jobs[@]} >= max_bin_jobs )); then
        wait_for_running_jobs
      fi
    else
      if ! run_bin "$pt_low" "$pt_high" "$y_low" "$y_high"; then
        bins_failed=$((bins_failed + 1))
        echo "[ERROR] bin failed: y=${y_low}-${y_high}, pt=${pt_low}-${pt_high}" >&2
      fi
    fi
  done
done

if [[ "$parallel_bins" == true ]]; then
  wait_for_running_jobs
fi

jobs_killed=$(find "${log_dir}" -maxdepth 1 -type f -name '*.log.killed' | wc -l)

summary_file="${log_dir}/summary.log"
{
  echo "[INFO] publish=${publish} dry_run=${dry_run} parallel_bins=${parallel_bins} max_bin_jobs=${max_bin_jobs}"
  echo "[INFO] steps mc_mass=${run_mc_mass} mass=${run_mass} ctau_pr=${run_ctau_pr} ctau_np=${run_ctau_np} subrange_mass=${run_subrange_mass} subrange_ctau=${run_subrange_ctau}"
  echo "[INFO] ct_min=${ct_min} ct_max=${ct_max} ct_step=${ct_step} max_slices=${max_slices} draw_mass_slices=${draw_mass_slices} trim_sparse_fwd=${trim_sparse_fwd} apply_quantile_cut=${apply_quantile_cut}"
  echo "[INFO] reuse_mass_slices=${reuse_mass_slices} center_veto=${center_veto_half_width} first_slice=${first_slice} n_slices=${n_slices} parallel_mass_slices=${parallel_mass_slices} max_parallel_jobs=${max_parallel_jobs} mass_slices_only=${mass_slices_only}"
  echo "[INFO] target_bin=${has_target_bin} pt_low=${target_pt_low} pt_high=${target_pt_high} y_low=${target_y_low} y_high=${target_y_high}"
  echo "[INFO] subrange_slice_list=${subrange_slice_list}"
  echo "[INFO] max_log_mb=${max_log_mb} log_check_interval=${log_check_interval}"
  echo "[INFO] total_bins=${bins_total}"
  echo "[INFO] failed_bins=${bins_failed}"
  echo "[INFO] killed_jobs=${jobs_killed}"
  echo "[INFO] generated_logs=${log_dir}"
} > "$summary_file"

echo "[INFO] summary: ${summary_file}"
if (( bins_total == 0 )); then
  echo "[ERROR] no bins matched the requested selection" >&2
  exit 1
fi
if (( bins_failed > 0 )); then
  exit 1
fi
