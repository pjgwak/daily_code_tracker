#!/usr/bin/env bash
#
# Usage:
#   bash run_ctau_binned_mass.sh [options]
#
# Examples:
#   bash run_ctau_binned_mass.sh
#   bash run_ctau_binned_mass.sh --jobs=8 --logy
#   bash run_ctau_binned_mass.sh --dry-run
#
# Options:
#   --jobs=N              Number of ROOT jobs to run in parallel. Default: 4
#   --logy                Save mass-slice plots with log y-axis.
#   --dry-run             Print commands without running ROOT.
#   -h, --help            Show help.
#
# Outputs:
#   logs/run_ctau_binned_mass_YYYYMMDD_HHMMSS/*.log
#   logs/run_ctau_binned_mass_YYYYMMDD_HHMMSS/*.cmd

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}"

jobs=4
dry_run=0
use_log_y=false
timestamp="$(date +"%Y%m%d_%H%M%S")"
log_dir="logs/run_ctau_binned_mass_${timestamp}"

usage() {
  cat <<USAGE
Usage: bash $(basename "$0") [options]

Runs draw_ctau_binned_mass.C for the Run 2 pp J/psi RAA pT bins.
Centrality rows are intentionally ignored here because the pp input dataset
has no centrality variable; bins are grouped by rapidity.

Options:
  --jobs=N              Number of ROOT jobs to run in parallel. Default: ${jobs}
  --logy                Save mass-slice plots with log y-axis.
  --dry-run             Print commands without running ROOT.
  -h, --help            Show this help.
USAGE
}

for arg in "$@"; do
  case "$arg" in
    --jobs=*)
      jobs="${arg#*=}"
      ;;
    --logy)
      use_log_y=true
      ;;
    --dry-run)
      dry_run=1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[ERROR] unknown option: ${arg}" >&2
      usage >&2
      exit 1
      ;;
  esac
done

require_positive_int() {
  local name="$1"
  local value="$2"
  if ! [[ "${value}" =~ ^[0-9]+$ ]] || (( value < 1 )); then
    echo "[ERROR] ${name} must be a positive integer: ${value}" >&2
    exit 1
  fi
}

require_positive_int "--jobs" "${jobs}"

# Format:
#   "yLow yHigh|ptLow1 ptHigh1;ptLow2 ptHigh2;..."
# Centrality is ignored for this pp-only diagnostic. These are the unique pT
# bins after collapsing the centrality rows in agents.md.
bin_groups=(
  "1.6 2.4|3.5 6.5;6.5 9.0;9.0 12.0;12.0 40.0;3.5 40.0"
  "0.0 1.6|6.5 9.0;9.0 12.0;12.0 15.0;15.0 20.0;20.0 25.0;25.0 40.0;6.5 40.0"
)

format_tag() {
  local value="$1"
  value="${value//./p}"
  value="${value//-/m}"
  echo "${value}"
}

mkdir -p "${log_dir}"

echo "[INFO] rapidity groups: ${#bin_groups[@]}"
echo "[INFO] parallel jobs: ${jobs}"
echo "[INFO] logs: ${log_dir}"

run_bin() {
  local y_low="$1"
  local y_high="$2"
  local pt_low="$3"
  local pt_high="$4"
  local label log_file
  label="y$(format_tag "${y_low}")_$(format_tag "${y_high}")_pt$(format_tag "${pt_low}")_$(format_tag "${pt_high}")"
  log_file="${log_dir}/${label}.log"

  local cmd=(root -l -b -q "draw_ctau_binned_mass.C(${pt_low},${pt_high},${y_low},${y_high},${use_log_y})")
  printf '%q ' "${cmd[@]}" > "${log_file}.cmd"
  printf '\n' >> "${log_file}.cmd"

  echo "[RUN] ${label} -> ${log_file}"
  if [[ "${dry_run}" -eq 1 ]]; then
    {
      printf '[DRY-RUN] '
      printf '%q ' "${cmd[@]}"
      printf '\n'
    } > "${log_file}"
    return 0
  fi

  "${cmd[@]}" > "${log_file}" 2>&1
}

active_pids=()
active_labels=()
active_logs=()
failed=0
bins_total=0

wait_for_oldest() {
  local pid="${active_pids[0]}"
  local label="${active_labels[0]}"
  local log_file="${active_logs[0]}"

  if wait "${pid}"; then
    echo "[DONE] ${label}"
  else
    echo "[ERROR] ${label} failed: ${log_file}" >&2
    failed=$((failed + 1))
  fi

  active_pids=("${active_pids[@]:1}")
  active_labels=("${active_labels[@]:1}")
  active_logs=("${active_logs[@]:1}")
}

for group in "${bin_groups[@]}"; do
  IFS='|' read -r y_bin pt_bin_list <<< "${group}"
  read -r y_low y_high <<< "${y_bin}"
  IFS=';' read -r -a pt_bins <<< "${pt_bin_list}"

  for pt_bin in "${pt_bins[@]}"; do
    [[ -z "${pt_bin// }" ]] && continue
    read -r pt_low pt_high <<< "${pt_bin}"
    label="y$(format_tag "${y_low}")_$(format_tag "${y_high}")_pt$(format_tag "${pt_low}")_$(format_tag "${pt_high}")"
    log_file="${log_dir}/${label}.log"

    bins_total=$((bins_total + 1))
    run_bin "${y_low}" "${y_high}" "${pt_low}" "${pt_high}" &
    active_pids+=("$!")
    active_labels+=("${label}")
    active_logs+=("${log_file}")

    if (( ${#active_pids[@]} >= jobs )); then
      wait_for_oldest
    fi
  done
done

while (( ${#active_pids[@]} > 0 )); do
  wait_for_oldest
done

summary_file="${log_dir}/summary.log"
{
  echo "[INFO] rapidity_groups=${#bin_groups[@]}"
  echo "[INFO] total_bins=${bins_total}"
  echo "[INFO] jobs=${jobs}"
  echo "[INFO] dry_run=${dry_run}"
  echo "[INFO] use_log_y=${use_log_y}"
  echo "[INFO] failed_bins=${failed}"
  echo "[INFO] generated_logs=${log_dir}"
} > "${summary_file}"

echo "[INFO] summary: ${summary_file}"
if (( failed > 0 )); then
  exit 1
fi
