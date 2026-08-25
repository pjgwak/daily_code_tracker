#!/usr/bin/env bash

set -uo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

max_jobs=8
max_log_mb=100
check_interval=2
dry_run=0
timestamp="$(date +"%Y%m%d_%H%M%S")"
log_dir="model_scan/logs_backup_bounded_${timestamp}"

for arg in "$@"; do
  case "$arg" in
    --max-jobs=*)
      max_jobs="${arg#*=}"
      ;;
    --max-log-mb=*)
      max_log_mb="${arg#*=}"
      ;;
    --check-interval=*)
      check_interval="${arg#*=}"
      ;;
    --dry-run)
      dry_run=1
      ;;
    *)
      echo "Unknown option: $arg" >&2
      echo "Usage: bash run_ctau_bkg_model_scan.sh [--max-jobs=N] [--max-log-mb=N] [--check-interval=N] [--dry-run]" >&2
      exit 1
      ;;
  esac
done

for value_name in max_jobs max_log_mb check_interval; do
  value="${!value_name}"
  if ! [[ "$value" =~ ^[1-9][0-9]*$ ]]; then
    echo "[ERROR] ${value_name} must be a positive integer: ${value}" >&2
    exit 1
  fi
done

# Columns: yLow yHigh ptLow ptHigh maxSS maxFL maxDS
# The maxima are read from backup_20.../roots ctau_bkg fit results.
# The scan only lowers complexity: SS=1..maxSS, FL=0..maxFL, DS=0..maxDS.
bin_specs=(
  "0 1.6 20 30 1 0 1"
  "0 1.6 16 20 2 2 1"
  "0 1.6 14 16 2 2 1"
  "0 1.6 12 14 2 2 1"
  "0 1.6 10 12 2 2 1"
  "0 1.6 9 10 2 2 1"
  "0 1.6 8 9 2 2 1"
  "0 1.6 7 8 2 2 1"
  "1.6 2.4 14 20 2 2 1"
  "1.6 2.4 12 14 2 2 1"
  "1.6 2.4 10 12 2 2 1"
  "1.6 2.4 9 10 2 2 1"
  "1.6 2.4 8 9 2 2 1"
  "1.6 2.4 7 8 2 2 1"
  "1.6 2.4 6 7 2 2 1"
  "1.6 2.4 5 6 2 2 1"
  "1.6 2.4 4 5 2 1 1"
  "1.6 2.4 3 4 2 1 1"
  "1.6 2.4 2 3 3 2 1"
  "1.6 2.4 1 2 3 2 1"
)

mkdir -p "$log_dir"
max_bytes=$((max_log_mb * 1024 * 1024))

format_tag() {
  printf '%.2f' "$1"
}

run_one() {
  local y_low="$1" y_high="$2" pt_low="$3" pt_high="$4"
  local ss="$5" fl="$6" ds="$7"
  local ytag="y$(format_tag "$y_low")_$(format_tag "$y_high")"
  local pttag="pt$(format_tag "$pt_low")_$(format_tag "$pt_high")"
  local model_tag="ss${ss}_fl${fl}_ds${ds}"
  local label="${ytag}_${pttag}_${model_tag}"
  local log_file="${log_dir}/fit_${label}.log"
  local macro_call="ctau_bkg.C(${pt_low},${pt_high},${y_low},${y_high},true,false,false,${ss},${fl},${ds},\"${model_tag}\")"

  if (( dry_run )); then
    echo "root -l -b -q '${macro_call}'"
    return 0
  fi

  root -l -b -q "$macro_call" >"$log_file" 2>&1 &
  local root_pid=$!
  local killed_for_size=0

  while kill -0 "$root_pid" 2>/dev/null; do
    if [[ -f "$log_file" ]]; then
      local size
      size="$(stat -c '%s' "$log_file" 2>/dev/null || echo 0)"
      if [[ "$size" =~ ^[0-9]+$ ]] && (( size > max_bytes )); then
        kill -TERM "$root_pid" 2>/dev/null || true
        killed_for_size=1
        break
      fi
    fi
    sleep "$check_interval"
  done

  local exit_code=0
  wait "$root_pid" || exit_code=$?
  if (( killed_for_size )); then
    exit_code=99
  fi

  printf '[RUNNER_STATUS] label=%s exitCode=%d killedForLogSize=%d log=%s\n' \
    "$label" "$exit_code" "$killed_for_size" "$log_file" >>"$log_file"

  if (( exit_code == 0 )); then
    echo "[DONE] $label"
  else
    echo "[FAIL] $label exit=${exit_code} log=${log_file}"
  fi
}

running=0
total=0
for spec in "${bin_specs[@]}"; do
  read -r y_low y_high pt_low pt_high max_ss max_fl max_ds <<<"$spec"
  total=$((total + max_ss * (max_fl + 1) * (max_ds + 1)))
  for ((ss=1; ss<=max_ss; ++ss)); do
    for ((fl=0; fl<=max_fl; ++fl)); do
      for ((ds=0; ds<=max_ds; ++ds)); do
        run_one "$y_low" "$y_high" "$pt_low" "$pt_high" "$ss" "$fl" "$ds" &
        running=$((running + 1))
        if (( running >= max_jobs )); then
          wait -n || true
          running=$((running - 1))
        fi
      done
    done
  done
done
wait || true

if (( dry_run )); then
  echo "[DRY-RUN] ${total} fits listed; no files created beyond ${log_dir}"
  exit 0
fi

summary_file="${log_dir}/model_scan_summary.log"
{
  echo "# ctau_bkg backup-bounded lower-complexity scan"
  echo "# total_expected=${total} max_jobs=${max_jobs} generated=${timestamp}"
  echo "# Send this file to Codex after all fits finish."
  for log_file in "${log_dir}"/fit_*.log; do
    result_line="$(sed -n '/^\[MODEL_SCAN_RESULT\]/p' "$log_file")"
    status_line="$(sed -n '/^\[RUNNER_STATUS\]/p' "$log_file")"
    if [[ -n "$result_line" ]]; then
      echo "$result_line"
    else
      echo "[MODEL_SCAN_MISSING] ${status_line}"
    fi
  done
} >"$summary_file"

echo "[COMPLETE] expected fits: ${total}"
echo "[COMPLETE] detailed logs: ${log_dir}"
echo "[COMPLETE] send this summary: ${summary_file}"
