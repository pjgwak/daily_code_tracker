#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

mode="both"
dry_run=0
max_jobs=24 # or 16
max_log_mb=100 # if a size of log file > 50 MB due to error messages, kill the process
# without it the log can take 10 GB!
log_check_interval=2 # log size monitoring period, 2 second
timestamp="$(date +"%Y%m%d_%H%M%S")"
log_dir="logs/run_publish_${timestamp}"

for arg in "$@"; do
  case "$arg" in
    --mc)
      mode="mc"
      ;;
    --mass)
      mode="mass"
      ;;
    --pr)
      mode="pr"
      ;;
    --np)
      mode="np"
      ;;
    --bkg)
      mode="bkg"
      ;;
    --err)
      mode="err"
      ;;
    --2d)
      mode="fit2d"
      ;;
    --dry-run)
      dry_run=1
      ;;
    --max-jobs=*)
      max_jobs="${arg#*=}"
      ;;
    --max-log-mb=*)
      max_log_mb="${arg#*=}"
      ;;
    --log-check-interval=*)
      log_check_interval="${arg#*=}"
      ;;
    *)
      echo "Unknown option: $arg" >&2
      echo "Usage: bash run_publish.sh [--mc|--mass|--pr|--np|--bkg|--err|--2d] [--dry-run] [--max-jobs=N] [--max-log-mb=N] [--log-check-interval=N]" >&2
      exit 1
      ;;
  esac
done

if ! [[ "$max_jobs" =~ ^[0-9]+$ ]] || (( max_jobs < 1 )); then
  echo "[ERROR] --max-jobs must be a positive integer" >&2
  exit 1
fi

if ! [[ "$max_log_mb" =~ ^[0-9]+$ ]] || (( max_log_mb < 1 )); then
  echo "[ERROR] --max-log-mb must be a positive integer" >&2
  exit 1
fi

if ! [[ "$log_check_interval" =~ ^[0-9]+$ ]] || (( log_check_interval < 1 )); then
  echo "[ERROR] --log-check-interval must be a positive integer" >&2
  exit 1
fi

mkdir -p "${log_dir}"
echo "[INFO] logs will be saved in: ${log_dir}"
echo "[INFO] log monitor: max_log_mb=${max_log_mb}, check_interval=${log_check_interval}s"

# Edit these bin groups as needed.
# Format:
#   "yLow yHigh|ptLow1 ptHigh1;ptLow2 ptHigh2;..."
bin_groups=(
  "0 1.6|20 35;16 20;14 16;12 14;10 12;9 10;8 9;7 8"
  "1.6 2.4|14 20;12 14;10 12;9 10;8 9;7 8;6 7;5 6;4 5;3 4;2 3;1 2"
  # "1.6 2.4|14 20;10 12"
)

run_root_macro() {
  local macro_name="$1"
  local pt_low="$2"
  local pt_high="$3"
  local y_low="$4"
  local y_high="$5"
  local args

  case "$macro_name" in
    mc_mass.C)
      args="${pt_low},${pt_high},${y_low},${y_high},false,false,true"
      ;;
    mass.C|ctau_pr.C|ctau_np.C|err2.C)
      args="${pt_low},${pt_high},${y_low},${y_high},false,true"
      ;;
    ctau_bkg.C)
      args="${pt_low},${pt_high},${y_low},${y_high},true,false,true"
      ;;
    fit2d.C)
      args="${pt_low},${pt_high},${y_low},${y_high},false,false,true"
      ;;
    *)
      echo "[ERROR] unsupported macro for publish mode: ${macro_name}" >&2
      return 1
      ;;
  esac

  local cmd="root -b -q '${macro_name}(${args})'"

  echo "$cmd"
  if [[ "$dry_run" -eq 0 ]]; then
    eval "$cmd"
  fi
}

run_bin() {
  local pt_low="$1"
  local pt_high="$2"
  local y_low="$3"
  local y_high="$4"

  echo "============================================================"
  echo "pt: ${pt_low}-${pt_high}, |y|: ${y_low}-${y_high}"
  echo "============================================================"

  if [[ "$mode" == "both" || "$mode" == "mc" ]]; then
    run_root_macro "mc_mass.C" "$pt_low" "$pt_high" "$y_low" "$y_high"
  fi

  if [[ "$mode" == "both" || "$mode" == "mass" ]]; then
    run_root_macro "mass.C" "$pt_low" "$pt_high" "$y_low" "$y_high"
  fi

  if [[ "$mode" == "both" || "$mode" == "pr" ]]; then
    run_root_macro "ctau_pr.C" "$pt_low" "$pt_high" "$y_low" "$y_high"
  fi

  if [[ "$mode" == "both" || "$mode" == "np" ]]; then
    run_root_macro "ctau_np.C" "$pt_low" "$pt_high" "$y_low" "$y_high"
  fi

  if [[ "$mode" == "both" || "$mode" == "err" ]]; then
    run_root_macro "err2.C" "$pt_low" "$pt_high" "$y_low" "$y_high"
  fi

  if [[ "$mode" == "both" || "$mode" == "bkg" ]]; then
    run_root_macro "ctau_bkg.C" "$pt_low" "$pt_high" "$y_low" "$y_high"
  fi

  if [[ "$mode" == "both" || "$mode" == "fit2d" ]]; then
    run_root_macro "fit2d.C" "$pt_low" "$pt_high" "$y_low" "$y_high"
  fi
}

format_tag() {
  local value="$1"
  value="${value//- /}"
  value="${value//-/_}"
  value="${value//./p}"
  echo "$value"
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

start_job() {
  local pt_low="$1"
  local pt_high="$2"
  local y_low="$3"
  local y_high="$4"
  local log_file="$5"
  local label="$6"
  local max_bytes="$7"
  local interval="$8"

  (
    run_bin "$pt_low" "$pt_high" "$y_low" "$y_high"
  ) > "${log_file}" 2>&1 &
  local job_pid=$!

  monitor_log_size "$job_pid" "$log_file" "$label" "$max_bytes" "$interval" &
  local monitor_pid=$!

  wait "$job_pid"
  local job_status=$?

  wait "$monitor_pid" 2>/dev/null || true

  return "$job_status"
}

jobs_started=0
jobs_killed=0
jobs_failed=0
max_log_bytes=$((max_log_mb * 1024 * 1024))

wait_for_running_jobs() {
  local pid
  local status
  local -a pids=()

  mapfile -t pids < <(jobs -pr)
  for pid in "${pids[@]}"; do
    if wait "$pid"; then
      continue
    fi

    status=$?
    jobs_failed=$((jobs_failed + 1))
    echo "[WARN] job pid=${pid} exited with status=${status}" >&2
  done
}

for group in "${bin_groups[@]}"; do
  IFS='|' read -r y_bin pt_bin_list <<< "$group"
  read -r y_low y_high <<< "$y_bin"

  IFS=';' read -r -a pt_bins <<< "$pt_bin_list"
  for pt_bin in "${pt_bins[@]}"; do
    read -r pt_low pt_high <<< "$pt_bin"
    y_low_tag="$(format_tag "$y_low")"
    y_high_tag="$(format_tag "$y_high")"
    pt_low_tag="$(format_tag "$pt_low")"
    pt_high_tag="$(format_tag "$pt_high")"
    log_file="${log_dir}/y${y_low_tag}_${y_high_tag}_pt${pt_low_tag}_${pt_high_tag}.log"
    label="pt ${pt_low}-${pt_high}, |y| ${y_low}-${y_high}"
    echo "[RUN] ${label} -> ${log_file}"
    start_job "$pt_low" "$pt_high" "$y_low" "$y_high" "$log_file" "$label" "$max_log_bytes" "$log_check_interval" &
    jobs_started=$((jobs_started + 1))
    if (( jobs_started % max_jobs == 0 )); then
      wait_for_running_jobs
    fi
  done
done

wait_for_running_jobs

jobs_killed=$(find "${log_dir}" -maxdepth 1 -type f -name '*.log.killed' | wc -l)

summary_file="${log_dir}/summary.log"
{
  echo "[INFO] mode=${mode} dry_run=${dry_run} max_jobs=${max_jobs}"
  echo "[INFO] max_log_mb=${max_log_mb} log_check_interval=${log_check_interval}"
  echo "[INFO] total_bins=${jobs_started}"
  echo "[INFO] killed_bins=${jobs_killed}"
  echo "[INFO] failed_jobs=${jobs_failed}"
  echo "[INFO] generated_logs=${log_dir}"
} > "${summary_file}"
echo "[INFO] summary: ${summary_file}"
