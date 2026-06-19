#!/bin/bash

# how to run
# N_EVT=-1 USE_GLOBAL=1 ./run_skim.sh
# N_EVT=-1 USE_GLOBAL=0 ./run_skim.sh

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
timestamp="$(date +"%Y%m%d_%H%M%S")"
log_dir="${script_dir}/logs/skim_${timestamp}"
mkdir -p "${log_dir}"
mkdir -p "${script_dir}/skim_roots"
mkdir -p "${script_dir}/roodataset_roots"
mkdir -p "${script_dir}/figs"

n_evt="${N_EVT:--1}"
use_global="${USE_GLOBAL:-1}"

if [[ "${use_global}" != "0" && "${use_global}" != "1" ]]; then
  echo "[ERROR] USE_GLOBAL must be 0 or 1"
  exit 1
fi

echo "[INFO] logs will be saved in: ${log_dir}"
echo "[INFO] nEvt = ${n_evt}"
echo "[INFO] useGlobal = ${use_global}"

required_macros=(
  "onia_to_skim_vector.C"
  "plot_muon_distribution.C"
  "skim_to_roodataset.C"
)

for macro in "${required_macros[@]}"; do
  if [[ ! -f "${script_dir}/${macro}" ]]; then
    echo "[ERROR] required macro not found: ${script_dir}/${macro}"
    exit 1
  fi
done

run_root_macro_bg() {
  local label="$1"
  local macro_call="$2"
  local log_path="$3"

  echo "[RUN] ${label} -> ${log_path}" >&2
  (
    cd "${script_dir}"
    root -l -b -q "${macro_call}"
  ) > "${log_path}" 2>&1 &
  LAST_PID=$!
}

declare -a pids=()

register_bg_job() {
  local pid="$1"
  pids+=("${pid}")
}

wait_for_jobs() {
  local status=0
  for pid in "${pids[@]}"; do
    wait "${pid}" || status=$?
  done
  pids=()
  return "${status}"
}

data_log="${log_dir}/skim_data.log"
pr_log="${log_dir}/skim_mc_pr.log"
np_log="${log_dir}/skim_mc_np.log"
plot_data_log="${log_dir}/plot_data.log"
plot_pr_log="${log_dir}/plot_mc_pr.log"
plot_np_log="${log_dir}/plot_mc_np.log"
roo_data_log="${log_dir}/roodataset_data.log"
roo_pr_log="${log_dir}/roodataset_mc_pr.log"
roo_np_log="${log_dir}/roodataset_mc_np.log"

run_root_macro_bg "OO data skim" \
  "onia_to_skim_vector.C(${n_evt},false,false,${use_global})" \
  "${data_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "OO prompt MC skim" \
  "onia_to_skim_vector.C(${n_evt},true,true,${use_global})" \
  "${pr_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "OO nonprompt MC skim" \
  "onia_to_skim_vector.C(${n_evt},true,false,${use_global})" \
  "${np_log}"
register_bg_job "${LAST_PID}"

if ! wait_for_jobs; then
  echo "[ERROR] skim job failed. Check logs in ${log_dir}"
  exit 1
fi

run_root_macro_bg "OO data muon plot" \
  "plot_muon_distribution.C(false,false,${use_global})" \
  "${plot_data_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "OO prompt MC muon plot" \
  "plot_muon_distribution.C(true,true,${use_global})" \
  "${plot_pr_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "OO nonprompt MC muon plot" \
  "plot_muon_distribution.C(true,false,${use_global})" \
  "${plot_np_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "OO data RooDataSet" \
  "skim_to_roodataset.C(${n_evt},false,0,${use_global},\"_OO25_HLT_OxyL1SingleMuOpen_v1\",false,false,false,false,0,200,2.6,3.5,true,0)" \
  "${roo_data_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "OO prompt MC RooDataSet" \
  "skim_to_roodataset.C(${n_evt},true,1,${use_global},\"_OO25_HLT_OxyL1SingleMuOpen_v1\",false,false,false,false,0,200,2.6,3.5,true,0)" \
  "${roo_pr_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "OO nonprompt MC RooDataSet" \
  "skim_to_roodataset.C(${n_evt},true,2,${use_global},\"_OO25_HLT_OxyL1SingleMuOpen_v1\",false,false,false,false,0,200,2.6,3.5,true,0)" \
  "${roo_np_log}"
register_bg_job "${LAST_PID}"

if ! wait_for_jobs; then
  echo "[ERROR] plot or RooDataSet job failed. Check logs in ${log_dir}"
  exit 1
fi

echo "[INFO] all jobs finished successfully"
echo "[INFO] skim logs: ${data_log}, ${pr_log}, ${np_log}"
echo "[INFO] plot logs: ${plot_data_log}, ${plot_pr_log}, ${plot_np_log}"
echo "[INFO] RooDataSet logs: ${roo_data_log}, ${roo_pr_log}, ${roo_np_log}"
