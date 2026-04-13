#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
timestamp="$(date +"%Y%m%d_%H%M%S")"
log_dir="${script_dir}/logs/eff_${timestamp}"
table_dir="${script_dir}/outputs/tables"
mkdir -p "${log_dir}"
mkdir -p "${table_dir}"

eff_n_evt="${EFF_N_EVT:--1}"
acc_n_evt="${ACC_N_EVT:--1}"

echo "[INFO] logs will be saved in: ${log_dir}"
echo "[INFO] efficiency nEvt = ${eff_n_evt}"
echo "[INFO] acceptance nEvt = ${acc_n_evt}"

required_macros=(
  "efficiency_1d.C"
  "efficiency_1d_pp.C"
  "acceptance_1d.C"
  "extract_efficiency_table.C"
  "extract_efficiency_table_pp.C"
  "extract_acceptance_table.C"
  "extract_acceptance_table_pp.C"
  "draw_acc.C"
  "draw_eff.C"
  "draw_eff_pp.C"
  "draw_eff_pt_cent.C"
  "draw_eff_pt_cent_ptw_compare.C"
)

for macro in "${required_macros[@]}"; do
  if [[ ! -f "${script_dir}/${macro}" ]]; then
    echo "[ERROR] required macro not found: ${script_dir}/${macro}"
    exit 1
  fi
done

run_root_macro() {
  local label="$1"
  local macro_call="$2"
  local log_path="$3"

  echo "[RUN] ${label} -> ${log_path}"
  (
    cd "${script_dir}"
    root -l -b -q "${macro_call}"
  ) > "${log_path}" 2>&1
}

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

pbpb_pr_ptw1_log="${log_dir}/eff_pbpb_pr_ptw1.log"
pbpb_np_ptw1_log="${log_dir}/eff_pbpb_np_ptw1.log"
pbpb_pr_ptw0_log="${log_dir}/eff_pbpb_pr_ptw0.log"
pbpb_np_ptw0_log="${log_dir}/eff_pbpb_np_ptw0.log"
pp_pr_ptw1_log="${log_dir}/eff_pp_pr_ptw1.log"
pp_np_ptw1_log="${log_dir}/eff_pp_np_ptw1.log"
pp_pr_ptw0_log="${log_dir}/eff_pp_pr_ptw0.log"
pp_np_ptw0_log="${log_dir}/eff_pp_np_ptw0.log"

acc_pbpb_pr_log="${log_dir}/acc_pbpb_pr.log"
acc_pbpb_np_log="${log_dir}/acc_pbpb_np.log"
acc_pp_pr_log="${log_dir}/acc_pp_pr.log"
acc_pp_np_log="${log_dir}/acc_pp_np.log"

draw_pbpb_pt_cent_log="${log_dir}/draw_eff_pt_cent.log"
draw_pbpb_ptw_compare_log="${log_dir}/draw_eff_pt_cent_ptw_compare.log"
draw_pbpb_summary_log="${log_dir}/draw_eff.log"
draw_pp_log="${log_dir}/draw_eff_pp.log"
draw_acc_pbpb_log="${log_dir}/draw_acc_pbpb.log"
draw_acc_pp_log="${log_dir}/draw_acc_pp.log"

table_pbpb_pr_log="${log_dir}/extract_eff_table_pbpb_pr.log"
table_pbpb_np_log="${log_dir}/extract_eff_table_pbpb_np.log"
table_pp_pr_log="${log_dir}/extract_eff_table_pp_pr.log"
table_pp_np_log="${log_dir}/extract_eff_table_pp_np.log"

acc_table_pbpb_pr_log="${log_dir}/extract_acc_table_pbpb_pr.log"
acc_table_pbpb_np_log="${log_dir}/extract_acc_table_pbpb_np.log"
acc_table_pp_pr_log="${log_dir}/extract_acc_table_pp_pr.log"
acc_table_pp_np_log="${log_dir}/extract_acc_table_pp_np.log"

# PbPb efficiency
run_root_macro_bg "PbPb prompt efficiency (ptW on)" \
  "efficiency_1d.C(${eff_n_evt},true,true,true,true,true)" \
  "${pbpb_pr_ptw1_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "PbPb nonprompt efficiency (ptW on)" \
  "efficiency_1d.C(${eff_n_evt},false,true,true,true,true)" \
  "${pbpb_np_ptw1_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "PbPb prompt efficiency (ptW off)" \
  "efficiency_1d.C(${eff_n_evt},true,true,true,true,false)" \
  "${pbpb_pr_ptw0_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "PbPb nonprompt efficiency (ptW off)" \
  "efficiency_1d.C(${eff_n_evt},false,true,true,true,false)" \
  "${pbpb_np_ptw0_log}"
register_bg_job "${LAST_PID}"

# pp efficiency
run_root_macro_bg "pp prompt efficiency (ptW on)" \
  "efficiency_1d_pp.C(${eff_n_evt},true,true,true,true,false,true)" \
  "${pp_pr_ptw1_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "pp nonprompt efficiency (ptW on)" \
  "efficiency_1d_pp.C(${eff_n_evt},false,true,true,true,false,true)" \
  "${pp_np_ptw1_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "pp prompt efficiency (ptW off)" \
  "efficiency_1d_pp.C(${eff_n_evt},true,true,true,false,false,true)" \
  "${pp_pr_ptw0_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "pp nonprompt efficiency (ptW off)" \
  "efficiency_1d_pp.C(${eff_n_evt},false,true,true,false,false,true)" \
  "${pp_np_ptw0_log}"
register_bg_job "${LAST_PID}"

if ! wait_for_jobs; then
  echo "[ERROR] efficiency job failed. Check logs in ${log_dir}"
  exit 1
fi

# Acceptance
run_root_macro_bg "PbPb prompt acceptance" \
  "acceptance_1d.C(${acc_n_evt},true,true,true,true)" \
  "${acc_pbpb_pr_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "PbPb nonprompt acceptance" \
  "acceptance_1d.C(${acc_n_evt},false,true,true,true)" \
  "${acc_pbpb_np_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "pp prompt acceptance" \
  "acceptance_1d.C(${acc_n_evt},true,false,true,true)" \
  "${acc_pp_pr_log}"
register_bg_job "${LAST_PID}"

run_root_macro_bg "pp nonprompt acceptance" \
  "acceptance_1d.C(${acc_n_evt},false,false,true,true)" \
  "${acc_pp_np_log}"
register_bg_job "${LAST_PID}"

if ! wait_for_jobs; then
  echo "[ERROR] acceptance job failed. Check logs in ${log_dir}"
  exit 1
fi

# Drawing
run_root_macro "draw_eff_pt_cent" \
  "draw_eff_pt_cent.C(true,true,true,true)" \
  "${draw_pbpb_pt_cent_log}"

run_root_macro "draw_eff_pt_cent_ptw_compare" \
  "draw_eff_pt_cent_ptw_compare.C(true,true,true)" \
  "${draw_pbpb_ptw_compare_log}"

run_root_macro "draw_eff" \
  "draw_eff.C(true,true,true,true)" \
  "${draw_pbpb_summary_log}"

run_root_macro "draw_eff_pp" \
  "draw_eff_pp.C(true,true,true,true)" \
  "${draw_pp_log}"

run_root_macro "draw_acc_pbpb" \
  "draw_acc.C(true,true,true)" \
  "${draw_acc_pbpb_log}"

run_root_macro "draw_acc_pp" \
  "draw_acc.C(false,true,true)" \
  "${draw_acc_pp_log}"

# Efficiency tables
run_root_macro "extract_efficiency_table_pbpb_pr" \
  "extract_efficiency_table.C(\"${script_dir}/skim_roots/eff_PbPb2018_isMC1_PR_ncollW1_genW1_ptW1_tnpW1_Dimuon_MiniAOD.root\",\"${table_dir}/efficiency_table_pbpb_PR.tsv\",false)" \
  "${table_pbpb_pr_log}"

run_root_macro "extract_efficiency_table_pbpb_np" \
  "extract_efficiency_table.C(\"${script_dir}/skim_roots/eff_PbPb2018_isMC1_NP_ncollW1_genW1_ptW1_tnpW1_Dimuon_MiniAOD.root\",\"${table_dir}/efficiency_table_pbpb_NP.tsv\",false)" \
  "${table_pbpb_np_log}"

run_root_macro "extract_efficiency_table_pp_pr" \
  "extract_efficiency_table_pp.C(\"${script_dir}/skim_roots/eff_pp5p02TeV_isMC1_PR_ncollW1_genW1_ptW1_tnpW1_Dimuon_MiniAOD.root\",\"${table_dir}/efficiency_table_pp_PR.tsv\",false)" \
  "${table_pp_pr_log}"

run_root_macro "extract_efficiency_table_pp_np" \
  "extract_efficiency_table_pp.C(\"${script_dir}/skim_roots/eff_pp5p02TeV_isMC1_NP_ncollW1_genW1_ptW1_tnpW1_Dimuon_MiniAOD.root\",\"${table_dir}/efficiency_table_pp_NP.tsv\",false)" \
  "${table_pp_np_log}"

# Acceptance tables
run_root_macro "extract_acceptance_table_pbpb_pr" \
  "extract_acceptance_table.C(\"${script_dir}/skim_roots/acc_PbPb2018_ppInput_isMC1_PR_ncollW0_genW1_ptW1.root\",\"${table_dir}/acceptance_table_pbpb_PR.tsv\",false)" \
  "${acc_table_pbpb_pr_log}"

run_root_macro "extract_acceptance_table_pbpb_np" \
  "extract_acceptance_table.C(\"${script_dir}/skim_roots/acc_PbPb2018_ppInput_isMC1_NP_ncollW0_genW1_ptW1.root\",\"${table_dir}/acceptance_table_pbpb_NP.tsv\",false)" \
  "${acc_table_pbpb_np_log}"

run_root_macro "extract_acceptance_table_pp_pr" \
  "extract_acceptance_table_pp.C(\"${script_dir}/skim_roots/acc_pp2018_ppInput_isMC1_PR_ncollW0_genW1_ptW1.root\",\"${table_dir}/acceptance_table_pp_PR.tsv\",false)" \
  "${acc_table_pp_pr_log}"

run_root_macro "extract_acceptance_table_pp_np" \
  "extract_acceptance_table_pp.C(\"${script_dir}/skim_roots/acc_pp2018_ppInput_isMC1_NP_ncollW0_genW1_ptW1.root\",\"${table_dir}/acceptance_table_pp_NP.tsv\",false)" \
  "${acc_table_pp_np_log}"

echo "[INFO] completed successfully"
echo "[INFO] PbPb PR efficiency log: ${pbpb_pr_ptw1_log}"
echo "[INFO] PbPb NP efficiency log: ${pbpb_np_ptw1_log}"
echo "[INFO] pp PR efficiency log: ${pp_pr_ptw1_log}"
echo "[INFO] pp NP efficiency log: ${pp_np_ptw1_log}"
echo "[INFO] PbPb PR acceptance log: ${acc_pbpb_pr_log}"
echo "[INFO] PbPb NP acceptance log: ${acc_pbpb_np_log}"
echo "[INFO] pp PR acceptance log: ${acc_pp_pr_log}"
echo "[INFO] pp NP acceptance log: ${acc_pp_np_log}"
echo "[INFO] PbPb draw log: ${draw_pbpb_summary_log}"
echo "[INFO] pp draw log: ${draw_pp_log}"
echo "[INFO] PbPb acceptance draw log: ${draw_acc_pbpb_log}"
echo "[INFO] pp acceptance draw log: ${draw_acc_pp_log}"
echo "[INFO] PbPb PR efficiency table log: ${table_pbpb_pr_log}"
echo "[INFO] PbPb NP efficiency table log: ${table_pbpb_np_log}"
echo "[INFO] pp PR efficiency table log: ${table_pp_pr_log}"
echo "[INFO] pp NP efficiency table log: ${table_pp_np_log}"
echo "[INFO] PbPb PR acceptance table log: ${acc_table_pbpb_pr_log}"
echo "[INFO] PbPb NP acceptance table log: ${acc_table_pbpb_np_log}"
echo "[INFO] pp PR acceptance table log: ${acc_table_pp_pr_log}"
echo "[INFO] pp NP acceptance table log: ${acc_table_pp_np_log}"
