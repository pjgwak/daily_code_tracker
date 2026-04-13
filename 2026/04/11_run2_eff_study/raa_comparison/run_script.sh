#!/usr/bin/env bash

set -euo pipefail

cd /data/users/pjgwak/work/daily_code_tracker/2026/04/11_run2_eff_study/raa_comparison

root -l -b -q 'compare_pT_Jpsi.C(true)'
root -l -b -q 'compare_Npart_Jpsi.C(true)'
