#!/usr/bin/env bash
# Wait for the first bkg_mass_study.py pass, then rerun the chains that failed (build-directory compile race on 10/31).
cd "$(dirname "$0")/.."
while pgrep -f "study/bkg_mass_study.py" >/dev/null; do sleep 60; done
grep FAILED logs/bkgmass_driver.log | awk '{print $1, $2}' | sort -u > logs/bkgmass_failed_first_pass.txt
echo "retrying $(wc -l < logs/bkgmass_failed_first_pass.txt) chains"
declare -A KV=([mlo2.7]="mass_lo=2.7" [mlo2.8]="mass_lo=2.8" [mhi3.4]="mass_hi=3.4" [cheb1]="mass_cheb=1" [cheb2]="mass_cheb=2" [cheb3]="mass_cheb=3" [cheb4]="mass_cheb=4")
run() { name=$1; b=$2; python3 study/prefit_variant.py "$name" "$b" "${KV[$name]}" --rule=FINAL_V2 > "logs/bkgmass_${name}_${b}.log" 2>&1; echo "$name $b rc=$?"; }
export -f run; export -A KV 2>/dev/null
# simple pool of 6
n=0
while read -r name b; do
  ( python3 study/prefit_variant.py "$name" "$b" "${KV[$name]}" --rule=FINAL_V2 > "logs/bkgmass_${name}_${b}.log" 2>&1; echo "$name $b rc=$?" ) &
  n=$((n+1)); if [ $((n % 6)) -eq 0 ]; then wait; fi
done < logs/bkgmass_failed_first_pass.txt
wait; echo "retry done"
