#!/bin/bash
# final-fit ctau low-edge scan on one prefit work dir: ctlow_scan.sh WORKDIR LOW... (rules: tail k 0.25/0.5/1)
w=$1; shift; ulimit -f 2000000
cp -n $w/config/prefit_settings.txt $w/config/prefit_settings_prefitwin.txt
for lo in "$@"; do
  awk -F'\t' -v lo=$lo 'BEGIN{OFS="\t"} $1=="fwd_10_12"{$2=lo}1' $w/config/prefit_settings_prefitwin.txt > $w/config/prefit_settings.txt
  echo "== low $lo"
  python3 scripts/refit_v2.py $w fwd_10_12 F_tail0.25_free_long_free1 F_tail0.5_np1 F_tail1_np1
  for r in F_tail0.25_free_long_free1 F_tail0.5_np1 F_tail1_np1; do rm -rf $w/rules/${r}_lo$lo; cp -r $w/rules/$r $w/rules/${r}_lo$lo; echo "$r lo$lo: $(python3 scripts/candidate_check.py $w ${r} | tail -1)"; done
done
