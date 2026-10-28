#!/bin/bash
# Follow a scan log; every fit within +-2% gets scripts/candidate_check.py. Output: $2 (one line per fit).
# usage: scripts/watch_scan.sh logs/scan_S3.log docs/evidence/fwd1012_scan_S3_candidates.txt scans/S3_fwd_10_12
log=$1; out=$2; top=$3; touch $out
while true; do
  grep "%" $log | awk '{p=$3; sub("%","",p); if (p+0>-2.0 && p+0<2.0) print $1, $2, $3}' | while read c r p; do
    grep -q "^$c $r " $out && continue
    lim=$(python3 scripts/candidate_check.py $top/$c $r 2>/dev/null)
    ch=$(grep "^$c .*$r " $log | head -1 | sed 's/.*worst/worst/')
    echo "$c $r $p $lim $ch" >> $out
  done
  grep -q "^done" $log && break
  sleep 60
done
