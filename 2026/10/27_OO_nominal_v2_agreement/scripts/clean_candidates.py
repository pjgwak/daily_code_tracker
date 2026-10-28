#!/usr/bin/env python3
"""From a scan results.tsv, keep fits with |pct| < CUT, st0/cq3, prefit worst chi2 < 1.5, then run limit_check on the
background prefit and the final fit; print only candidates with no parameter at a bound (NsignalSS2 / signal2 lifetime
of the degenerate 2nd NP component are reported separately since production has them too).
  python3 scripts/clean_candidates.py scans/S2_fwd_10_12 [CUT=2.0]"""
import csv, json, subprocess, sys
from pathlib import Path
top = Path(sys.argv[1]); cut = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
rows = [r for r in csv.DictReader(open(top / 'results.tsv'), delimiter='\t') if r['pct'] and r['status'] == '0' and r['covQual'] == '3'
        and abs(float(r['pct'])) < cut and float(r['prefit_worst_chi2']) < 1.5]
IGN = ('NsignalSS2', 'signal2_log_lifetime_offset', 'NsignalSS3', 'signal3_log_lifetime_offset')
for r in sorted(rows, key=lambda r: abs(float(r['pct']))):
    out = subprocess.run(['python3', str(Path(__file__).parent / 'limit_check.py'), str(top / f"c{r['combo']}"), r['rule']], capture_output=True, text=True).stdout
    lim = [l.split()[1] for l in out.splitlines() if 'LIMIT' in l]
    hard = [x for x in lim if x not in IGN]
    tag = 'CLEAN' if not hard else 'atlimit:' + ','.join(sorted(set(hard)))
    print(f"c{r['combo']:<4} {r['rule']:<22} {float(r['pct']):+6.2f}%  wc2 {float(r['prefit_worst_chi2']):.2f}  {tag:<60} {r['changes']}", flush=True)
