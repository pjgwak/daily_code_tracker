#!/usr/bin/env python3
"""Limit status of one scan combo x rule: AT-LIMIT flags of every prefit stage (prefit_audit.txt) that production (c0) does not
already have, plus parameters at a bound in the background prefit and the final fit (limit_check.py).
Accepted as in production: signal_mass_frac_ratio1, the degenerate 2nd NP component (NsignalSS2, signal2 lifetime).
  python3 scripts/candidate_check.py COMBO_DIR RULE    -> 'CLEAN' or comma list"""
import subprocess, sys
from pathlib import Path
ACCEPT = {'signal_mass_frac_ratio1', 'NsignalSS2', 'signal2_log_lifetime_offset'}
w, rule = Path(sys.argv[1]), sys.argv[2]
bad = set()
a = w / 'prefit_audit.txt'
if a.exists():
    for l in a.read_text().splitlines():
        if l.strip().startswith('!') and 'AT-LIMIT' in l: bad.add(l.split()[1])
out = subprocess.run(['python3', str(Path(__file__).parent / 'limit_check.py'), str(w), rule], capture_output=True, text=True).stdout
bad |= {l.split()[1] for l in out.splitlines() if 'LIMIT' in l}
bad -= ACCEPT
print(','.join(sorted(bad)) or 'CLEAN')
