#!/usr/bin/env python3
"""Rewrite config/pr_gap_ratios.txt from the MC PR prefits in roots/ (q3=(s3-s2)/s1, q4=(s4-s3)/s1).

  python3 study/gap_ratios.py [BASE]    (BASE defaults to this folder)
"""
from pathlib import Path
import csv, subprocess, sys
BASE = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
ROOT = str(Path.home() / 'build/root-6.36.06/bin/root')
rows = []
for r in csv.DictReader(open(BASE / 'config/prefit_settings.txt'), delimiter='\t'):
    reg, lo, hi = r['bin'].split('_')
    y = 'y1.60_2.40' if reg == 'fwd' else 'y0.00_1.60'
    f = BASE / 'roots' / y / 'ctau_pr' / f'ctau_resolution_{y}_pt{float(lo):.2f}_{float(hi):.2f}.root'
    out = subprocess.run([ROOT, '-l', '-b', '-q', '-e', f'TFile f("{f}"); int n=f.Get<TParameter<int>>("nResolutionComponents")->GetVal(); printf("N %d",n); for(int i=1;i<=4;++i){{auto*p=f.Get<TParameter<double>>(Form("ctauTime%dScale",i)); printf(" %.9g", p?p->GetVal():0.);}} printf("\\n");'], capture_output=True, text=True).stdout
    v = next(l.split()[1:] for l in out.splitlines() if l.startswith('N '))
    n, s = int(v[0]), [float(x) for x in v[1:]]
    q3 = (s[2] - s[1]) / s[0] if n >= 3 else 0.0
    q4 = (s[3] - s[2]) / s[0] if n >= 4 else 0.0
    rows.append((r['bin'], n, s[0], q3, q4))
with open(BASE / 'config/pr_gap_ratios.txt', 'w') as fh:
    fh.write('bin\tnres\ts1\tq3\tq4\n')
    for b, n, s1, q3, q4 in rows: fh.write(f'{b}\t{n}\t{s1:.6f}\t{q3:.6f}\t{q4:.6f}\n')
print('wrote', BASE / 'config/pr_gap_ratios.txt')
