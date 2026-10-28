#!/usr/bin/env python3
"""Re-score every B value recorded in ../25_OO_sept_acc_general_rule/docs/evidence/*.txt against nominal ver2.
Writes docs/evidence/old_candidates_vs_v2.txt (per bin, sorted by |%|). These were fitted on the prefits of their time;
anything used is refitted with run_v2.py / prefit_v2.py before adoption."""
from pathlib import Path
import re, sys
HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / 'scripts'))
from nominal_v2 import NOM
P1 = re.compile(r'^(\S+)\s+(?:(\S+)\s+)?((?:fwd|mid)_\d+_\d+)\s+B=([\d.]+)±([\d.]+).*?(st\d.*)$')
P2 = re.compile(r'^(\S+?)_((?:fwd|mid)_\d+_\d+)\s+(\S+)\s+B=([\d.]+)±([\d.]+).*?(st\d.*)$')
rows = set()
for f in sorted((HERE.parent / '25_OO_sept_acc_general_rule/docs/evidence').glob('*.txt')):
    for l in open(f):
        if m := P1.match(l):
            a, x, b, B, e, st = m.groups(); rows.add((b, a + ('|' + x if x else ''), float(B), float(e), f.stem, st.rsplit(' ', 1)[0]))
        elif m := P2.match(l):
            a, b, x, B, e, st = m.groups(); rows.add((b, a + '|' + x, float(B), float(e), f.stem, st.rsplit(' ', 1)[0]))
out = []
for b, (v, ve) in NOM.items():
    rs = sorted((r for r in rows if r[0] == b), key=lambda r: abs(r[2] / v - 1))
    out.append(f'== {b}  ver2 {v:.5f} ± {ve:.5f}  ({len(rs)} recorded fits)')
    out += [f'  {100 * (r[2] / v - 1):+6.2f}% {(r[2] - v) / ve:+5.2f}σ  B={r[2]:.4f}±{r[3]:.4f}  {r[1]:<48} {r[4]:<28} {r[5]}' for r in rs]
(HERE / 'docs/evidence/old_candidates_vs_v2.txt').write_text('\n'.join(out) + '\n')
print(len(rows), 'rows')
