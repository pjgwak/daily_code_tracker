#!/usr/bin/env python3
"""Build a collect.C-style summary table for a per-bin choice file from the scan logs (no refit).

  python3 study/choice_summary.py config/final_v3_choice.txt docs/evidence/summary_FINAL_V3_draft.txt
"""
import sys, re
sys.path.insert(0, 'study'); import v3_summary as S
res = S.load()
text = open('utils/nominal_bfraction.h').read(); fwd, mid = text.split('midBins()', 1); NOM = {}
for reg, part in (('fwd', fwd), ('mid', mid)):
    for lo, hi, v, e in re.findall(r'\{\s*([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\d.]+)', part): NOM[f'{reg}_{lo}_{hi}'] = (float(v), float(e))
rows = [l.rstrip('\n').split('\t') for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
with open(sys.argv[2], 'w') as o:
    o.write('region\tptLow\tptHigh\tB\tBerr\tBnom\tBnomErr\tpull\tpercent\trule\n')
    order = {'fwd': 0, 'mid': 1}
    for b, r, *_ in sorted(rows, key=lambda x: (order[x[0].split('_')[0]], float(x[0].split('_')[1]))):
        v = res.get((r, b)); n, ne = NOM[b]
        if not v: print('missing', b, r); continue
        reg, lo, hi = b.split('_')
        o.write(f'{reg}\t{lo}\t{hi}\t{v["B"]:.6f}\t{v["Be"]:.6f}\t{n}\t{ne}\t{(v["B"]-n)/ne:.4f}\t{100*(v["B"]/n-1):.3f}\t{r}\n')
print('wrote', sys.argv[2])
