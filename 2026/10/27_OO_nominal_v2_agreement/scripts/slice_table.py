#!/usr/bin/env python3
"""Summarise an OO_SLICE_DUMP csv: per mass slice x ctau region, data / model (PR, NP, bkg) and (data-model)/sqrt(model)."""
import csv, sys, collections
rows = list(csv.DictReader(open(sys.argv[1])))
g = collections.OrderedDict()
for r in rows: g.setdefault((r['m_lo'], r['m_hi']), []).append(r)
regions = [('t<0', lambda a, b: b <= 0), ('0-0.1', lambda a, b: a >= 0 and b <= 0.1), ('0.1-0.3', lambda a, b: a >= 0.1 and b <= 0.3),
           ('0.3-1', lambda a, b: a >= 0.3 and b <= 1.0), ('>1', lambda a, b: a >= 1.0)]
print('mass slice'.ljust(14) + ''.join(f'{n:>26}' for n, _ in regions))
chi = 0
for k, rs in g.items():
    line = f'{k[0]}-{k[1]}'.ljust(14)
    for n, f in regions:
        sel = [r for r in rs if f(float(r['t_lo']), float(r['t_hi']))]
        d = sum(float(r['data']) for r in sel); m = sum(float(r['pr']) + float(r['np']) + float(r['bkg']) for r in sel)
        b = sum(float(r['bkg']) for r in sel); p = (d - m) / max(m, 1e-9) ** .5; chi += p * p
        line += f'{d:6.0f}/{m:6.1f} bkg{b:5.0f} {p:+5.1f}'
    print(line)
print(f'sum pull^2 = {chi:.1f} over {len(g) * len(regions)} cells')
