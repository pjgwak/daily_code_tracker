#!/usr/bin/env python3
"""Compose a final rule directory from per-bin choices among scanned rules.

  python3 study/compose_final.py FINAL_V3 config/final_v3_choice.txt

The choice file has lines "<bin>\t<rule>[\t<note>]" (rule = study/rules/<rule>/ directory name, or
"prefit:<variant>" for a bin whose prefits come from study/prefit_variants/<variant>_<bin>/ — the knob file
is then taken from that variant's rule). Writes study/rules/<NAME>/<bin file>.txt with the source noted.
"""
import sys, shutil
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
name, choice = sys.argv[1], Path(sys.argv[2])
out = BASE / 'study/rules' / name
out.mkdir(parents=True, exist_ok=True)
EDGES = {'fwd': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 20], 'mid': [7, 8, 9, 10, 12, 14, 16, 20, 30]}
done = set()
for line in open(choice):
    if not line.strip() or line.startswith('#'): continue
    f = line.rstrip('\n').split('\t'); b, rule = f[0], f[1]; note = f[2] if len(f) > 2 else ''
    reg, lo, hi = b.split('_'); y = 'y1.60_2.40' if reg == 'fwd' else 'y0.00_1.60'
    fn = f'{y}_pt{float(lo):.2f}_{float(hi):.2f}.txt'
    src = BASE / 'study/rules' / rule / fn
    text = src.read_text().splitlines()
    text[0] = f'rule {name}'
    text.insert(1, f'# {name} ({b}): from study/rules/{rule}' + (f' — {note}' if note else ''))
    (out / fn).write_text('\n'.join(text) + '\n'); done.add(b)
missing = [f'{r}_{a}_{c}' for r, e in EDGES.items() for a, c in zip(e, e[1:]) if f'{r}_{a}_{c}' not in done]
print('wrote', out, '| missing bins:', missing or 'none')
