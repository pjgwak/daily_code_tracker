#!/usr/bin/env python3
"""Summarize the final-fit treatment (fixed / constrained / free) per PDF group.

  python3 study/treatment.py   (reads roots/*/fit2d, writes drawing/final_summary/treatment.{csv,md})
Group names follow the prefit that supplies each PDF.
"""
import csv, subprocess, sys
from collections import defaultdict
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
ROOT = str(Path.home() / 'build/root-6.36.06/bin/root')

def group(p):
    if p in ('bFraction', 'Nsig', 'Nbkg'): return 'yield/B'
    if p.startswith('signal_mass'): return 'sig mass (MC mass prefit)'
    if p.startswith('bkg_mass'): return 'bkg mass (DATA mass prefit)'
    if p.startswith('ctauTime') or p.startswith('ctauRecursiveFrac'): return 'PR (MC PR prefit)'
    if p.startswith('signal') or p.startswith('NsignalSS') or p.startswith('npComposition'): return 'NP (MC NP prefit)'
    if p.startswith('bkg') or p[:3] in ('Nss', 'Nds', 'Nfl', 'Npl'): return 'bkg ctau (sideband prefit)'
    return 'other'

def main(rule='FINAL_v3'):
    out = BASE / 'drawing' / 'final_summary'
    subprocess.run([ROOT, '-l', '-b', '-q', f'study/dump_treatment.C("{BASE}","{out / "treatment.csv"}")'], cwd=BASE, check=True, stdout=subprocess.DEVNULL)
    rows = list(csv.DictReader(open(out / 'treatment.csv')))
    table = defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))
    for r in rows:
        table[r['bin']][group(r['parameter'])][['fixed', 'constrained', 'free'].index(r['treatment'])] += 1
    groups = ['sig mass (MC mass prefit)', 'bkg mass (DATA mass prefit)', 'PR (MC PR prefit)', 'NP (MC NP prefit)', 'bkg ctau (sideband prefit)', 'yield/B', 'other']
    lines = ['| bin | ' + ' | '.join(groups) + ' |', '|---' * (len(groups) + 1) + '|']
    for b in table:
        lines.append(f'| {b} | ' + ' | '.join('/'.join(map(str, table[b][g])) for g in groups) + ' |')
    text = f'고정/제약/자유 개수 ({rule})\n\n' + '\n'.join(lines) + '\n'
    (out / 'treatment.md').write_text(text); print(text)

if __name__ == '__main__':
    main()
