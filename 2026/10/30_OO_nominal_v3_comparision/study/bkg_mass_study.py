#!/usr/bin/env python3
"""Background-mass study (2026-10-31): how much B moves per bin when the dimuon mass treatment changes,
with the final-fit rule fixed (default FINAL_V2). Variants per bin, each a prefit_variant.py chain:
  mlo2.7 / mlo2.8 : mass window lower edge 2.6 -> 2.7 / 2.8 GeV (whole chain refitted)
  mhi3.4          : mass window upper edge 3.5 -> 3.4 GeV
  cheb<n>         : mass background Chebychev order n (1-4, except the production order)

  python3 study/bkg_mass_study.py [--jobs 6] [--bins fwd_1_2,...] [--rule FINAL_V2]
  python3 study/prefit_variant.py --summary     (results; vs nominal ver3)
"""
from concurrent.futures import ThreadPoolExecutor
import argparse, csv, subprocess, sys
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser(); ap.add_argument('--jobs', type=int, default=6); ap.add_argument('--bins', default=''); ap.add_argument('--rule', default='FINAL_V2')
ap.add_argument('--variants', default='mlo2.7,mlo2.8,mhi3.4,cheb')
a = ap.parse_args()
rows = {r['bin']: r for r in csv.DictReader(open(BASE / 'config/prefit_settings.txt'), delimiter='\t')}
bins = [b for b in a.bins.split(',') if b] or list(rows)
tasks = []
for b in bins:
    for v in a.variants.split(','):
        if v == 'mlo2.7': tasks.append(('mlo2.7', b, ['mass_lo=2.7']))
        elif v == 'mlo2.8': tasks.append(('mlo2.8', b, ['mass_lo=2.8']))
        elif v == 'mhi3.4': tasks.append(('mhi3.4', b, ['mass_hi=3.4']))
        elif v == 'cheb':
            for n in (1, 2, 3, 4):
                if str(n) != rows[b]['mass_cheb']: tasks.append((f'cheb{n}', b, [f'mass_cheb={n}']))
# heavy bins first
order = {'fwd_1_2': 0, 'fwd_2_3': 1, 'fwd_3_4': 2, 'fwd_4_5': 3}
tasks.sort(key=lambda t: order.get(t[1], 9))
print(f'{len(tasks)} variant chains, {a.jobs} parallel', flush=True)


def one(t):
    name, b, kv = t
    log = BASE / 'logs' / f'bkgmass_{name}_{b}.log'
    with open(log, 'w') as fh:
        rc = subprocess.run(['python3', str(BASE / 'study/prefit_variant.py'), name, b, *kv, f'--rule={a.rule}'], cwd=BASE, stdout=fh, stderr=subprocess.STDOUT).returncode
    reg, lo, hi = b.split('_'); y = 'y1.60_2.40' if reg == 'fwd' else 'y0.00_1.60'
    res = BASE / 'study/prefit_variants' / f'{name}_{b}' / 'roots' / y / 'fit2d' / f'fit2d_result_{y}_pt{float(lo):.2f}_{float(hi):.2f}.root'
    print(f'{name:<8} {b:<10} ' + ('ok' if rc == 0 and res.exists() else f'FAILED rc={rc} result={res.exists()} ({log.relative_to(BASE)})'), flush=True)


with ThreadPoolExecutor(a.jobs) as pool:
    list(pool.map(one, tasks))
print('done; python3 study/prefit_variant.py --summary')
