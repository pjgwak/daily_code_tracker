#!/usr/bin/env python3
"""fwd 2-3 ctau-window study (2026-10-31): refit mass screen + sideband background + final 2D fit for several ctau windows,
with the final-fit rule structure fixed (default: FINAL_V2 fwd 2-3 file; ct_low/ct_high replaced per window).

  python3 study/fwd23_window_study.py [--jobs 8] [--base FINAL_V2] [--prefix W23] WINDOW [WINDOW ...]
  WINDOW = lo:hi[:ss/flip/ds[:wide]]   e.g. -6:10:1/1/1   -20:20:2/0/1:1   (default composition: the production one)

Each window becomes a rule directory study/rules/<prefix>_w<lo>_<hi>[_b<ss><flip><ds>[w]] and a prefit_variant.py chain
(work dir study/prefit_variants/<name>_fwd_2_3, log logs/<name>_fwd_2_3.log).  Results: python3 study/prefit_variant.py --summary
"""
from concurrent.futures import ThreadPoolExecutor
import argparse, csv, subprocess, sys
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser(); ap.add_argument('windows', nargs='+'); ap.add_argument('--jobs', type=int, default=8)
ap.add_argument('--base', default='FINAL_V2'); ap.add_argument('--prefix', default='W23'); ap.add_argument('--bin', default='fwd_2_3')
a = ap.parse_args()
reg, lo_, hi_ = a.bin.split('_'); y = 'y1.60_2.40' if reg == 'fwd' else 'y0.00_1.60'
tag = f'{y}_pt{float(lo_):.2f}_{float(hi_):.2f}'
src = (BASE / 'study/rules' / a.base / f'{tag}.txt').read_text().splitlines()
body = [l for l in src if l and not l.startswith(('rule', '#', 'ct_low', 'ct_high'))]
row = {r['bin']: r for r in csv.DictReader(open(BASE / 'config/prefit_settings.txt'), delimiter='\t')}[a.bin]
tasks = []
for w in a.windows:
    parts = w.split(':'); lo, hi = parts[0], parts[1]
    name = f'{a.prefix}_w{lo}_{hi}'.replace('-', 'm').replace('.', 'p')
    kv = [f'ct_low={lo}', f'ct_high={hi}', f'bkg_saved_ct={lo},{hi}']
    if len(parts) > 2 and parts[2]:
        ss, fl, ds = parts[2].split('/'); wide = parts[3] if len(parts) > 3 else row.get('bkg_wide_tau', '0')
        kv += [f'bkg_ss={ss}', f'bkg_flip={fl}', f'bkg_ds={ds}', f'bkg_wide_tau={wide}']
        name += f'_b{ss}{fl}{ds}' + ('w' if wide == '1' else '')
    d = BASE / 'study/rules' / name; d.mkdir(parents=True, exist_ok=True)
    (d / f'{tag}.txt').write_text('\n'.join([f'rule {name}', f'# {a.bin} ctau-window study 2026-10-31: {a.base} structure, ctau window [{lo}, {hi}] mm, changes {kv}', f'ct_low {lo}', f'ct_high {hi}'] + body) + '\n')
    tasks.append((name, kv))
print(f'{len(tasks)} window chains, {a.jobs} parallel', flush=True)


def one(t):
    name, kv = t
    log = BASE / 'logs' / f'{name}_{a.bin}.log'
    with open(log, 'w') as fh:
        rc = subprocess.run(['python3', str(BASE / 'study/prefit_variant.py'), name, a.bin, *kv, f'--rule={name}'], cwd=BASE, stdout=fh, stderr=subprocess.STDOUT).returncode
    res = BASE / 'study/prefit_variants' / f'{name}_{a.bin}' / 'roots' / y / 'fit2d' / f'fit2d_result_{tag}.root'
    if rc == 0 and res.exists():
        sys.path.insert(0, str(BASE / 'study')); import run_rules as rr
        if not rr.NOM: rr.nominal()
        r = rr.result_at(res, a.bin); print(rr.fmt(dict(r, rule=name)) if r else f'{name} no result', flush=True)
    else:
        print(f'{name:<34} FAILED rc={rc} result={res.exists()} ({log.relative_to(BASE)})', flush=True)


with ThreadPoolExecutor(a.jobs) as pool:
    list(pool.map(one, tasks))
print('done; python3 study/prefit_variant.py --summary')
