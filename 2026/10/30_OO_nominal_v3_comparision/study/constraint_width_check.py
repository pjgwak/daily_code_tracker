#!/usr/bin/env python3
"""Is a joint strength k a 'large' number? (--bins=fwd_10_12,fwd_12_14,fwd_14_20 for the forward joint x16 bins) Compare, per mid bin and PR parameter: the MC prefit error, the joint-constraint
width sqrt(k) x MC error used in the final fit, the parameter value, and the data-only error from a fit with the PR widths free
(study/runs/V3_free_free_fix_lf_free). Usage: python3 study/constraint_width_check.py [k ...]"""
import csv, subprocess, sys
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]; ROOT = str(Path.home() / 'build/root-6.36.06/bin/root')
args = [a for a in sys.argv[1:] if not a.startswith('--bins=')]
ks = [float(x) for x in args] or [64, 384, 1024]
sel = next((a.split('=', 1)[1].split(',') for a in sys.argv[1:] if a.startswith('--bins=')), None)  # --bins=fwd_10_12,... (default: mid)
err = {}
for r in csv.DictReader(open(BASE / 'config/pr_width_err.txt'), delimiter='\t'):
    err.setdefault(r['bin'], {})[r['parameter']] = (float(r['value']), float(r['error']))
bins = sel or ['mid_7_8', 'mid_8_9', 'mid_9_10', 'mid_10_12', 'mid_12_14', 'mid_14_16', 'mid_16_20', 'mid_20_30']
pars = ['ctauTime1Mean', 'ctauTime1Scale', 'ctauTime2Delta', 'ctauTime3Delta', 'ctauTime4Delta']
def free_errors(b):
    reg, lo, hi = b.split('_'); y = 'y1.60_2.40' if reg == 'fwd' else 'y0.00_1.60'
    f = BASE / 'study/runs/V3_free_free_fix_lf_free/roots' / y / 'fit2d' / f'fit2d_result_{y}_pt{float(lo):.2f}_{float(hi):.2f}.root'
    if not f.exists(): return {}
    expr = f'TFile f("{f}"); auto*r=f.Get<RooFitResult>("modelResult"); for(auto*a:r->floatParsFinal()){{auto*v=(RooRealVar*)a; printf("PAR %s %.6g %.6g\\n",v->GetName(),v->getVal(),v->getError());}}'
    out = subprocess.run([ROOT, '-n', '-l', '-b', '-q', '-e', expr], capture_output=True, text=True).stdout
    return {l.split()[1]: (float(l.split()[2]), float(l.split()[3])) for l in out.splitlines() if l.startswith('PAR')}
print('bin        parameter       MC value   MC err   ' + ' '.join(f'sqrt({k:g})xerr' for k in ks) + '   data-only err (widths free)  data/MC-err')
for b in bins:
    fe = free_errors(b)
    for p in pars:
        if p not in err.get(b, {}): continue
        v, e = err[b][p]
        d = fe.get(p)
        print(f'{b:<10} {p:<15} {v:9.5f} {e:9.5f}   ' + ' '.join(f'{(k**0.5)*e:12.5f}' for k in ks) + (f'   {d[1]:9.5f} (val {d[0]:.5f})  {d[1]/e:6.1f}' if d else '   (free fit not available)'))
