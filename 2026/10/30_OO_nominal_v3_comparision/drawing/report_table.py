#!/usr/bin/env python3
"""20-bin result table: drawing/final_summary/report_table.txt (from summary.txt and roots/*/fit2d).

Columns as in the 24 folder. pr_out / np_out are the signal PR / NP probabilities outside the fitted ctau
window, computed from the final-fit parameters (diagnostics_*/FINAL_MODEL_PARAMETERS.csv):
  PR: Gaussian mixture (widths sigma1 + gaps, recursive fractions), both edges;
  NP: exponential mixture beyond ct_high (weights from the native NP yields; resolution smearing ignored).
"""
import csv, math, subprocess
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
ROOT = str(Path.home() / 'build/root-6.36.06/bin/root')
D = BASE / 'drawing' / 'final_summary'


def phi(x): return 0.5 * math.erfc(-x / math.sqrt(2))


def extra(tag, y):
    f = BASE / 'roots' / y / 'fit2d' / f'fit2d_result_{tag}.root'
    out = subprocess.run([ROOT, '-l', '-b', '-q', '-e', f'TFile f("{f}"); for(auto n:{{"appliedCtauLow","appliedCtauHigh","npFirstLifetimeFinal","npSecondLifetimeFinal"}}){{auto*p=f.Get<TParameter<double>>(n); printf("X %s %.10g\\n",n,p?p->GetVal():NAN);}}'],
                         capture_output=True, text=True).stdout
    return {l.split()[1]: float(l.split()[2]) for l in out.splitlines() if l.startswith('X ')}


def outside(tag, y, x):
    p = {r['name']: float(r['value']) for r in csv.DictReader(open(BASE / 'roots' / y / 'fit2d' / f'diagnostics_{tag}' / 'FINAL_MODEL_PARAMETERS.csv'))}
    lo, hi = x['appliedCtauLow'], x['appliedCtauHigh']
    widths, s = [], 0.0
    for i, k in enumerate(['ctauTime1Scale', 'ctauTime2Delta', 'ctauTime3Delta', 'ctauTime4Delta']):
        if k in p: s += p[k]; widths.append(s)
    fr, rem = [], 1.0
    for i in range(len(widths) - 1):
        f = p.get(f'ctauRecursiveFrac{i + 1}', 0.0); fr.append(rem * f); rem *= 1 - f
    fr.append(rem)
    mean = p.get('ctauTime1Mean', 0.0)
    pr = sum(w * (phi((lo - mean) / sg) + 1 - phi((hi - mean) / sg)) for w, sg in zip(fr, widths))
    taus = [t for t in (x.get('npFirstLifetimeFinal'), x.get('npSecondLifetimeFinal')) if t and math.isfinite(t) and t > 0]
    ys = [p.get(f'NsignalSS{i + 1}', 0.0) for i in range(len(taus))]
    tot = sum(ys) or 1.0
    np_ = sum(yv / tot * math.exp(-hi / t) for yv, t in zip(ys, taus)) if len(ys) == len(taus) and sum(ys) > 0 else (math.exp(-hi / taus[0]) if taus else float('nan'))
    return lo, hi, pr, np_


rows = list(csv.DictReader(open(D / 'summary.txt'), delimiter='\t'))
cols = ['bin', 'B', 'Berr', 'Bnom', 'percent', 'pull', 'status', 'covQual', 'massChi2', 'massNdf', 'massP', 'timeChi2', 'timeNdf', 'timeP', 'ct_low', 'ct_high', 'pr_out', 'np_out']
with open(D / 'report_table.txt', 'w') as fh:
    fh.write('\t'.join(cols) + '\n')
    for r in rows:
        y = 'y1.60_2.40' if r['region'] == 'fwd' else 'y0.00_1.60'
        tag = f'{y}_pt{float(r["ptLow"]):.2f}_{float(r["ptHigh"]):.2f}'
        lo, hi, pr, np_ = outside(tag, y, extra(tag, y))
        vals = [f'{r["region"]}_{r["ptLow"]}_{r["ptHigh"]}', r['B'], r['Berr'], r['Bnom'], r['percent'], r['pull'], r['status'], r['covQual'],
                r['massChi2'], r['massNdf'], r['massPvalue'], r['timeChi2'], r['timeNdf'], r['timePvalue'], lo, hi, pr, np_]
        fh.write('\t'.join(str(v) for v in vals) + '\n')
print('wrote', D / 'report_table.txt')
