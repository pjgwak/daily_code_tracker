#!/usr/bin/env python3
"""B fraction vs nominal for FINAL_V4: drawing/final_summary/nominal_agreement.{pdf,png} (from summary.txt)."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D = Path(__file__).resolve().parent / 'final_summary'
rows = list(csv.DictReader(open(D / 'summary.txt'), delimiter='\t'))
fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex='col', gridspec_kw=dict(height_ratios=[2, 1]))
for col, (reg, label) in enumerate([('fwd', '1.6 < |y| < 2.4'), ('mid', '|y| < 1.6')]):
    rs = [r for r in rows if r['region'] == reg]
    x = [(float(r['ptLow']) + float(r['ptHigh'])) / 2 for r in rs]; xe = [(float(r['ptHigh']) - float(r['ptLow'])) / 2 for r in rs]
    B = [float(r['B']) for r in rs]; Be = [float(r['Berr']) for r in rs]
    N = [float(r['Bnom']) for r in rs]; Ne = [float(r['BnomErr']) for r in rs]
    ax = axes[0][col]
    ax.errorbar(x, N, xerr=xe, yerr=Ne, fmt='s', color='0.45', mfc='none', label='nominal')
    ax.errorbar(x, B, xerr=xe, yerr=Be, fmt='o', color='#2a4ea0', label='this fit (FINAL_V4)')
    ax.set_title(f'OO 5.36 TeV, {label}'); ax.set_ylabel('nonprompt fraction B'); ax.legend(frameon=False)
    ax = axes[1][col]
    ax.axhspan(-2, 2, color='#e9f1ec')
    ax.errorbar(x, [100 * (b / n - 1) for b, n in zip(B, N)], xerr=xe, yerr=[100 * e / n for e, n in zip(Be, N)], fmt='o', color='#2a4ea0')
    ax.axhline(0, color='0.3', lw=0.8); ax.set_ylim(-12, 12); ax.set_ylabel('B/B$_{nom}$ − 1 [%]'); ax.set_xlabel('$p_T$ [GeV/c]')
    if reg == 'mid': ax.set_xscale('log'); axes[0][col].set_xscale('log')
fig.tight_layout()
for ext in ('pdf', 'png'): fig.savefig(D / f'nominal_agreement.{ext}', dpi=130)
print('wrote', D / 'nominal_agreement.pdf')
