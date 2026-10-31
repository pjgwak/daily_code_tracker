#!/usr/bin/env python3
"""B vs nominal and the derived PR / NP yield systematics per bin, fwd and mid.

  python3 drawing/v3_systematics.py [--ymax=5] OUT_PREFIX SUMMARY[:label] [SUMMARY2[:label2] ...]

--ymax fixes the y axis of the systematic panels (wide view, e.g. 0-5 %); without it the axis is zoomed to the data.

SUMMARY is a workflow/collect.C table (B, Berr, Bnom, BnomErr per bin). Writes OUT_PREFIX_bfrac.png,
OUT_PREFIX_sys_np.png, OUT_PREFIX_sys_pr.png (and .pdf). NP systematic = |B-Bnom|/Bnom, PR = |B-Bnom|/(1-Bnom).
"""
import csv, sys
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
args = [a for a in sys.argv[1:] if not a.startswith('--ymax=')]
ymax = next((float(a.split('=', 1)[1]) for a in sys.argv[1:] if a.startswith('--ymax=')), None)
out = Path(args[0]); series = []
for a in args[1:]:
    f, _, lab = a.partition(':'); series.append((lab or Path(f).stem, list(csv.DictReader(open(f), delimiter='\t'))))
colors = ['#c0392b', '#2a4ea0', '#2e8b57', '#8e44ad']
REG = [('fwd', '1.6 < |y| < 2.4'), ('mid', '|y| < 1.6')]


def xs(rs): return [(float(r['ptLow']) + float(r['ptHigh'])) / 2 for r in rs], [(float(r['ptHigh']) - float(r['ptLow'])) / 2 for r in rs]


fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
for ax, (reg, label) in zip(axes, REG):
    rs0 = [r for r in series[0][1] if r['region'] == reg]; x, xe = xs(rs0)
    ax.errorbar(x, [float(r['Bnom']) for r in rs0], xerr=xe, yerr=[float(r['BnomErr']) for r in rs0], fmt='o', color='k', label='Nominal (ver3)')
    for (lab, rows), c in zip(series, colors):
        rs = [r for r in rows if r['region'] == reg]; x, xe = xs(rs)
        ax.errorbar([v + 0.08 * i for i, v in enumerate(x)] if False else x, [float(r['B']) for r in rs], xerr=xe, yerr=[float(r['Berr']) for r in rs], fmt='s', ms=4, color=c, ls=':', label=f'2D fit ({lab})')
    ax.set_title(f'OO 5.36 TeV, {label}'); ax.set_xlabel('$p_T$ [GeV/c]'); ax.set_ylabel('Non-prompt fraction'); ax.legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig(f'{out}_bfrac.png', dpi=130); fig.savefig(f'{out}_bfrac.pdf')

for kind, title, fn in (('np', 'Non-prompt yield', lambda b, n: abs(b - n) / n), ('pr', 'Prompt yield', lambda b, n: abs(b - n) / (1 - n))):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    for ax, (reg, label) in zip(axes, REG):
        for (lab, rows), c in zip(series, colors):
            rs = [r for r in rows if r['region'] == reg]
            for i, r in enumerate(rs):
                lo, hi = float(r['ptLow']), float(r['ptHigh']); v = 100 * fn(float(r['B']), float(r['Bnom']))
                ax.hlines(v, lo, hi, color=c, lw=2.2, label=f'2D fit ({lab})' if i == 0 else None)
        ax.axhspan(0, 2, color='#e9f1ec', zorder=0); ax.set_ylim(0, ymax); ax.axhline(2, color='#2f7d5b', lw=0.8, ls='--'); ax.set_xlabel('$p_T$ [GeV/c]'); ax.set_ylabel('Relative systematic uncertainty [%]')
        ax.set_title(f'{title}, {label}'); ax.legend(frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(f'{out}_sys_{kind}.png', dpi=130); fig.savefig(f'{out}_sys_{kind}.pdf')
print('wrote', f'{out}_bfrac.png', f'{out}_sys_np.png', f'{out}_sys_pr.png')
