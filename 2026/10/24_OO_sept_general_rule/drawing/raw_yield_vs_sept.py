#!/usr/bin/env python3
"""Raw yield, new acceptance vs September OO folder: drawing/raw_yield_vs_sept/raw_yield_vs_sept.{pdf,png}
(from raw_yield.txt written by dump_raw_yield_vs_sept.C)."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D = Path(__file__).resolve().parent / 'raw_yield_vs_sept'
rows = list(csv.DictReader(open(D / 'raw_yield.txt'), delimiter='\t'))
OLD, NEW, CAND, INK = '#2a78d6', '#eb6834', '#1baf7a', '#0b0b0b'
fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex='col', gridspec_kw=dict(height_ratios=[2, 1]))
for col, (reg, label) in enumerate([('fwd', '1.6 < |y| < 2.4'), ('mid', '|y| < 1.6')]):
    rs = [r for r in rows if r['region'] == reg]
    f = lambda k: [float(r[k]) for r in rs]
    x = [(a + b) / 2 for a, b in zip(f('ptLow'), f('ptHigh'))]; xe = [(b - a) / 2 for a, b in zip(f('ptLow'), f('ptHigh'))]
    no, noe, nn, nne = f('nsigOld'), f('nsigOldErr'), f('nsigNew'), f('nsigNewErr')
    co, cn = f('candOld'), f('candNew')
    ax = axes[0][col]
    ax.errorbar(x, no, xerr=xe, yerr=noe, fmt='s', ms=5, color=OLD, label='Sept (old acc.), N$_{sig}$')
    ax.errorbar(x, nn, xerr=xe, yerr=nne, fmt='o', ms=5, color=NEW, label='FINAL_v3 (new acc.), N$_{sig}$')
    ax.set_yscale('log'); ax.set_ylabel('fitted J/ψ yield N$_{sig}$ (raw)')
    ax.set_title(f'OO 5.36 TeV, {label}'); ax.legend(frameon=False)
    ax.grid(axis='y', color='0.9', lw=0.6)
    ax = axes[1][col]
    rsig = [b / a for a, b in zip(no, nn)]
    rsig_e = [r * ((ea / a) ** 2 + (eb / b) ** 2) ** 0.5 for r, a, b, ea, eb in zip(rsig, no, nn, noe, nne)]
    rc = [b / a for a, b in zip(co, cn)]
    ax.errorbar(x, rsig, xerr=xe, yerr=rsig_e, fmt='o', ms=5, color=INK, label='N$_{sig}$ new / Sept')
    ax.errorbar(x, rc, xerr=xe, fmt='D', ms=4, mfc='none', color=CAND, label='dimuon candidates new / Sept')
    for xi, r in zip(x, rsig):
        ax.annotate(f'{r:.2f}', (xi, r), textcoords='offset points', xytext=(0, 7), ha='center', fontsize=7, color='0.3')
    ax.axhline(1, color='0.3', lw=0.8); ax.set_ylim(0, 1.15); ax.set_ylabel('new / Sept')
    ax.set_xlabel('$p_T$ [GeV/c]'); ax.grid(axis='y', color='0.9', lw=0.6)
    if col == 0: ax.legend(frameon=False, fontsize=8, loc='lower right')
    if reg == 'mid':
        for a in (ax, axes[0][col]):
            a.set_xscale('log'); a.set_xticks([7, 8, 10, 12, 14, 16, 20, 30]); a.set_xticklabels(['7', '8', '10', '12', '14', '16', '20', '30']); a.minorticks_off()
        yt = [600, 1000, 2000, 3000, 4000, 6000]; axes[0][col].set_yticks(yt); axes[0][col].set_yticklabels([str(v) for v in yt])
fig.tight_layout()
for ext in ('pdf', 'png'): fig.savefig(D / f'raw_yield_vs_sept.{ext}', dpi=130)
print('wrote', D / 'raw_yield_vs_sept.pdf')
