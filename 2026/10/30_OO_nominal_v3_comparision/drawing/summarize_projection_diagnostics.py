#!/usr/bin/env python3
"""Collect bin-integrated projection checks; never choose binning by best chi2."""
import argparse, csv, json, math, re
from pathlib import Path
import subprocess

p=argparse.ArgumentParser()
p.add_argument('directory',nargs='?',default=str(Path(__file__).resolve().parent/'projection_diagnostics'))
a=p.parse_args(); directory=Path(a.directory)
rows=[]
pattern=re.compile(r'(.*?)_y([\d.]+)_([\d.]+)_pt([\d.]+)_([\d.]+)(?:_bkgOff)?$')
for f in sorted(directory.glob('*.csv')):
    match=pattern.fullmatch(f.stem)
    if not match: continue
    stage,yl,yh,pl,ph=match.groups()
    for r in csv.DictReader(f.open()):
        r.update(stage=stage,yLow=yl,yHigh=yh,ptLow=pl,ptHigh=ph,pdf=f.with_suffix('.pdf').name)
        ndf=int(r['grouped_ndf'])
        r['grouped_pearson_ndf']=float(r['grouped_pearson'])/ndf if ndf>0 else math.nan
        r['grouped_deviance_ndf']=float(r['grouped_deviance'])/ndf if ndf>0 else math.nan
        rows.append(r)
if not rows: raise SystemExit('No projection diagnostics found')
with (directory/'summary.csv').open('w') as out:
    writer=csv.DictWriter(out,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
selected=[r for r in rows if r['view']==('signal_display' if 'ctau' in r['stage'] else 'display_full')]
selected.sort(key=lambda r:(r['stage'],float(r['yLow']),float(r['ptLow'])))
lines=['# Projection diagnostics','',
       'Saved fits are replayed without minimization. Full-fit normalization is preserved in every window.',
       'Signal windows contain the union of the central 99% intervals of PR and NP separately. Core windows contain the central 99% of PR; they need not contain 99% of NP.',
       'Prefit lifetime windows use the central 99% of their own PDF. Mass uses the full fit window.',
       '',
       'Display bin count is max(400, 2*sqrt(N_eff)), rounded up to a multiple of 50; full displays also retain at least the native bin count. The rule depends on statistics, never on chi-square. N_eff=(sum weights)^2/sum weights^2.',
       'Dense signal views use twice the display bins. Native/full 60/120/240 and signal 60/120/240/480 scans remain in `summary.csv`, so fine binning cannot hide a broad residual. The table reports the dense display and the fixed 120-bin comparison together.',
       'Expected bin contents are exact PDF integrals. max_center_bias_pull measures the largest bin-center approximation error in units of sqrt(expected); that approximation is not used for residuals.',
       'Roo chi-square excludes empty observed bins and uses asymmetric data errors. Pearson/deviance include empty bins. Grouping merges adjacent expected counts to at least 5, with the final remainder merged backwards.',
       'NDF = bins minus projection-relevant floating parameters is a descriptive approximation for constrained fits and restricted projections. No calibrated p-value or automatic fit acceptance is inferred.',
       '', '| Stage | y | pT | Window | Display bins | Roo chi2/ndf | Roo / 120 bins | Grouped Pearson/ndf | Center bias (sigma) | Sparse bins | PDF |',
       '|---|---|---|---|---:|---:|---:|---:|---:|---:|---|']
for r in selected:
    val=float(r['roo_chi2_ndf']);marked=f'**{val:.3f}**' if val>=2 else f'{val:.3f}'
    coarse=next(q for q in rows if q['pdf']==r['pdf'] and q['view']==('signal' if 'ctau' in r['stage'] else 'full') and q['bins']=='120')
    lines.append(f"| {r['stage']} | {r['yLow']}–{r['yHigh']} | {r['ptLow']}–{r['ptHigh']} | {float(r['low']):.3f}–{float(r['high']):.3f} | {r['bins']} | {marked} | {float(coarse['roo_chi2_ndf']):.3f} | {r['grouped_pearson_ndf']:.3f} | {float(r['max_center_bias_pull']):.3g} | {r['expected_lt5']}/{r['bins']} | [plot]({r['pdf']}) |")
(directory/'README.md').write_text('\n'.join(lines)+'\n')
# Fixed-window/fixed-binning overview: report full and focused results together.
final=[r for r in rows if r['stage']=='fit2d_ctau']
keys=sorted({(float(r['yLow']),float(r['ptLow']),float(r['ptHigh'])) for r in final})
if keys:
    specs=[('native',None),('full','120'),('signal','60'),('signal','120'),('signal','240'),('signal','480'),('core','120'),('signal_display',None),('signal_dense',None)]
    values=[[math.nan for _ in specs] for _ in keys]
    for i,key in enumerate(keys):
        for j,(view,bins) in enumerate(specs):
            rs=[r for r in final if (float(r['yLow']),float(r['ptLow']),float(r['ptHigh']))==key and r['view']==view and (bins is None or r['bins']==bins)]
            if rs: values[i][j]=float(rs[0]['roo_chi2_ndf'])
    code=['#include "TCanvas.h"','#include "TH2D.h"','#include "TStyle.h"','#include "TColor.h"','void make_binning_overview() {',
          'gStyle->SetOptStat(0); double stops[]={0.,.4,1.}, red[]={.70,.98,1.}, green[]={.83,.98,.53}, blue[]={1.,.98,.43}; TColor::CreateGradientColorTable(3,stops,red,green,blue,50); gStyle->SetPaintTextFormat(".2f");',
          'TCanvas c("c","Projection binning comparison",1200,1000); c.SetLeftMargin(.15); c.SetBottomMargin(.15); c.SetRightMargin(.14);',
          f'TH2D h("h","Lifetime projection: same fit, different binning; ; ; Roo chi-square / approximate NDF",{len(specs)},0,{len(specs)},{len(keys)},0,{len(keys)});',
          'h.SetMinimum(0.); h.SetMaximum(2.5); h.SetMarkerSize(1.0); h.GetXaxis()->SetLabelSize(.025); h.GetYaxis()->SetLabelSize(.027);']
    labels=['Native full','Full / 120','Signal / 60','Signal / 120','Signal / 240','Signal / 480','PR / 120','Signal / N','Signal / 2N']
    for j,label in enumerate(labels): code.append(f'h.GetXaxis()->SetBinLabel({j+1},{json.dumps(label)});')
    for i,key in enumerate(keys):
        label=f"{'Mid' if key[0]==0 else 'Fwd'} {key[1]:g}-{key[2]:g}"
        code.append(f'h.GetYaxis()->SetBinLabel({len(keys)-i},{json.dumps(label)});')
        for j,value in enumerate(values[i]):
            if math.isfinite(value): code.append(f'h.SetBinContent({j+1},{len(keys)-i},{value:.17g});')
    code+=['h.Draw("COLZ TEXT");',f'c.Print({json.dumps(str((directory/"binning_overview.pdf").resolve()))});','}']
    macro=directory/'make_binning_overview.C'
    macro.write_text('\n'.join(code)+'\n')
    subprocess.run(['root','-n','-l','-b','-q',str(macro)],check=True)

print(json.dumps(dict(files=len({r['pdf'] for r in rows}),scenarios=len(rows),table_rows=len(selected),display_ge2=[r['pdf'] for r in selected if float(r['roo_chi2_ndf'])>=2], coarse_signal120_ge2=[r['pdf'] for r in rows if r['view']=='signal' and r['bins']=='120' and float(r['roo_chi2_ndf'])>=2])))
