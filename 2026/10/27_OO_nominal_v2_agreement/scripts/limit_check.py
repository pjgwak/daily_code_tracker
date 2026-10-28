#!/usr/bin/env python3
"""List parameters within 1% of a bound (or within 0.5 error) in a work dir's background ctau prefit and final fit.
  python3 scripts/limit_check.py WORKDIR [RULE]"""
import subprocess, sys
from pathlib import Path
ROOT = str(Path.home() / 'build/root-6.36.06/bin/root')
EXPR = r'''TFile f("%s"); auto*r=f.Get<RooFitResult>("%s"); if(!r){printf("no %s\n");return;} printf("status %%d cq %%d\n",r->status(),r->covQual());
for(auto*a:r->floatParsFinal()){auto*v=(RooRealVar*)a; double lo=v->getMin(),hi=v->getMax(),x=v->getVal(),e=v->getError();
 if(std::isfinite(lo)&&std::isfinite(hi)&&(x-lo<0.01*(hi-lo)||hi-x<0.01*(hi-lo)||(e>0&&(x-lo<0.5*e||hi-x<0.5*e)))) printf("  LIMIT %%-34s %%.5g [%%.5g, %%.5g] err %%.3g\n",v->GetName(),x,lo,hi,e);
 if(e>fabs(x)&&TString(v->GetName()).BeginsWith("N")) printf("  ERR>VAL %%-32s %%.5g +- %%.3g\n",v->GetName(),x,e);}'''
def show(f, key):
    out = subprocess.run([ROOT, '-n', '-l', '-b', '-q', '-e', EXPR % (f, key, key)], capture_output=True, text=True).stdout
    print(f'{f}:'); print('\n'.join(l for l in out.splitlines() if l.strip() and not l.startswith(('Processing', ' ****'))))
w = Path(sys.argv[1])
for f in sorted(w.glob('roots/y*/ctau_bkg/ctau_bkg_fitresult_*.root')): show(f, 'timeResult')
rule = sys.argv[2] if len(sys.argv) > 2 else None
for f in sorted((w / 'rules' / rule if rule else w).glob('roots/y*/fit2d/fit2d_result_*.root')): show(f, 'modelResult')
