#!/usr/bin/env python3
"""fwd 10-12 (or any bin): final fit with B fixed (fit2d hook OO_B_FIX) on the production prefits, knob dir of a rule.
  python3 scripts/bfix_scan.py BIN RULE B1 B2 ...     -> runs/BFIX_<RULE>_<B>/  (one fit per B, parallel)
  python3 scripts/bfix_scan.py --compare FILE_FREE FILE_FIXED   -> parameter shifts (fixed - free) / free error, and minNll difference
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json, os, subprocess, sys
HERE = Path(__file__).resolve().parents[1]
P25 = HERE.parent / '25_OO_sept_acc_general_rule'
sys.path[:0] = [str(P25 / 'study'), str(HERE / 'scripts')]
import run_rules as rr  # noqa: E402
from run_v2 import knob_dir  # noqa: E402
ROOT = rr.ROOT

DUMP = r'''TFile f("%s"); auto*r=f.Get<RooFitResult>("modelResult"); if(!r)throw std::runtime_error("Missing modelResult");
auto*raw=f.Get<TParameter<double>>("rawFinalNll"); if(!raw)throw std::runtime_error("Missing rawFinalNll; offset minNll values cannot be compared across fits");
printf("NLL %%.17g st %%d cq %%d\n",raw->GetVal(),r->status(),r->covQual());
for(auto*a:r->floatParsFinal()){auto*v=(RooRealVar*)a; printf("P %%s %%.8g %%.8g %%.8g %%.8g\n",v->GetName(),v->getVal(),v->getError(),v->getMin(),v->getMax());}
for(auto*a:r->constPars()){auto*v=dynamic_cast<RooRealVar*>(a); if(v) printf("C %%s %%.8g\n",v->GetName(),v->getVal());}'''


def dump(f):
    out = subprocess.run([ROOT, '-n', '-l', '-b', '-q', '-e', DUMP % f], capture_output=True, text=True).stdout
    d = {'P': {}, 'C': {}}
    for l in out.splitlines():
        p = l.split()
        if p[:1] == ['NLL']: d['nll'], d['st'], d['cq'] = float(p[1]), int(p[3]), int(p[5])
        elif p[:1] == ['P']: d['P'][p[1]] = tuple(map(float, p[2:]))
        elif p[:1] == ['C']: d['C'][p[1]] = float(p[2])
    return d


def compare(free, fixed):
    a, b = dump(free), dump(fixed)
    if 'nll' not in a or 'nll' not in b:
        raise ValueError('Both inputs must contain modelResult and rawFinalNll; offset minNll is not a fallback.')
    print(f"rawFinalNll free {a['nll']:.6f} (st{a['st']} cq{a['cq']})  fixed {b['nll']:.6f} (st{b['st']} cq{b['cq']})  2dNLL = {2 * (b['nll'] - a['nll']):.6f}")
    B = a['P'].get('bFraction'); print(f"B free {B[0]:.5f} ± {B[1]:.5f}  fixed {b['C'].get('bFraction', float('nan')):.5f}")
    rows = []
    for n, (v, e, lo, hi) in a['P'].items():
        if n in b['P']:
            v2, e2, lo2, hi2 = b['P'][n]
            lim = ' AT-LIMIT' if min(v2 - lo2, hi2 - v2) < 0.01 * (hi2 - lo2) else ''
            rows.append(((v2 - v) / e if e > 0 else 0, n, v, e, v2, lim))
    for s, n, v, e, v2, lim in sorted(rows, key=lambda r: -abs(r[0])):
        print(f'  {n:<34} free {v:12.6g} ± {e:<10.3g} fixed {v2:12.6g}  shift {s:+6.2f} err{lim}')


def main():
    if sys.argv[1] == '--compare': return compare(*sys.argv[2:4])
    b, rule, Bs = sys.argv[1], sys.argv[2], sys.argv[3:]
    kd = knob_dir(rule)
    lib = HERE / 'build_bfix/fit2d_C.so'  # macro/fit2d.C copy: OO_B_FIX keeps the B domain
    reg, lo, hi = b.split('_'); ylo, yhi = ('1.6', '2.4') if reg == 'fwd' else ('0', '1.6')

    def one(B):
        out = HERE / 'runs' / f'BFIX_{rule}_{B}'; (out / 'logs').mkdir(parents=True, exist_ok=True)
        env = {k: v for k, v in os.environ.items() if not k.startswith('OO_')}
        env.update(ROOT_HIST='0', OMP_NUM_THREADS='1', OO_RULE_DIR=str(kd), OO_B_FIX=B)
        call = f'fit2d({lo},{hi},{ylo},{yhi},false,false,false,\\"{out}\\")'
        expr = f'gSystem->Load("{lib}");Int_t e=0;gROOT->ProcessLine("{call};",&e);gSystem->Exit(e?11:0);'
        with open(out / 'logs' / f'{b}.log', 'w') as log:
            rc = subprocess.run(['bash', '-c', f'ulimit -f 2000000; exec {ROOT} -n -l -b -q -e {json.dumps(expr)}'], cwd=P25, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
        print(B, 'rc', rc, flush=True)

    with ThreadPoolExecutor(len(Bs)) as pool: list(pool.map(one, Bs))


if __name__ == '__main__':
    main()
