#!/usr/bin/env python3
"""Run the final 2D fit for rule knob sets on chosen bins, using the prefits in roots/.

  python3 study/run_rules.py --bins mid_7_8,fwd_3_4 RULE_PATTERN [...]   (fnmatch patterns over study/rules.py RULES)
  python3 study/run_rules.py --summary [RULE_PATTERN ...]                 (table of finished runs only)

Results: study/runs/<rule>/roots/<yTag>/fit2d/ (one directory per rule); log next to it.
Each finished fit prints one line immediately (B, % and sigma vs nominal, status).
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse, csv, fnmatch, json, os, subprocess, sys
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / 'study'))
import rules  # noqa: E402
ROOT = str(Path.home() / 'build/root-6.36.06/bin/root')
NOM = {}


def nominal():
    """Nominal B per bin from utils/nominal_bfraction.h (forward bins() then midBins())."""
    import re
    text = (BASE / 'utils/nominal_bfraction.h').read_text()
    fwd, mid = text.split('midBins()', 1)
    for reg, part in (('fwd', fwd), ('mid', mid)):
        for lo, hi, v, e in re.findall(r'\{\s*([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\d.]+)', part):
            NOM[f'{reg}_{lo}_{hi}'] = (float(v), float(e))


def ytag(b):
    r, lo, hi = b.split('_')
    y = 'y1.60_2.40' if r == 'fwd' else 'y0.00_1.60'
    return y, f'{y}_pt{float(lo):.2f}_{float(hi):.2f}'


def done(rule, b):
    y, t = ytag(b)
    return (BASE / 'study/runs' / rule / 'roots' / y / 'fit2d' / f'fit2d_result_{t}.root').exists()


def result(rule, b):
    y, t = ytag(b)
    r = result_at(BASE / 'study/runs' / rule / 'roots' / y / 'fit2d' / f'fit2d_result_{t}.root', b)
    return dict(r, rule=rule) if r else None


def result_at(f, b):
    if not f.exists(): return None
    out = subprocess.run([ROOT, '-n', '-l', '-b', '-q', '-e', f'TFile f("{f}"); auto*r=f.Get<RooFitResult>("modelResult"); auto*b=r?(RooRealVar*)r->floatParsFinal().find("bFraction"):nullptr;'
                          'auto P=[&](const char*n){auto*p=f.Get<TParameter<double>>(n);return p?p->GetVal():NAN;};'
                          'if(b) printf("RES %.8g %.8g %d %d %.4g %.4g\\n",b->getVal(),b->getError(),r->status(),r->covQual(),P("massPvalue"),P("timePvalue"));'],
                         capture_output=True, text=True).stdout
    l = next((l.split()[1:] for l in out.splitlines() if l.startswith('RES')), None)
    if not l: return None
    B, Be, st, cq, mp, tp = float(l[0]), float(l[1]), int(l[2]), int(l[3]), float(l[4]), float(l[5])
    v, e = NOM[b]
    pct, pull = 100 * (B / v - 1), (B - v) / e
    ok = st == 0 and cq == 3 and (abs(pct) <= 2 or abs(pull) <= 0.138)
    return dict(rule='', bin=b, B=B, Berr=Be, pct=pct, pull=pull, status=st, covQual=cq, massP=mp, timeP=tp, pass_=ok)


def fmt(r):
    return (f'{r["rule"]:<34} {r["bin"]:<10} B={r["B"]:.4f}±{r["Berr"]:.4f} {r["pct"]:+6.2f}% {r["pull"]:+5.2f}σ '
            f'st{r["status"]} cq{r["covQual"]} mp {r["massP"]:.2f} tp {r["timeP"]:.2f} {"PASS" if r["pass_"] else "fail"}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('patterns', nargs='*'); ap.add_argument('--bins', default='')
    ap.add_argument('--jobs', type=int, default=40); ap.add_argument('--summary', action='store_true')
    ap.add_argument('--force', action='store_true', help='refit bins that already have a result')
    a = ap.parse_args()
    nominal()
    names = sorted({n for p in a.patterns or ['*'] for n in fnmatch.filter(rules.RULES, p)})
    # hand-written knob directories (study/rules/<name>/, not defined in rules.py) are run as they are
    custom = sorted({d.name for p in a.patterns for d in (BASE / 'study/rules').glob(p) if d.is_dir() and d.name not in rules.RULES})
    if a.summary:
        runs = BASE / 'study/runs'
        names = sorted(d.name for d in runs.iterdir() if d.is_dir() and any(fnmatch.fnmatch(d.name, p) for p in a.patterns or ['*']))
        rows = [r for n in names for b in NOM if (r := result(n, b))]
        for r in sorted(rows, key=lambda r: (r['bin'], abs(r['pct']))): print(fmt(r))
        return
    bins = [b for b in a.bins.split(',') if b]
    if not (names or custom) or not bins: sys.exit('need rule patterns and --bins')
    for n in names: rules.write_final_s1(n) if n.startswith('FINAL_S') else rules.write(n)  # FINAL_S* carry the hand-written fwd 2-3 file
    names += custom
    build = BASE / 'build'
    lib = next(build.rglob('fit2d_C.so'))
    tasks = [(n, b) for n in names for b in bins if a.force or not done(n, b)]
    print(f'{len(tasks)} fits ({len(names)} rules x {len(bins)} bins), {a.jobs} parallel', flush=True)

    def one(t):
        n, b = t
        reg, lo, hi = b.split('_'); ylo, yhi = ('1.6', '2.4') if reg == 'fwd' else ('0', '1.6')
        out = BASE / 'study/runs' / n; (out / 'logs').mkdir(parents=True, exist_ok=True)
        env = {k: v for k, v in os.environ.items() if not k.startswith('OO_')}
        env.update(ROOT_HIST='0', OMP_NUM_THREADS='1', OO_RULE_DIR=str(BASE / 'study/rules' / n))
        call = f'fit2d({lo},{hi},{ylo},{yhi},false,false,false,\\"{out}\\")'
        expr = f'gSystem->Load("{lib}");Int_t e=0;gROOT->ProcessLine("{call};",&e);gSystem->Exit(e?11:0);'
        with open(out / 'logs' / f'{b}.log', 'w') as log:
            rc = subprocess.run(['bash', '-c', f'ulimit -f 2000000; exec {ROOT} -n -l -b -q -e {json.dumps(expr)}'], cwd=BASE, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
        r = result(n, b) if rc == 0 else None
        print(fmt(r) if r else f'{n:<34} {b:<10} FAILED rc={rc}', flush=True)

    with ThreadPoolExecutor(a.jobs) as pool:
        list(pool.map(one, tasks))


if __name__ == '__main__':
    main()
