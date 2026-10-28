#!/usr/bin/env python3
"""Final 2D fit of knob sets on chosen bins, on the production prefits of ../25_OO_sept_acc_general_rule (roots/),
scored against nominal ver2 (scripts/nominal_v2.py).

  python3 scripts/run_v2.py --bins fwd_4_5,mid_7_8 RULE [RULE ...]   RULE = name in 25/study/rules.py RULES,
                                                                     a knob dir in rules/ (this folder) or in 25/study/rules/
  python3 scripts/run_v2.py --summary [PATTERN ...]
Results: runs/<rule>/roots/<yTag>/fit2d/, logs in runs/<rule>/logs/.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse, fnmatch, json, os, subprocess, sys
HERE = Path(__file__).resolve().parents[1]
P25 = HERE.parent / '25_OO_sept_acc_general_rule'
sys.path[:0] = [str(P25 / 'study'), str(HERE / 'scripts')]
import rules, run_rules as rr, nominal_v2  # noqa: E402
rr.NOM.update(nominal_v2.NOM)
ROOT = rr.ROOT


def knob_dir(name):
    """Own rules/ first, then 25/study/rules (FINAL_S3 etc.), else generate from 25 rules.py into own rules/."""
    for d in (HERE / 'rules' / name, P25 / 'study/rules' / name):
        if d.is_dir(): return d
    if name not in rules.RULES: sys.exit(f'unknown rule {name}')
    d = HERE / 'rules' / name; d.mkdir(parents=True)
    for r, e in rules.EDGES.items():
        for lo, hi in zip(e, e[1:]):
            y = 'y1.60_2.40' if r == 'fwd' else 'y0.00_1.60'
            c = rules.choice(name, r, lo); b = f'{r}_{lo}_{hi}'
            c.update(bin=b, purity=rules.PURITY[b], gaps=rules.GAPS[b])
            lines = [f'rule {name}', '# ' + ' '.join(f'{a}={v}' for a, v in c.items())] + [f'{a} {v}' for a, v in rules.knobs(c).items()]
            lines += [f'ct_low {rules.window(r, lo)[0]}', f'ct_high {rules.window(r, lo)[1]}']
            (d / f'{y}_pt{lo:.2f}_{hi:.2f}.txt').write_text('\n'.join(lines) + '\n')
    return d


def result(rule, b):
    y, t = rr.ytag(b)
    r = rr.result_at(HERE / 'runs' / rule / 'roots' / y / 'fit2d' / f'fit2d_result_{t}.root', b)
    return dict(r, rule=rule) if r else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rules', nargs='*'); ap.add_argument('--bins', default='')
    ap.add_argument('--jobs', type=int, default=32); ap.add_argument('--summary', action='store_true'); ap.add_argument('--force', action='store_true')
    ap.add_argument('--env', action='append', default=[], help='KEY=VALUE passed to fit2d (e.g. OO_MASS_LO=2.7); use with --tag')
    ap.add_argument('--tag', default='', help='suffix of the runs/ directory name')
    ap.add_argument('--lib', default='', help='fit2d_C.so to use (default: 25/build); build/ here has the macro/fit2d.C study hooks')
    a = ap.parse_args()
    if a.summary:
        names = sorted(d.name for d in (HERE / 'runs').iterdir() if d.is_dir() and any(fnmatch.fnmatch(d.name, p) for p in a.rules or ['*']))
        rows = [r for n in names for b in rr.NOM if (r := result(n, b))]
        for r in sorted(rows, key=lambda r: (r['bin'], abs(r['pct']))): print(rr.fmt(r))
        return
    own = [d.name for d in (HERE / 'rules').iterdir() if d.is_dir()]
    names = sorted({n for p in a.rules for n in (fnmatch.filter(own, p) + fnmatch.filter(rules.RULES, p) or [p])})
    bins = [b for b in a.bins.split(',') if b]
    dirs = {n: knob_dir(n) for n in names}
    lib = Path(a.lib).resolve() if a.lib else next((P25 / 'build').rglob('fit2d_C.so'))
    tasks = [(n, b) for n in names for b in bins if a.force or not result(n + a.tag, b)]
    print(f'{len(tasks)} fits, {a.jobs} parallel', flush=True)

    def one(t):
        n, b = t
        reg, lo, hi = b.split('_'); ylo, yhi = ('1.6', '2.4') if reg == 'fwd' else ('0', '1.6')
        out = HERE / 'runs' / (n + a.tag); (out / 'logs').mkdir(parents=True, exist_ok=True)
        env = {k: v for k, v in os.environ.items() if not k.startswith('OO_')}
        env.update(ROOT_HIST='0', OMP_NUM_THREADS='1', OO_RULE_DIR=str(dirs[n]), **dict(e.split('=', 1) for e in a.env))
        call = f'fit2d({lo},{hi},{ylo},{yhi},false,false,false,\\"{out}\\")'
        expr = f'gSystem->Load("{lib}");Int_t e=0;gROOT->ProcessLine("{call};",&e);gSystem->Exit(e?11:0);'
        with open(out / 'logs' / f'{b}.log', 'w') as log:
            rc = subprocess.run(['bash', '-c', f'ulimit -f 2000000; ulimit -v 4000000; exec {ROOT} -n -l -b -q -e {json.dumps(expr)}'],
                                cwd=P25, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
        r = result(n + a.tag, b) if rc == 0 else None
        print(rr.fmt(r) if r else f'{n:<34} {b:<10} FAILED rc={rc}', flush=True)

    with ThreadPoolExecutor(a.jobs) as pool:
        list(pool.map(one, tasks))


if __name__ == '__main__':
    main()
