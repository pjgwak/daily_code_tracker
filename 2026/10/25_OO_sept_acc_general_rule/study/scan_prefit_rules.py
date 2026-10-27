#!/usr/bin/env python3
"""Random scan of prefit choices x MC-centred final-fit rules for one bin (2026-10-27, forward 10-14 study).

  python3 study/scan_prefit_rules.py NAME BIN --n 150 --seed 1 --jobs 32 --rules R1,R2,...
  python3 study/scan_prefit_rules.py --summary NAME BIN          (results.tsv sorted by |pct|)

Each combo is a prefit_variant-style work directory study/scans/<NAME>_<BIN>/c<k>/ (prefits refitted from the
first changed stage, then prefit_audit.txt), on which every rule's fit2d is run (rules/<rule>/roots/.../fit2d).
Combo 0 is the current config/prefit_settings.txt row. Results: study/scans/<NAME>_<BIN>/results.tsv.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse, csv, json, os, random, shutil, subprocess, sys
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / 'study'))
import rules, run_rules as rr  # noqa: E402
from prefit_variant import FIRST, ORDER, LINKS, stages_from  # noqa: E402
ROOT = rr.ROOT

AXES = {  # (value, weight); '' = keep baseline
    'pr_nres': [('4', 3), ('3', 1)],
    'pr_quant': [('0.0005,0.9995', 2), ('0.0001,0.9999', 1), ('0.002,0.998', 1)],
    'np_count': [('1', 1), ('2', 1)],
    'np_range': [('', 2), ('0,3', 1), ('0,10', 1)],
    'mass_cheb': [('2', 2), ('3', 1), ('1', 1)],
    'mc_sig': [('2,1', 3), ('1,1', 1), ('1,2', 1)],
    'bkg': [('2,1,0', 1), ('1,1,0', 1), ('1,0,1', 1), ('2,0,1', 1), ('1,1,1', 1), ('2,1,1', 1)],
    'bkg_wide_tau': [('0', 1), ('1', 1)],
    'bkg_core_max': [('', 1), ('3', 1)],
    'bkg_profile_frac': [('0', 2), ('1', 1)],
    'window': [('-0.5,3.0', 2), ('-1.0,4.0', 1), ('-2.0,6.0', 1), ('-5.0,10.0', 1), ('-20.0,20.0', 1)],
}


def draw(rng):
    kv = {}
    for k, opts in AXES.items():
        v = rng.choices([o for o, _ in opts], [w for _, w in opts])[0]
        if v == '': continue
        if k == 'mc_sig': kv['mc_cb'], kv['mc_gauss'] = v.split(',')
        elif k == 'bkg': kv['bkg_ss'], kv['bkg_flip'], kv['bkg_ds'] = v.split(',')
        elif k == 'window': kv['ct_low'], kv['ct_high'] = v.split(','); kv['bkg_saved_ct'] = v
        else: kv[k] = v
    return kv


def make_work(work, b, kv):
    if work.exists(): shutil.rmtree(work)
    (work / 'study').mkdir(parents=True)
    for l in LINKS: (work / l).symlink_to(BASE / l)
    (work / 'study/rules').symlink_to(BASE / 'study/rules')
    shutil.copy(BASE / 'study/production.py', work / 'study/production.py')
    shutil.copytree(BASE / 'config', work / 'config')
    y, tag = rr.ytag(b)
    for st in ('mc_mass', 'ctau_pr', 'ctau_np', 'mass', 'ctau_bkg'):
        (work / 'roots' / y / st).mkdir(parents=True)
        for f in (BASE / 'roots' / y / st).glob(f'*{tag}*'):
            if f.is_file(): shutil.copy(f, work / 'roots' / y / st / f.name)
    (work / 'logs').mkdir()
    path = work / 'config/prefit_settings.txt'
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    base = next(r for r in rows if r['bin'] == b)
    changed = {k: v for k, v in kv.items() if base.get(k, '') != v}
    for r in rows:
        if r['bin'] == b: r.update(kv)
    fields = list(rows[0]) + [k for k in kv if k not in rows[0]]
    with open(path, 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n', restval=''); w.writeheader(); w.writerows(rows)
    return changed


def audit(work, b):
    """Per-stage (status, covQual, chi2/ndf) from prefit_audit.txt."""
    f = work / 'prefit_audit.txt'
    if not f.exists(): return {}
    lines = f.read_text().splitlines()
    head = f'{b.split("_")[0]} {b.split("_")[1]}-{b.split("_")[2]}'
    if head not in lines: return {}
    out = {}
    for l in lines[lines.index(head) + 1:]:
        if l and not l.startswith(' '): break
        p = l.split()
        if len(p) > 3 and 'chi2/ndf' in l:
            out[p[0]] = (p[1], p[2], float(l.split('chi2/ndf')[1].split()[0]))
    return out


def prefit(work, b, kv):
    changed = make_work(work, b, kv)
    (work / 'CHANGES.json').write_text(json.dumps(changed))
    if not changed:
        stages = []
    else:
        todo = {s for k in changed for s in stages_from(FIRST[k])}  # union over every changed key (not only the earliest stage)
        stages = [s for s in ORDER if s in todo and s != '2d']
    reg, lo, hi = b.split('_')
    rc = 0
    if stages:
        env = dict(os.environ, OO_MACRO_BASE=str(BASE))
        with open(work / 'prefit.log', 'w') as log:
            rc = subprocess.run(['python3', str(work / 'study/production.py'), *stages, '--rapidity', reg, '--pt', f'{lo}:{hi}', '--jobs', '1'],
                                cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    subprocess.run([ROOT, '-l', '-b', '-q', f'{BASE}/study/prefit_audit.C+("roots","prefit_audit.txt")'], cwd=work, capture_output=True)
    return rc


def fit(work, b, rule, lib):
    reg, lo, hi = b.split('_'); ylo, yhi = ('1.6', '2.4') if reg == 'fwd' else ('0', '1.6')
    out = work / 'rules' / rule; (out / 'logs').mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith('OO_')}
    # the final-fit ctau window follows the combo (ct_low/ct_high in the work config), not rules.window()
    rdir = work / 'rules_knobs' / rule; rdir.mkdir(parents=True, exist_ok=True)
    y, t = rr.ytag(b)
    row = next(r for r in csv.DictReader(open(work / 'config/prefit_settings.txt'), delimiter='\t') if r['bin'] == b)
    text = ''.join(l for l in (BASE / 'study/rules' / rule / f'{t}.txt').read_text().splitlines(True) if not l.startswith(('ct_low ', 'ct_high ')))
    (rdir / f'{t}.txt').write_text(text + f"ct_low {row['ct_low']}\nct_high {row['ct_high']}\n")
    env.update(ROOT_HIST='0', OMP_NUM_THREADS='1', OO_RULE_DIR=str(rdir))
    call = f'fit2d({lo},{hi},{ylo},{yhi},false,false,false,\\"{out}\\")'
    expr = f'gSystem->Load("{lib}");Int_t e=0;gROOT->ProcessLine("{call};",&e);gSystem->Exit(e?11:0);'
    with open(out / 'logs' / f'{b}.log', 'w') as log:
        rc = subprocess.run(['bash', '-c', f'ulimit -f 2000000; exec {ROOT} -n -l -b -q -e {json.dumps(expr)}'], cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    y, t = rr.ytag(b)
    r = rr.result_at(out / 'roots' / y / 'fit2d' / f'fit2d_result_{t}.root', b) if rc == 0 else None
    if r is None:
        why = next((l.strip()[:80] for l in open(out / 'logs' / f'{b}.log', errors='replace') if 'what():' in l), f'rc={rc}')
        return dict(rule=rule, bin=b, fail=why)
    return dict(r, rule=rule)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('name'); ap.add_argument('bin')
    ap.add_argument('--n', type=int, default=100); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--jobs', type=int, default=24); ap.add_argument('--rules', default='F_tail0.1_free_long_free1')
    ap.add_argument('--summary', action='store_true')
    a = ap.parse_args()
    rr.nominal()
    b = a.bin; top = BASE / 'study/scans' / f'{a.name}_{b}'
    res_path = top / 'results.tsv'
    if a.summary:
        rows = list(csv.DictReader(open(res_path), delimiter='\t'))
        for r in sorted(rows, key=lambda r: abs(float(r['pct'])) if r['pct'] else 999): print('\t'.join(r.values()))
        return
    rule_names = [r for r in a.rules.split(',') if r]
    for n in rule_names:
        if n in rules.RULES: rules.write(n)
    lib = next((BASE / 'build').rglob('fit2d_C.so'))
    rng = random.Random(a.seed)
    combos = [{}] + [draw(rng) for _ in range(a.n)]
    top.mkdir(parents=True, exist_ok=True)
    (top / 'combos.json').write_text(json.dumps(combos, indent=0))
    print(f'{len(combos)} combos x {len(rule_names)} rules, {a.jobs} parallel', flush=True)

    def one(i):
        work = top / f'c{i}'
        rc = prefit(work, b, combos[i])
        au = audit(work, b)
        changed = json.loads((work / 'CHANGES.json').read_text())
        rows = []
        for rule in rule_names:
            r = fit(work, b, rule, lib) if rc == 0 else dict(rule=rule, bin=b, fail=f'prefit rc={rc}')
            worst = max((v[2] for v in au.values()), default=float('nan'))
            row = dict(combo=i, changes=json.dumps(changed, sort_keys=True), rule=rule,
                       B=f'{r.get("B", float("nan")):.5f}', Berr=f'{r.get("Berr", float("nan")):.5f}',
                       pct=f'{r["pct"]:.3f}' if 'pct' in r else '', pull=f'{r["pull"]:.3f}' if 'pull' in r else '',
                       status=r.get('status', ''), covQual=r.get('covQual', ''), massP=r.get('massP', ''), timeP=r.get('timeP', ''),
                       prefit_worst_chi2=f'{worst:.3f}', prefit=';'.join(f'{k}:{v[0]}{v[1]}:{v[2]:.2f}' for k, v in au.items()), fail=r.get('fail', ''))
            rows.append(row)
            print(f'c{i:<4} {rule:<32} ' + (f'{r["pct"]:+7.2f}% {r["pull"]:+5.2f}σ st{r["status"]} cq{r["covQual"]}' if 'pct' in r else f'FAIL {r.get("fail")}') + f'  worst chi2 {worst:.2f}  {changed}', flush=True)
        return rows

    fields = ['combo', 'changes', 'rule', 'B', 'Berr', 'pct', 'pull', 'status', 'covQual', 'massP', 'timeP', 'prefit_worst_chi2', 'prefit', 'fail']
    with open(res_path, 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n'); w.writeheader(); fh.flush()
        with ThreadPoolExecutor(a.jobs) as pool:
            for rows in pool.map(one, range(len(combos))):
                w.writerows(rows); fh.flush()
    print('done', res_path)


if __name__ == '__main__':
    main()
