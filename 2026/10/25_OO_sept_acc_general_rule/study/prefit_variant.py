#!/usr/bin/env python3
"""Refit one bin with changed prefit choices, from the first affected stage through the final 2D fit,
in an isolated work directory (the main roots/ and config/ are not touched).

  python3 study/prefit_variant.py NAME BIN key=value [key=value ...] [--rule=FINAL_S1]
      keys are config/prefit_settings.txt columns (pr_quant, pr_nres, np_count, mass_cheb, bkg_ss, ...)
  python3 study/prefit_variant.py --summary

Work directory: study/prefit_variants/<NAME>_<BIN>/ (macros, utils, inputs linked; own config, roots, logs).
"""
from pathlib import Path
import csv, os, shutil, subprocess, sys
BASE = Path(__file__).resolve().parents[1]
ORDER = ['mc', 'pr', 'np', 'mass', 'bkg', '2d']
FIRST = {'mc_cb': 'mc', 'mc_gauss': 'mc', 'pr_nres': 'pr', 'pr_quant': 'pr', 'np_count': 'np', 'mass_cheb': 'mass',
         'ct_low': 'mass', 'ct_high': 'mass', 'bkg_ss': 'bkg', 'bkg_flip': 'bkg', 'bkg_ds': 'bkg', 'bkg_profile_frac': 'bkg', 'bkg_saved_ct': 'bkg', 'bkg_wide_tau': 'bkg', 'bkg_core_max': 'bkg'}
# stage -> stages that read its output
NEEDS = {'mc': {'mass'}, 'pr': {'np', 'bkg'}, 'np': set(), 'mass': {'bkg'}, 'bkg': set()}
LINKS = ['mc_mass.C', 'ctau_pr.C', 'ctau_np.C', 'mass_window_screen.C', 'mass.C', 'ctau_bkg.C', 'fit2d.C', 'utils', 'include', 'inputs', 'workflow', 'build']


def stages_from(first):
    todo, frontier = {first}, [first]
    while frontier:
        for n in NEEDS.get(frontier.pop(), ()):
            if n not in todo: todo.add(n); frontier.append(n)
    return [s for s in ORDER if s in todo] + ['2d']


def main():
    if sys.argv[1:2] == ['--summary']:
        sys.path.insert(0, str(BASE / 'study')); import run_rules as rr
        rr.nominal()
        for d in sorted(x for x in (BASE / 'study/prefit_variants').glob('*') if x.is_dir()):
            b = '_'.join(d.name.split('_')[-3:])
            y, t = rr.ytag(b)
            r = rr.result_at(d / 'roots' / y / 'fit2d' / f'fit2d_result_{t}.root', b)
            if r: print(rr.fmt(dict(r, rule=d.name))); continue
            why = [l.strip() for f in d.glob('logs/*/*.log') for l in open(f, errors='replace') if 'what():' in l]
            if not why: print(f'{d.name:<34} {b:<10} running'); continue
            audit = d / 'prefit_audit.txt'
            bad = ''
            if audit.exists():
                lines = audit.read_text().splitlines()
                head = f'{b.split("_")[0]} {b.split("_")[1]}-{b.split("_")[2]}'
                i = lines.index(head)
                bad = '; '.join(l.split('  res')[0].strip() for l in lines[i + 1:i + 6] if 'chi2/ndf' in l and float(l.split('chi2/ndf')[1].split()[0]) >= 2)
            print(f'{d.name:<34} {b:<10} REJECTED {why[0][:90]}' + (f'  [{bad}]' if bad else ''))
        return
    args = [a for a in sys.argv[1:] if not a.startswith('--rule')]
    rule = next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--rule=')), 'FINAL_S1')
    name, b, kv = args[0], args[1], dict(a.split('=', 1) for a in args[2:])
    first = min((FIRST[k] for k in kv), key=ORDER.index)
    work = BASE / 'study/prefit_variants' / f'{name}_{b}'
    if work.exists(): shutil.rmtree(work)
    (work / 'study').mkdir(parents=True)
    for l in LINKS: (work / l).symlink_to(BASE / l)
    (work / 'study/rules').symlink_to(BASE / 'study/rules')
    shutil.copy(BASE / 'study/production.py', work / 'study/production.py')
    shutil.copytree(BASE / 'config', work / 'config')
    reg, lo, hi = b.split('_'); y = 'y1.60_2.40' if reg == 'fwd' else 'y0.00_1.60'
    tag = f'{y}_pt{float(lo):.2f}_{float(hi):.2f}'
    for st in ('mc_mass', 'ctau_pr', 'ctau_np', 'mass', 'ctau_bkg'):  # this bin's current prefits
        (work / 'roots' / y / st).mkdir(parents=True)
        for f in (BASE / 'roots' / y / st).glob(f'*{tag}*'):
            if f.is_file(): shutil.copy(f, work / 'roots' / y / st / f.name)
    (work / 'logs').mkdir()
    path = work / 'config/prefit_settings.txt'
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    for r in rows:
        if r['bin'] == b: r.update(kv)
    fields = list(rows[0]) + [k for k in kv if k not in rows[0]]
    with open(path, 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n', restval=''); w.writeheader(); w.writerows(rows)
    stages = stages_from(first)
    (work / 'VARIANT.txt').write_text(f'bin {b}\nrule {rule}\nchanges {kv}\nstages {stages}\n')
    env = dict(os.environ, OO_MACRO_BASE=str(BASE), OO_FINAL_RULE_DIR=str(BASE / 'study/rules' / rule))
    rc = subprocess.run(['python3', str(work / 'study/production.py'), *stages, '--rapidity', reg, '--pt', f'{lo}:{hi}', '--jobs', '1'],
                        cwd=work, env=env).returncode
    subprocess.run([str(Path.home() / 'build/root-6.36.06/bin/root'), '-l', '-b', '-q', f'{BASE}/study/prefit_audit.C+("roots","prefit_audit.txt")'],
                   cwd=work, capture_output=True)
    sys.exit(rc)


if __name__ == '__main__':
    main()
