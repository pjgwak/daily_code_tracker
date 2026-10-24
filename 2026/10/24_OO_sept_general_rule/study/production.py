#!/usr/bin/env python3
"""Production runner (September folder layout): prefits and the final 2D fit into roots/, figs/, figs_publish/.

Called by run_ctau_bins.sh / run_publish.sh. Per-bin prefit choices come from config/prefit_settings.tsv;
the 2D fit uses the knob files in study/rules/FINAL_v3 (OO_RULE_DIR).
  stages : mc, pr, np, mass, bkg, 2d (in this order when several are given)
  --redraw  : redraw saved results into figs/ without refitting
  --publish : redraw saved results into figs_publish/ without refitting
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse, csv, datetime, hashlib, json, os, subprocess, sys
BASE = Path(__file__).resolve().parents[1]
ROOT = str(Path.home() / 'build/root-6.36.06/bin/root')
MACRO = {'mc': 'mc_mass.C', 'pr': 'ctau_pr.C', 'np': 'ctau_np.C', 'mass': 'mass_window_screen.C', 'bkg': 'ctau_bkg.C', '2d': 'fit2d.C'}
ORDER = ['mc', 'pr', 'np', 'mass', 'bkg', '2d']
LOG_LIMIT_KB = 2_000_000


def call(stage, s, draw, publish):
    lo, hi = s['lo'], s['hi']; ylo, yhi = s['ylo'], s['yhi']
    d, p = str(draw).lower(), str(publish).lower()
    return {'mc': f'mc_mass({lo},{hi},{ylo},{yhi},false,{d},{p})',
            'pr': f'ctau_pr({lo},{hi},{ylo},{yhi},{d},{p})',
            'np': f'ctau_np({lo},{hi},{ylo},{yhi},{d},{p})',
            'mass': f'mass_window_screen({lo},{hi},{ylo},{yhi},{d},{p},{s["ct_low"]},{s["ct_high"]})',
            'bkg': f'ctau_bkg({lo},{hi},{ylo},{yhi},true,{d},{p})',
            '2d': f'fit2d({lo},{hi},{ylo},{yhi},false,{d},{p},"")'}[stage]


def stage_env(stage, s):
    e = {'ROOT_HIST': '0', 'OMP_NUM_THREADS': '1'}
    if stage == 'mc': e['OO_MC_SIG'] = f'{s["mc_cb"]},{s["mc_gauss"]}'
    if stage == 'pr': e['OO_PR_NRES'] = s['pr_nres']
    if stage == 'np': e['OO_NP_COUNT'] = s['np_count']
    if stage == 'mass': e['OO_MASS_CHEB'] = s['mass_cheb']
    if stage == 'bkg': e.update(OO_BKG_COUNTS=f'{s["bkg_ss"]},{s["bkg_flip"]},{s["bkg_ds"]}', OO_BKG_CT_LO=s['ct_low'], OO_BKG_CT_HI=s['ct_high'], OO_BKG_ROBUST='1')
    if stage == '2d': e['OO_RULE_DIR'] = str(BASE / 'study/rules/FINAL_v3')
    return e


def update_pr_map(settings):
    """Refitted PR prefits change their md5: refresh config/pr_calibration_windows.tsv (and its window if saved)."""
    path = BASE / 'config/pr_calibration_windows.tsv'
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    for r in rows:
        y, _ = r['tag'].split('_pt')
        f = BASE / 'roots' / y / 'ctau_pr' / f'ctau_resolution_{r["tag"]}.root'
        if not f.exists(): continue
        r['prefit_md5'] = hashlib.md5(f.read_bytes()).hexdigest()
        out = subprocess.run([ROOT, '-l', '-b', '-q', '-e', f'TFile f("{f}"); auto*a=f.Get<TParameter<double>>("prCalibrationLow"); auto*b=f.Get<TParameter<double>>("prCalibrationHigh"); if(a&&b) printf("RNG %.17g %.17g\\n",a->GetVal(),b->GetVal());'], capture_output=True, text=True).stdout
        rng = next((l.split()[1:] for l in out.splitlines() if l.startswith('RNG')), None)
        if rng: r['ct_low'], r['ct_high'] = rng
    with open(path, 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('stages', nargs='+', choices=ORDER)
    ap.add_argument('--redraw', action='store_true'); ap.add_argument('--publish', action='store_true')
    ap.add_argument('--rapidity', default='both', choices=['both', 'fwd', 'mid']); ap.add_argument('--pt', help='LOW:HIGH, e.g. 8:9')
    ap.add_argument('--jobs', type=int, default=12)
    a = ap.parse_args()
    settings = []
    for r in csv.DictReader(open(BASE / 'config/prefit_settings.tsv'), delimiter='\t'):
        reg, lo, hi = r['bin'].split('_')
        r.update(region=reg, lo=lo, hi=hi, ylo='1.6' if reg == 'fwd' else '0', yhi='2.4' if reg == 'fwd' else '1.6')
        if a.rapidity != 'both' and reg != a.rapidity: continue
        if a.pt and a.pt != f'{lo}:{hi}': continue
        settings.append(r)
    draw = a.redraw or a.publish
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    build = BASE / 'build'; build.mkdir(exist_ok=True)
    for stage in [s for s in ORDER if s in a.stages]:
        macro = BASE / MACRO[stage]
        r = subprocess.run([ROOT, '-n', '-l', '-b', '-q', '-e', f'gSystem->SetBuildDir("{build}",true);gSystem->Exit(gSystem->CompileMacro("{macro}","kO")?0:10);'], cwd=BASE, capture_output=True, text=True)
        if r.returncode: sys.exit(f'compile failed: {macro}\n{r.stdout[-2000:]}')
        lib = next(build.rglob(macro.name.replace('.C', '_C.so')))
        logdir = BASE / 'logs' / f'{stage}{"_publish" if a.publish else "_redraw" if a.redraw else ""}_{stamp}'; logdir.mkdir(parents=True)

        def one(s):
            env = {k: v for k, v in os.environ.items() if not k.startswith('OO_')}
            env.update(stage_env(stage, s))
            expr = f'gSystem->Load("{lib}");Int_t e=0;gROOT->ProcessLine("{call(stage, s, draw, a.publish).replace(chr(34), chr(92) + chr(34))};",&e);gSystem->Exit(e?11:0);'
            with open(logdir / f'{s["bin"]}.log', 'w') as log:
                rc = subprocess.run(['bash', '-c', f'ulimit -f {LOG_LIMIT_KB}; exec {ROOT} -n -l -b -q -e {json.dumps(expr)}'], cwd=BASE, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
            return s['bin'], rc

        with ThreadPoolExecutor(a.jobs) as pool:
            res = list(pool.map(one, settings))
        bad = [b for b, rc in res if rc]
        print(f'[{stage}] {len(res) - len(bad)}/{len(res)} ok' + (f', failed: {bad}' if bad else '') + f'  logs: {logdir.relative_to(BASE)}', flush=True)
        if stage == 'pr' and not draw: update_pr_map(settings)
    if '2d' in a.stages and not draw:
        subprocess.run([ROOT, '-n', '-l', '-b', '-q', f'workflow/collect.C(".","drawing/final_summary/summary.tsv")'], cwd=BASE, capture_output=True)
        print('summary: drawing/final_summary/summary.tsv')


if __name__ == '__main__':
    main()
