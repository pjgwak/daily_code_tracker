#!/usr/bin/env python3
"""Production runner (September folder layout): prefits and the final 2D fit into roots/, figs/, figs_publish/.

Called by run_ctau_bins.sh / run_publish.sh. Per-bin prefit choices come from config/prefit_settings.txt;
the 2D fit uses the knob files in study/rules/FINAL_S3 (OO_RULE_DIR).
  stages : mc, pr, np, mass, bkg, 2d (in this order when several are given). Bins run in parallel and each
           bin moves to its next stage as soon as its own previous stage finishes (a failed stage stops that bin).
  --redraw  : redraw saved results into figs/ without refitting
  --publish : redraw saved results into figs_publish/ without refitting
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse, csv, datetime, hashlib, json, os, subprocess, sys, threading
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
    if s.get('mass_lo'): e['OO_MASS_LO'] = s['mass_lo']  # dimuon mass fit window (utils/mass_window.h), every stage
    if s.get('mass_hi'): e['OO_MASS_HI'] = s['mass_hi']
    if stage == 'mc': e['OO_MC_SIG'] = f'{s["mc_cb"]},{s["mc_gauss"]}'
    if stage == 'pr': e.update(OO_PR_NRES=s['pr_nres'], OO_PR_QUANT=s.get('pr_quant', '0.0005,0.9995'))
    if stage == 'np' and s.get('np_par'): e['OO_NP_PAR'] = s['np_par']  # per-bin start value / range overrides of the truth-ctau fit
    if stage == 'pr' and s.get('pr_par'): e['OO_PR_PAR'] = s['pr_par']  # per-bin start value / range overrides of the MC resolution fit (name=init,min,max;...)
    if stage == 'np': e['OO_NP_COUNT'] = s['np_count']
    if stage == 'np' and s.get('np_range'): e['OO_NP_RANGE'] = s['np_range']  # truth-ctau fit range, default 0,6
    if stage == 'mass': e['OO_MASS_CHEB'] = s['mass_cheb']
    if stage == 'bkg': e.update(OO_BKG_COUNTS=f'{s["bkg_ss"]},{s["bkg_flip"]},{s["bkg_ds"]}', OO_BKG_CT_LO=s['ct_low'], OO_BKG_CT_HI=s['ct_high'], OO_BKG_ROBUST='1', OO_BKG_PROFILE_FRAC=s.get('bkg_profile_frac', '0'))
    if stage == 'bkg' and s.get('bkg_wide_tau', '0') == '1': e['OO_BKG_WIDE_TAU'] = '1'
    if stage == 'bkg' and s.get('bkg_core_max'): e['OO_BKG_CORE_MAX'] = s['bkg_core_max']
    if stage == 'bkg' and s.get('bkg_par'): e['OO_BKG_PAR'] = s['bkg_par']  # per-bin start value / range overrides (name=init,min,max;...)
    if stage == 'bkg' and s.get('bkg_plot_bins'): e['OO_BKG_PLOT_BINS'] = s['bkg_plot_bins']  # display / chi2 binning of the sideband ctau fit
    if stage == '2d': e['OO_RULE_DIR'] = os.environ.get('OO_FINAL_RULE_DIR', str(BASE / 'study/rules/FINAL_S3'))  # prefit_variant.py --rule
    return e


PR_MAP_LOCK = threading.Lock()


def tag(s):
    return f'{"y1.60_2.40" if s["region"] == "fwd" else "y0.00_1.60"}_pt{float(s["lo"]):.2f}_{float(s["hi"]):.2f}'


def update_pr_map(settings):
    """Refitted PR prefits change their md5: refresh config/pr_calibration_windows.txt (and its window if saved)."""
    with PR_MAP_LOCK:
        _update_pr_map({tag(s) for s in settings})


def _update_pr_map(tags):
    path = BASE / 'config/pr_calibration_windows.txt'
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    for r in rows:
        if r['tag'] not in tags: continue
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
    for r in csv.DictReader(open(BASE / 'config/prefit_settings.txt'), delimiter='\t'):
        reg, lo, hi = r['bin'].split('_')
        r.update(region=reg, lo=lo, hi=hi, ylo='1.6' if reg == 'fwd' else '0', yhi='2.4' if reg == 'fwd' else '1.6')
        if a.rapidity != 'both' and reg != a.rapidity: continue
        if a.pt and a.pt != f'{lo}:{hi}': continue
        settings.append(r)
    draw = a.redraw or a.publish
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    build = BASE / 'build'; build.mkdir(exist_ok=True)
    stages = [s for s in ORDER if s in a.stages]
    libs, logdirs = {}, {}
    for stage in stages:
        macro = Path(os.environ.get('OO_MACRO_BASE', BASE)) / MACRO[stage]  # isolated work dirs compile against the main sources
        r = subprocess.run([ROOT, '-n', '-l', '-b', '-q', '-e', f'gSystem->SetBuildDir("{build}",true);gSystem->Exit(gSystem->CompileMacro("{macro}","kO")?0:10);'], cwd=BASE, capture_output=True, text=True)
        if r.returncode: sys.exit(f'compile failed: {macro}\n{r.stdout[-2000:]}')
        libs[stage] = next(build.rglob(macro.name.replace('.C', '_C.so')))
        logdirs[stage] = BASE / 'logs' / f'{stage}{"_publish" if a.publish else "_redraw" if a.redraw else ""}_{stamp}'; logdirs[stage].mkdir(parents=True, exist_ok=True)

    def run(stage, s):
        env = {k: v for k, v in os.environ.items() if not k.startswith('OO_') or k == 'OO_BKG_CHI2_SCAN' or k == 'OO_PREFIT_CHI2_MAX' or k == 'OO_DIAGNOSTICS_DIR' or k.startswith('OO_PREFIT_JITTER')}  # diagnostics and restart flags pass through
        env.update(stage_env(stage, s))
        expr = f'gSystem->Load("{libs[stage]}");Int_t e=0;gROOT->ProcessLine("{call(stage, s, draw, a.publish).replace(chr(34), chr(92) + chr(34))};",&e);gSystem->Exit(e?11:0);'
        with open(logdirs[stage] / f'{s["bin"]}.log', 'w') as log:
            rc = subprocess.run(['bash', '-c', f'ulimit -f {LOG_LIMIT_KB}; exec {ROOT} -n -l -b -q -e {json.dumps(expr)}'], cwd=BASE, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
        if stage == 'pr' and not draw and rc == 0: update_pr_map([s])  # this bin's 2d checks the PR fingerprint
        return rc

    def chain(s):
        done = []
        for stage in stages:
            rc = run(stage, s)
            t = datetime.datetime.now().strftime('%H:%M:%S')
            print(f'{t} {s["bin"]} [{stage}] ' + ('ok' if rc == 0 else f'FAILED (rc {rc}), log {logdirs[stage].relative_to(BASE)}/{s["bin"]}.log'), flush=True)
            if rc: return s['bin'], done, stage
            done.append(stage)
        return s['bin'], done, None

    with ThreadPoolExecutor(a.jobs) as pool:
        res = list(pool.map(chain, settings))
    for stage in stages:
        ok = sum(stage in d for _, d, _ in res)
        bad = [b for b, _, f in res if f == stage]
        print(f'[{stage}] {ok}/{len(res)} ok' + (f', failed: {bad}' if bad else '') + f'  logs: {logdirs[stage].relative_to(BASE)}', flush=True)
    if '2d' in a.stages and not draw:
        subprocess.run([ROOT, '-n', '-l', '-b', '-q', f'workflow/collect.C(".","drawing/final_summary/summary.txt")'], cwd=BASE, capture_output=True)
        print('summary: drawing/final_summary/summary.txt')


if __name__ == '__main__':
    main()
