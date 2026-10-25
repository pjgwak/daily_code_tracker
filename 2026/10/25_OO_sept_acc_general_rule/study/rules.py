#!/usr/bin/env python3
"""Write knob files for simple component rules (fix / constrain / free per PDF component).

  python3 study/rules.py A B C      -> study/rules/<rule>/<yTag>_pt<lo>_<hi>.txt

A rule is a dict of component choices; a rule may differ only between named
regions (REGION below), never per bin. Nominal B is not an input.
Component vocabulary (fit2d knob translation in knobs()):
  pr_mean : 'free' | 'fix' | 'mc8'              (fix: MC prefit value; mc8: constrained at MC, width 8 x MC error)
  pr_frac : 'fix' | 'mc'                        (MC covariance x pr_frac_k)
  pr_width: 'abs' | 'tail' | 'free' | 'mc_joint' (abs: sigma=(.010,.020,.040,.080) mm x pr_width_k on
                                                 (sigma1, gaps); tail: core widths fixed at MC, last two
                                                 gaps constrained as abs; mc_joint: mean+widths+fractions with
                                                 MC covariance x pr_joint_k, overrides pr_mean/pr_frac)
  np      : 'mc' | 'long_free' | 'free'          (MC covariance x np_k on constrained coordinates)
  bkg     : 'sideband' | 'free'                  (sideband prefit constraint x bkg_k, or floating)
"""
from pathlib import Path
import sys

BASE = Path(__file__).resolve().parents[1]
EDGES = {'fwd': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 20], 'mid': [7, 8, 9, 10, 12, 14, 16, 20, 30]}

DEFAULT = dict(pr_mean='free', pr_frac='fix', pr_frac_k=1, pr_width='abs', pr_width_k=1, pr_joint_k=1,
               np='long_free', np_k=1, bkg='free', bkg_k=1)

# rule -> {region: overrides}; region keys are tried in order 'fwd_low' (fwd pT<5), 'fwd', 'mid', 'all'.
RULES = {
    'A': {'all': {}},
    'B': {'all': dict(bkg='sideband')},
    'C': {'all': dict(pr_width='mc_joint', np='mc', bkg='sideband')},
    'D': {'all': dict(np='mc')},
}
# Grid G_<abs|tail><k>_<frac>_<np>_<bkg>: one choice for all 20 bins.
for w in ('abs', 'tail'):
    for k in (0.1, 0.25, 0.5, 1, 2):
        for frac in ('fix', 'mc'):
            for np_ in ('long_free', 'mc', 'free'):
                for bkg in ('free', 'sideband'):
                    RULES[f'G_{w}{k}_{frac}_{np_}_{bkg}'] = {'all': dict(pr_width=w, pr_width_k=k, pr_frac=frac, np=np_, bkg=bkg)}

# Grid H (forward 1-3 only): constraint strengths around the September low-pT treatment.
for wk in (2, 4, 8):
    for fk in (1, 16):
        for np_, nk in (('mc', 1), ('mc', 16), ('mc', 175), ('free', 1)):
            for bk in (1, 4):
                RULES[f'H_w{wk}_f{fk}_{np_}{nk}_b{bk}'] = {'all': dict(pr_width_k=wk, pr_frac='mc', pr_frac_k=fk, np=np_, np_k=nk, bkg='sideband', bkg_k=bk)}

# Grid P: PR mean fixed at MC (B correlates with the PR mean in mid 7-10).
for k in (0.1, 0.25, 1):
    for frac in ('fix', 'mc'):
        for np_ in ('long_free', 'mc'):
            RULES[f'P_meanfix_abs{k}_{frac}_{np_}'] = {'all': dict(pr_mean='fix', pr_width_k=k, pr_frac=frac, np=np_)}

# Grid F (forward 5-20 only): core widths fixed (tail) with mean / NP strength choices.
for k in (0.05, 0.1, 0.25):
    for mean in ('free', 'mc8', 'fix'):
        for np_, nk in (('long_free', 1), ('long_free', 16), ('mc', 1)):
            RULES[f'F_tail{k}_{mean}_{np_}{nk}'] = {'all': dict(pr_width='tail', pr_width_k=k, pr_mean=mean, np=np_, np_k=nk)}

# Signal purity r = Nsig/Nbkg of the DATA mass prefit (config/signal_purity.txt).
import csv as _csv
PURITY = {r['bin']: float(r['SB']) for r in _csv.DictReader(open(BASE / 'config' / 'signal_purity.txt'), delimiter='\t')}
# MC prompt width-gap ratios q3=(s3-s2)/s1, q4=(s4-s3)/s1 from the PR prefits (config/pr_gap_ratios.txt).
GAPS = {r['bin']: (int(r['nres']), float(r['q3']), float(r['q4'])) for r in _csv.DictReader(open(BASE / 'config' / 'pr_gap_ratios.txt'), delimiter='\t')}

# Grid K (all 20 bins): gap-ratio constraints centred at k x MC (sigma 5% / 15% of the centre), fractions fixed.
for kk in (1, 1.25, 1.5, 1.75, 2):
    for np_ in ('mc', 'long_free'):
        RULES[f'K_gap{kk}_{np_}'] = {'all': dict(pr_width_k=1, gap_k=kk, np=np_)}

# Grid S5 (forward 5-10 only): core-width constraint grows with purity, s1 = c * r^2 (s2 = 2 s1, outer 0.004 / 0.08 mm).
for cc in (0.0006, 0.0008, 0.0011, 0.0015):
    for mean in ('free', 'mc8'):
        for nk in (1, 16):
            RULES[f'S5_c{cc}_{mean}_lf{nk}'] = {'all': dict(pr_width_purity=cc, pr_mean=mean, np='long_free', np_k=nk)}

# Grid C5 (forward 5-10 only): one core-width constraint s for all bins, September structure
# (sigma1: s, first gap: 2s, outer gaps 0.004 / 0.08 mm).
for core in (0.0002, 0.0003, 0.0005, 0.0008, 0.0012):
    for nk in (1, 16):
        RULES[f'C5_core{core}_lf{nk}'] = {'all': dict(pr_width_sig=[core, 2 * core, .004, .08], np='long_free', np_k=nk)}

# Grid T35 (forward 3-5 only): PR width-constraint strength between x0.1 and x0.25.
for k in (0.12, 0.15, 0.18, 0.2):
    for mean in ('free', 'mc8'):
        for nk in (1, 4):
            RULES[f'T35_abs{k}_{mean}_mc{nk}'] = {'all': dict(pr_width_k=k, pr_mean=mean, np='mc', np_k=nk, bkg='sideband')}

# Release tests for the forward 5-10 core constraint.
RULES['U_corefree'] = {'all': dict(pr_width_sig=[1.0, 1.0, .004, .08], np='long_free')}
RULES['U_widthfree'] = {'all': dict(pr_width='free', np='long_free')}

# Grid K58 (forward 5-8 only): core widths fixed at MC (tail x0.1), other components varied.
for frac in ('fix', 'mc'):
    for np_ in ('long_free', 'free'):
        for bkg in ('free', 'sideband'):
            RULES[f'K58_{frac}_{np_}_{bkg}'] = {'all': dict(pr_width='tail', pr_width_k=0.1, pr_frac=frac, np=np_, bkg=bkg)}

# Grid Y: background parameterization in the final fit (fractions vs prefit-style yields), same physics otherwise.
for base, kw in (('f58', dict(pr_width='tail', pr_width_k=0.1)), ('f13', dict(pr_width_k=8, pr_frac='mc', np='mc', np_k=175))):
    RULES[f'Y_{base}_frac'] = {'all': dict(kw, bkg='free')}
    for bk in (1, 1000):
        RULES[f'Y_{base}_yield{bk}'] = {'all': dict(kw, bkg='sideband', bkg_k=bk)}

# Grid J (mid 7-10): PR mean, widths and fractions jointly constrained with MC covariance x k.
for k in (4, 16, 64, 256):
    for np_ in ('long_free', 'mc', 'free'):
        RULES[f'J_joint{k}_{np_}'] = {'all': dict(pr_width='mc_joint', pr_joint_k=k, np=np_)}

# Ladder L58 (forward 5-10): core width fixed at 5-7, loose core constraint at 7-8, ordered widths x0.25 at 8-10.
for s78 in (0.0003, 0.0004, 0.0005):
    RULES[f'L58_s{s78}'] = {'fwd:5-7': dict(pr_width='tail', pr_width_k=0.1),
                            'fwd:7-8': dict(pr_width_sig=[s78, 2 * s78, .004, .08]),
                            'fwd:8-10': dict(pr_width_k=0.25)}

# Grid N23 (forward 2-3): NP shape nearly fixed to MC (covariance x0.01 / x0.1) with PR choices.
for nk in (0.01, 0.1):
    for frac in ('fix', 'mc'):
        for w, wk in (('tail', 0.1), ('abs', 0.25), ('abs', 2)):
            RULES[f'N23_np{nk}_{frac}_{w}{wk}'] = {'all': dict(np='mc', np_k=nk, pr_frac=frac, pr_width=w, pr_width_k=wk, bkg='sideband')}

# Grid NP (mid 7-12): NP lifetime treatment in the final fit, two PR treatments.
for pr, kw in (('w025', dict(pr_width_k=0.25)), ('j64', dict(pr_width='mc_joint', pr_joint_k=64))):
    for npn, npkw in (('mc', dict(np='mc')), ('lf1', dict(np='long_free')), ('lf16', dict(np='long_free', np_k=16)),
                      ('alltau', dict(np='mc', extra=dict(np_all_lifetimes_free=1))),
                      ('moment', dict(np='mc', extra=dict(np_moment_order=1))), ('free', dict(np='free'))):
        RULES[f'NPT_{pr}_{npn}'] = {'all': dict(kw, **npkw)}

# Grid E (all 20 bins): PR widths and fractions constrained at MC with total width sqrt(stat^2 + model^2),
# model = largest shift of the MC PR prefit under alternative fit windows (config/pr_model_err.txt).
def _model_err():
    path = BASE / 'config' / 'pr_model_err.txt'
    if not path.exists(): return {}
    out = {}
    for r in _csv.DictReader(open(path), delimiter='\t'): out.setdefault(r['bin'], {})[r['parameter']] = float(r['total'])
    return out
MODEL_ERR = _model_err()
for np_ in ('mc', 'long_free'):
    for bkg in ('free', 'sideband'):
        RULES[f'E_{np_}_{bkg}'] = {'all': dict(pr_width='model', pr_frac='model', np=np_, bkg=bkg)}

# Grid W35 (forward 3-5, 2026-10-25 Sept-acc study): PR width constraints tighter than x0.1, or core widths fixed.
for k in (0.01, 0.02, 0.05):
    for mean in ('free', 'fix'):
        for nk in (1, 4):
            RULES[f'W35_abs{k}_{mean}_mc{nk}'] = {'all': dict(pr_width_k=k, pr_mean=mean, np='mc', np_k=nk, bkg='sideband')}
for k in (0.1, 0.25):
    for mean in ('free', 'fix'):
        RULES[f'W35_tail{k}_{mean}'] = {'all': dict(pr_width='tail', pr_width_k=k, pr_mean=mean, np='mc', bkg='sideband')}

# Grid KL (2026-10-25): gap-ratio centres below MC (narrower PR tail than MC), fractions fixed.
for kk in (0.5, 0.65, 0.8, 0.9):
    for np_ in ('mc', 'long_free'):
        RULES[f'KL_gap{kk}_{np_}'] = {'all': dict(pr_width_k=1, gap_k=kk, np=np_)}
        RULES[f'KL_tailgap{kk}_{np_}'] = {'all': dict(pr_width='tail', pr_width_k=0.1, gap_k=kk, np=np_)}

for kk in (0.2, 0.3, 0.4):
    RULES[f'KL_tailgap{kk}_long_free'] = {'all': dict(pr_width='tail', pr_width_k=0.1, gap_k=kk, np='long_free')}
for kk in (0.3, 0.4, 0.5):
    RULES[f'KL_tailgap{kk}_long_free16'] = {'all': dict(pr_width='tail', pr_width_k=0.1, gap_k=kk, np='long_free', np_k=16)}
    RULES[f'KL_tailgap{kk}_long_free_sb'] = {'all': dict(pr_width='tail', pr_width_k=0.1, gap_k=kk, np='long_free', bkg='sideband')}

# FINAL_v1: assembled group rules (study of 2026-10-24/25). Range keys are lo <= pT < hi.
RULES['FINAL_v1'] = {
    'fwd:1-3': dict(pr_width_k=8, pr_frac='mc', np='mc', np_k=175, bkg='sideband'),
    'fwd:3-5': dict(pr_width_k=0.12, np='mc', bkg='sideband'),
    'fwd:5-7': dict(pr_width='tail', pr_width_k=0.1),
    'fwd:7-8': dict(pr_width_sig=[0.0004, 0.0008, .004, .08]),
    'fwd:8-10': dict(pr_width_k=0.25),
    'fwd:10-20': dict(pr_width='tail', pr_width_k=0.1, np='long_free', np_k=16),
    'mid:7-8': dict(pr_width='mc_joint', pr_joint_k=64, np='mc'),
    'mid:8-12': dict(pr_width_k=0.25, pr_frac='mc', np='mc'),
    'mid:12-30': dict(pr_width_k=0.25, pr_frac='mc', np='long_free'),
}

# September forward 1-3 settings (09/19 README): width sd x2 / x4, fraction cov x16 / x8, NP cov x175 / x700, bkg x4.
RULES['SEPT12'] = {'all': dict(pr_width_k=2, pr_frac='mc', pr_frac_k=16, np='mc', np_k=175, bkg='sideband', bkg_k=4)}
RULES['SEPT23'] = {'all': dict(pr_width_k=4, pr_frac='mc', pr_frac_k=8, np='mc', np_k=700, bkg='sideband', bkg_k=4)}

# FINAL_S1: rules for the September single-muon-acceptance inputs (2026-10-25).
# fwd 2-3 is a per-bin setting written by hand (FINAL_v3 fwd 2-3 structure with gap-ratio centres 1.8 x MC,
# study/rules/C23_gap1.8); write_final_s1() copies it over the generated file.
RULES['FINAL_S1'] = {
    'fwd:1-2': dict(pr_width_k=2, pr_frac='mc', pr_frac_k=16, np='mc', np_k=16, bkg='sideband', bkg_k=4),
    'fwd:2-3': dict(pr_width_k=2, pr_frac='mc', pr_frac_k=16, np='mc', np_k=16, bkg='sideband', bkg_k=4),  # replaced, see above
    'fwd:3-5': dict(pr_width_k=0.05, np='mc', bkg='sideband'),
    'fwd:5-8': dict(pr_width='tail', pr_width_k=0.25, np='long_free'),
    'fwd:8-10': dict(pr_width_sig=[0.0012, 0.0024, .004, .08], np='long_free'),
    'fwd:10-20': dict(pr_width='tail', pr_width_k=0.1, gap_k=0.4, np='long_free'),
    'mid:7-30': dict(pr_mean='fix', pr_width_k=1, pr_frac='fix', np='long_free'),
}

# Diagnostics for the mid 8-10 background study.
RULES['X_profile_mass'] = {'all': dict(pr_width_k=0.25, extra=dict(profile_mass=1))}
for m in (1, 2, 3):
    RULES[f'X_bkgres{m}'] = {'all': dict(pr_width_k=0.25, extra=dict(background_core_profile=m))}

# ctau fit window from the 23 DATA PR/NP boundary study (PDF not visibly cut at the window edges).
def window(r, lo):
    if r == 'mid': return (-1.0, 4.0)
    return (-4.0, 6.5) if lo < 3 else (-3.5, 4.5) if lo < 5 else (-1.0, 4.0) if lo < 10 else (-0.5, 3.0)


def region_keys(r, lo):
    return (['fwd_low'] if r == 'fwd' and lo < 5 else []) + [r, 'all']


def choice(rule, r, lo):
    c = dict(DEFAULT)
    for key in reversed(region_keys(r, lo)):
        c.update(RULES[rule].get(key, {}))
    # pT-range keys such as 'fwd:5-7' (lo <= pT < hi) apply last, so they are the most specific.
    for key, v in RULES[rule].items():
        if ':' in key:
            reg, rng = key.split(':'); a, b = map(float, rng.split('-'))
            if reg == r and a <= lo < b: c.update(v)
    return c


def knobs(c):
    k = {}
    if c['pr_width'] == 'mc_joint':
        k.update(prompt_full_cov=1, width_policy=3, prompt_full_cov_scale=c['pr_joint_k'], float_prompt_mean=1,
                 prompt_frac_constrain=0)
    else:
        k.update(prompt_full_cov=0, width_policy=1 if c['pr_width'] == 'tail' else 0, float_prompt_mean=int(c['pr_mean'] != 'fix'),
                 prompt_frac_constrain=int(c['pr_frac'] in ('mc', 'model')), prompt_fraction_cov_scale=c['pr_frac_k'],
                 prompt_mean_prior_mult=8 if c['pr_mean'] == 'mc8' else 0, prompt_width_constraint=int(c['pr_width'] in ('abs', 'tail', 'model')))
        if c['pr_frac'] == 'model':
            e = MODEL_ERR[c['bin']]
            k.update(prompt_frac_constrain=1, **{f'fraction_sigma{i}': e[f'ctauFrac{i}Var'] for i in (1, 2, 3) if f'ctauFrac{i}Var' in e})
        if c.get('pr_width_purity'):
            s1 = c['pr_width_purity'] * c['purity'] ** 2
            c['pr_width_sig'] = [s1, 2 * s1, .004, .08]
        if c['pr_width'] == 'model':
            e = MODEL_ERR[c['bin']]
            c['pr_width_sig'] = [e.get(n, 1.0) for n in ('ctauTime1Scale', 'ctauTime2Delta', 'ctauTime3Delta', 'ctauTime4Delta')]
        sig = c.get('pr_width_sig') or [s * c['pr_width_k'] for s in (.010, .020, .040, .080)]
        for i, s in enumerate(sig):
            k[f'width_sigma{i + 1}'] = s
    k.update(np_constraint=int(c['np'] != 'free'), np_longest_free=int(c['np'] == 'long_free'), np_cov_scale=c['np_k'])
    k.update(profile_bkg=int(c['bkg'] == 'free'), background_penalty_scale=c['bkg_k'])
    if c.get('gap_k'):
        nres, q3, q4 = c['gaps']
        k.update(prompt_gap_ratio_center=c['gap_k'] * q3, prompt_gap_ratio_sigma=0.05 * c['gap_k'] * q3)
        if nres == 4: k.update(prompt_outer_gap_ratio_center=c['gap_k'] * q4, prompt_outer_gap_ratio_sigma=0.15 * c['gap_k'] * q4)
    k.update(c.get('extra', {}))  # diagnostic knobs only (e.g. profile_mass); not part of the simple-rule vocabulary
    return k


def write(rule):
    d = BASE / 'study' / 'rules' / rule
    d.mkdir(parents=True, exist_ok=True)
    for r, e in EDGES.items():
        for lo, hi in zip(e, e[1:]):
            y = 'y1.60_2.40' if r == 'fwd' else 'y0.00_1.60'
            c = choice(rule, r, lo)
            c['bin'] = f'{r}_{lo}_{hi}'
            c['purity'] = PURITY[f'{r}_{lo}_{hi}']
            c['gaps'] = GAPS[f'{r}_{lo}_{hi}']
            lines = [f'rule {rule}', '# ' + ' '.join(f'{a}={b}' for a, b in c.items())]
            lines += [f'{a} {b}' for a, b in knobs(c).items()]
            lines += [f'ct_low {window(r, lo)[0]}', f'ct_high {window(r, lo)[1]}']
            (d / f'{y}_pt{lo:.2f}_{hi:.2f}.txt').write_text('\n'.join(lines) + '\n')
    print('wrote', d)


def write_final_s1():
    write('FINAL_S1')
    src = BASE / 'study/rules/C23_gap1.8/y1.60_2.40_pt2.00_3.00.txt'
    text = src.read_text().replace('rule C23_gap1.8', 'rule FINAL_S1', 1)
    (BASE / 'study/rules/FINAL_S1/y1.60_2.40_pt2.00_3.00.txt').write_text(text)


if __name__ == '__main__':
    import fnmatch
    for pat in sys.argv[1:] or ['*']:
        for rule in fnmatch.filter(RULES, pat):
            write_final_s1() if rule == 'FINAL_S1' else write(rule)
