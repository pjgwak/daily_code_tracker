#!/usr/bin/env python3
"""fwd 10-12 final-fit-only knob grid (2026-10-28 night) on the production prefits: rules/G10_<pw>_<mean>_<np>_<bkg>_<lo>/ (fwd 10-12 file only).
pw: tail k (core fixed, tails x k) or fixed (prompt_width_fixed, needs build_fixw); np: long (longest lifetime free), all (all free), mc;
bkg: free (profile_bkg 1) or sb<k> (constrained to the sideband prefit, penalty scale k)."""
from pathlib import Path
import itertools
H = Path(__file__).resolve().parents[1]
base = (H / 'rules/FINAL_V2a/y1.60_2.40_pt10.00_12.00.txt').read_text().splitlines()
def kv(lines):
    return [l for l in lines if l and not l.startswith(('rule', '#'))]
names = []
for pw, mean, np_, bkg, lo in itertools.product(['t0.1', 't0.25', 't1', 'fix'], ['mf', 'mx'], ['long', 'all', 'mc'], ['free', 'sb1', 'sb0.25'], ['-0.5', '-0.3', '-0.15']):
    d = {l.split()[0]: l.split()[1] for l in kv(base)}
    if pw == 'fix': d['prompt_width_fixed'] = '1'
    else:
        k = float(pw[1:]); d.update(width_policy='1', width_sigma1=f'{0.01*k:g}', width_sigma2=f'{0.02*k:g}', width_sigma3=f'{0.04*k:g}', width_sigma4=f'{0.08*k:g}')
    d['float_prompt_mean'] = '1' if mean == 'mf' else '0'
    d['np_longest_free'] = '1' if np_ == 'long' else '0'
    if np_ == 'all': d['np_all_lifetimes_free'] = '1'
    if bkg == 'free': d['profile_bkg'] = '1'
    else: d.update(profile_bkg='0', background_penalty_scale=bkg[2:])
    d['ct_low'] = lo
    n = f'G10_{pw}_{mean}_{np_}_{bkg}_{lo}'
    (H / 'rules' / n).mkdir(parents=True, exist_ok=True)
    (H / 'rules' / n / 'y1.60_2.40_pt10.00_12.00.txt').write_text(f'rule {n}\n' + ''.join(f'{k} {v}\n' for k, v in d.items()))
    names.append(n)
print(' '.join(names))
