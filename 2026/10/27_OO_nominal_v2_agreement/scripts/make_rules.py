#!/usr/bin/env python3
"""Write knob dirs rules/<name>/ for fwd 10-20 from the 25/study/rules.py vocabulary (rules.choice/knobs), name -> overrides.
  python3 scripts/make_rules.py   (definitions below; existing dirs are overwritten)"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE.parent / '25_OO_sept_acc_general_rule/study'))
import rules  # noqa: E402

DEFS = {}
# Grid AT (2026-10-28): NP lifetimes all free, NP fractions constrained at MC (x np_k), with PR treatments; NP2 prefits of FINAL_S3.
for prn, pr in (('tail0.1', dict(pr_width='tail', pr_width_k=0.1)), ('abs0.25', dict(pr_width_k=0.25)), ('abs1', dict(pr_width_k=1)),
                ('tail0.25', dict(pr_width='tail', pr_width_k=0.25)), ('j64', dict(pr_width='mc_joint', pr_joint_k=64))):
    for nk in (1, 2, 16):
        for bkg in ('free', 'sideband'):
            DEFS[f'AT_{prn}_np{nk}_{bkg}'] = dict(pr, np='mc', np_k=nk, bkg=bkg, extra=dict(np_all_lifetimes_free=1))

if __name__ == '__main__':
    for name, ov in DEFS.items():
        d = HERE / 'rules' / name; d.mkdir(parents=True, exist_ok=True)
        for lo, hi in ((10, 12), (12, 14), (14, 20)):
            b = f'fwd_{lo}_{hi}'
            c = dict(rules.DEFAULT); c.update(ov); c.update(bin=b, purity=rules.PURITY[b], gaps=rules.GAPS[b])
            lines = [f'rule {name}', '# ' + ' '.join(f'{a}={v}' for a, v in c.items())] + [f'{a} {v}' for a, v in rules.knobs(c).items()]
            lines += [f'ct_low {rules.window("fwd", lo)[0]}', f'ct_high {rules.window("fwd", lo)[1]}']
            (d / f'y1.60_2.40_pt{lo:.2f}_{hi:.2f}.txt').write_text('\n'.join(lines) + '\n')
    print(len(DEFS), 'rules')
