#!/usr/bin/env python3
"""Run extra final-fit rules on the prefits of an existing prefit-variant work directory.

  python3 study/refit_rules.py WORKDIR BIN RULE [RULE ...]     -> WORKDIR/rules/<RULE>/roots/.../fit2d, one summary line each
The final-fit ctau window follows WORKDIR/config/prefit_settings.txt (ct_low/ct_high) as in scan_prefit_rules.py.
"""
import sys
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / 'study'))
import rules, run_rules as rr, scan_prefit_rules as sp  # noqa: E402

def main():
    work, b, names = Path(sys.argv[1]).resolve(), sys.argv[2], sys.argv[3:]
    rr.nominal()
    lib = next((BASE / 'build').rglob('fit2d_C.so'))
    for n in names:
        if n in rules.RULES and not (BASE / 'study/rules' / n).exists(): rules.write(n)
        r = sp.fit(work, b, n, lib)
        print(f'{work.name:<30} {n:<32} ' + (f'B={r["B"]:.4f}±{r["Berr"]:.4f} {r["pct"]:+6.2f}% {r["pull"]:+5.2f}σ st{r["status"]} cq{r["covQual"]} mp {r["massP"]:.2f} tp {r["timeP"]:.2f} {"PASS" if r["pass_"] else "fail"}' if 'pct' in r else f'FAIL {r.get("fail")}'), flush=True)

if __name__ == '__main__':
    main()
