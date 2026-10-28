#!/usr/bin/env python3
"""Run extra final-fit rules on the prefits of a prefit_variants/ or scans/ work dir (scored vs nominal ver2).
  python3 scripts/refit_v2.py WORKDIR BIN RULE [RULE ...]"""
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, str(Path(__file__).resolve().parent))
import scan_v2 as sp  # noqa: E402
sp.rr.nominal()
work, b, names = Path(sys.argv[1]).resolve(), sys.argv[2], sys.argv[3:]
lib = next((sp.BASE / 'build').rglob('fit2d_C.so'))
def one(n):
    r = sp.fit(work, b, n, lib)
    print(f'{work.name:<26} {n:<28} ' + (f'B={r["B"]:.4f}±{r["Berr"]:.4f} {r["pct"]:+6.2f}% {r["pull"]:+5.2f}σ st{r["status"]} cq{r["covQual"]} mp {r["massP"]:.2f} tp {r["timeP"]:.2f}' if 'pct' in r else f'FAIL {r.get("fail")}'), flush=True)
with ThreadPoolExecutor(len(names)) as pool: list(pool.map(one, names))
