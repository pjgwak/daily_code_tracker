#!/usr/bin/env python3
"""Replay the saved final 2D fits (FINAL_S3, no minimization) and require projection diagnostics.

Only fit2d.C writes diagnostics in this folder; the prefit macros have no diagnostics hook.
The replay goes through run_publish.sh so each bin gets its production settings (OO_RULE_DIR etc.).
"""
import argparse, hashlib, json, os, subprocess, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument('--jobs', type=int, default=4)
p.add_argument('--output', default='drawing/projection_diagnostics')
a = p.parse_args()
if not 1 <= a.jobs <= 12: p.error("--jobs must be between 1 and 12")
def hashes():
    return {str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted((ROOT/"roots").rglob("*.root"))}
before = hashes()
out = (ROOT/a.output).resolve()
out.mkdir(parents=True, exist_ok=True)
bins = [(0,1.6,l,h) for l,h in zip([7,8,9,10,12,14,16,20],[8,9,10,12,14,16,20,30])]
bins += [(1.6,2.4,l,h) for l,h in zip([1,2,3,4,5,6,7,8,9,10,12,14],[2,3,4,5,6,7,8,9,10,12,14,20])]
start = time.time()
with (out/'replay.log').open('w') as log:
    rc = subprocess.run(['bash','run_publish.sh','--2d',f'--max-jobs={a.jobs}'],cwd=ROOT,
                        env=dict(os.environ,OO_DIAGNOSTICS_DIR=str(out)),stdout=log,stderr=subprocess.STDOUT).returncode
results = []
for yl,yh,lo,hi in bins:
    tag = f'y{yl:.2f}_{yh:.2f}_pt{lo:.2f}_{hi:.2f}'
    names = [f'fit2d_{v}_{tag}' for v in ['mass','ctau']]
    ok = all((out/f'{n}.{e}').is_file() and (out/f'{n}.{e}').stat().st_mtime >= start for n in names for e in ('csv','pdf'))
    results.append(dict(stage='fit2d',tag=tag,ok=ok))
(out/'run_results.json').write_text(json.dumps(dict(returncode=rc,bins=results),indent=2)+'\n')
after = hashes()
(out/'input_integrity.json').write_text(json.dumps(dict(unchanged=before==after, before=before, after=after),indent=2)+'\n')
if before != after: raise RuntimeError('Saved ROOT files changed during replay')
print(f'returncode={rc}, bins ok {sum(r["ok"] for r in results)}/{len(results)}')
raise SystemExit(0 if rc == 0 and all(r['ok'] for r in results) else 1)
