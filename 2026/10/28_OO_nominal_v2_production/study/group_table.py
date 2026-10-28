#!/usr/bin/env python3
"""Per-rule pass count over a bin group from study/runs logs (lines printed by run_rules.py).
  python3 study/group_table.py fwd_5_6,fwd_6_7,... [log ...]"""
import re, sys, glob, collections
bins = sys.argv[1].split(','); logs = sys.argv[2:] or glob.glob('study/runs/batch*.log')
res = collections.defaultdict(dict)
pat = re.compile(r'^(\S+)\s+(\S+)\s+B=([\d.]+)±([\d.]+)\s+([+-][\d.]+)%\s+([+-][\d.]+)σ\s+st(-?\d+) cq(-?\d+).*?(PASS|fail)')
for l in logs:
    for line in open(l):
        m = pat.match(line)
        if m and m[2] in bins: res[m[1]][m[2]] = (float(m[5]), m[9] == 'PASS', float(m[4]))
rows = sorted(res.items(), key=lambda kv: (-sum(v[1] for v in kv[1].values()), max(abs(v[0]) for v in kv[1].values())))
print(f'{"rule":<34} pass  ' + ' '.join(f'{b:>10}' for b in bins))
for r, d in rows[:25]:
    print(f'{r:<34} {sum(v[1] for v in d.values())}/{len(d)}  ' + ' '.join(f'{d[b][0]:+9.2f}{"*" if d[b][1] else " "}' if b in d else f'{"-":>10}' for b in bins))
