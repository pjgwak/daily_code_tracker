#!/usr/bin/env python3
"""Summarise rule-scan logs (logs/scan_*.log) against nominal ver3, per group.

  python3 study/v3_summary.py [--group fwd_1_5|fwd_5_10|fwd_10_20|mid|all] [--top N] [--rule PATTERN]

For every rule, over the bins of the group that have finished: max |%|, RMS %, number of bins within 2 %,
and the per-bin % (NP systematic = |%|; PR systematic = |%| x B/(1-B)).
"""
import argparse, fnmatch, glob, math, re, collections
PAT = re.compile(r'^(\S+)\s+(\S+)\s+B=([\d.]+)±([\d.]+)\s+([+-][\d.]+)%\s+([+-][\d.]+)σ\s+st(\d) cq(\d) mp ([\d.]+) tp ([\d.]+)')
GROUPS = {'fwd_1_3': ['fwd_1_2', 'fwd_2_3'], 'fwd_3_5': ['fwd_3_4', 'fwd_4_5'], 'fwd_5_10': ['fwd_5_6', 'fwd_6_7', 'fwd_7_8', 'fwd_8_9', 'fwd_9_10'],
          'fwd_10_20': ['fwd_10_12', 'fwd_12_14', 'fwd_14_20'], 'mid': ['mid_7_8', 'mid_8_9', 'mid_9_10', 'mid_10_12', 'mid_12_14', 'mid_14_16', 'mid_16_20', 'mid_20_30']}


def load():
    res = {}
    for f in glob.glob('logs/scan_*.log'):
        for line in open(f):
            m = PAT.match(line)
            if not m: continue
            rule, b, B, Be, pct, pull, st, cq, mp, tp = m.groups()
            res[(rule, b)] = dict(B=float(B), Be=float(Be), pct=float(pct), pull=float(pull), st=int(st), cq=int(cq), mp=float(mp), tp=float(tp))
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--group', default='all'); ap.add_argument('--top', type=int, default=15); ap.add_argument('--rule', default='*')
    ap.add_argument('--bin', default='', help='rank rules for one bin instead of a group')
    a = ap.parse_args(); res = load()
    if a.bin:
        rows = sorted(((r, v) for (r, b), v in res.items() if b == a.bin and fnmatch.fnmatch(r, a.rule)), key=lambda x: abs(x[1]['pct']))
        print(f'## {a.bin}: {len(rows)} fits, within 2%: {sum(abs(v["pct"]) <= 2 and v["st"] == 0 and v["cq"] == 3 for _, v in rows)}')
        for r, v in rows[:a.top]:
            print(f'  {r:<40} B={v["B"]:.4f}±{v["Be"]:.4f} {v["pct"]:+6.2f}% {v["pull"]:+5.2f}σ st{v["st"]} cq{v["cq"]} mp {v["mp"]:.2f} tp {v["tp"]:.2f}')
        return
    for g, bins in GROUPS.items():
        if a.group not in ('all', g): continue
        rules = sorted({r for (r, b) in res if b in bins and fnmatch.fnmatch(r, a.rule)})
        table = []
        for r in rules:
            vs = [(b, res.get((r, b))) for b in bins]
            done = [(b, v) for b, v in vs if v]
            if not done: continue
            bad = [b for b, v in done if v['st'] != 0 or v['cq'] != 3]
            pcts = [v['pct'] for b, v in done]
            table.append((r, len(done), len(bins), max(abs(p) for p in pcts), math.sqrt(sum(p * p for p in pcts) / len(pcts)), sum(abs(p) <= 2 for p in pcts), bad, {b: v['pct'] for b, v in done}))
        table.sort(key=lambda t: (-(t[1] == t[2]), t[3]))
        print(f'\n## {g}  ({len(table)} rules; sorted by max|%| over finished bins; complete rules first)')
        print(f'{"rule":<40} {"done":>5} {"max|%|":>7} {"rms%":>6} {"≤2%":>4}  ' + ' '.join(f'{b.split("_",1)[1]:>7}' for b in bins))
        for r, n, nb, mx, rms, ok, bad, d in table[:a.top]:
            print(f'{r:<40} {n:>2}/{nb:<2} {mx:7.2f} {rms:6.2f} {ok:>4}  ' + ' '.join(f'{d[b]:+7.2f}' if b in d else f'{"-":>7}' for b in bins) + (f'  BAD:{",".join(bad)}' if bad else ''))


if __name__ == '__main__':
    main()
