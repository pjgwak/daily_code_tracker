#!/usr/bin/env python3
import os
import subprocess
import sys


# Usage:
#   python3 run_script.py [nEvt] fill -1
#   python3 run_script.py ratio
#   python3 run_script.py fit
#   python3 run_script.py [nEvt] all
# Example:
#   python3 run_script.py -1 fill
#   python3 run_script.py 1000 fill ratio

def main(argv):
    n_evt = 1000
    step_tokens = []
    for token in argv:
        try:
            n_evt = int(token)
            continue
        except ValueError:
            step_tokens.append(token)

    if not step_tokens:
        print("No step specified. Nothing to run.")
        return 0

    selected = set()
    for token in step_tokens:
        key = token.lower()
        if key == "all":
            selected.update(["fill", "ratio", "fit"])
            continue
        if key in ("fill", "ratio", "fit"):
            selected.add(key)
            continue
        print(f"Unknown step: {token}", file=sys.stderr)
        return 2

    cwd = os.path.dirname(os.path.abspath(__file__))
    commands = []
    if "fill" in selected:
        if n_evt == 1000:
            print("Only handle 1000 events. Please add '-1' if you want to loop to all events")
        commands.extend(
            [
                ["root", "-l", "-b", "-q", f"make_pt_bins.C({n_evt},0,1)"],  # fil bins for data
                ["root", "-l", "-b", "-q", f"make_pt_bins.C({n_evt},1,1)"],  # fil bins for mc pr
                ["root", "-l", "-b", "-q", f"make_pt_bins.C({n_evt},1,0)"],  # fil bins for mc np
            ]
        )
    if "ratio" in selected:
        commands.append(["root", "-l", "-b", "-q", "make_ratio_bins.C"])  # ratio
    if "fit" in selected:
        commands.append(["root", "-l", "-b", "-q", "fit_ratio.C"])  # fit

    for cmd in commands:
        print("+", " ".join(cmd))
        subprocess.Popen(cmd, cwd=cwd)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
