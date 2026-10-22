#!/usr/bin/env python3
"""OO analysis entry point; resolve all paths relative to this file."""
from pathlib import Path
import argparse
import json
import subprocess
import sys

BASE = Path(__file__).resolve().parent


def selection():
    return json.loads((BASE / "config/selection.json").read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("bins", "status", "results", "verify", "run-bin", "skim"))
    parser.add_argument("arguments", nargs=argparse.REMAINDER, help="command options; use COMMAND --help")
    args = parser.parse_args()
    if args.command in ("bins", "status"):
        if args.arguments:
            parser.error("unrecognized arguments: " + " ".join(args.arguments))
        bins = selection()["bins"]
        if args.command == "bins":
            for row in bins:
                lo, hi = row["ctau_window"]
                print(f'{row["id"]:12} {row["rule"]:30} ctau [{lo:g}, {hi:g}] mm')
        else:
            complete = 0
            for row in bins:
                final = (BASE / row["output_root"]).is_file()
                prefits = sum((BASE / p["path"]).is_file() for p in row["prefits"].values())
                complete += final and prefits == 5
                print(f'{row["id"]:12} final={"present" if final else "missing":7} prefits={prefits}/5')
            print(f"{complete}/{len(bins)} bins have saved inputs and final outputs; run verify for integrity.")
        return 0
    script = {"run-bin": "fit", "results": "results", "verify": "verify", "skim": "skim"}[args.command]
    if args.command == "run-bin" and args.arguments and not args.arguments[0].startswith("-"):
        args.arguments.insert(0, "--bin")
    return subprocess.run(
        [sys.executable, str(BASE / "tools" / f"{script}.py"), *args.arguments],
        cwd=BASE,
    ).returncode


if __name__ == "__main__":
    raise SystemExit(main())
