#!/usr/bin/env python3
"""Regenerate accepted inputs, keeping partial runs separate from production."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


BASE = Path(__file__).resolve().parents[1]
PRODUCTION = BASE / "output/inputs"
OUTPUT_NAMES = {
    "data": "roodataset_roots/RooDataSet_miniAOD_isMC0_globalOff_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root",
    "prompt": "roodataset_roots/RooDataSet_miniAOD_isMC1_PR_Hijing_globalOff_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root",
    "nonprompt": "np_gen_singlemu.root",
}


def root_executable():
    requested = os.environ.get("ROOT_EXECUTABLE") or os.environ.get("ROOT_BIN")
    if requested:
        root = shutil.which(requested)
        if not root:
            raise RuntimeError(f"ROOT_BIN is not executable: {requested}")
        return root
    root = shutil.which("root")
    if root:
        return root
    fallback = Path.home() / "build/root-6.36.06/bin/root"
    if fallback.is_file() and os.access(fallback, os.X_OK):
        return str(fallback)
    raise RuntimeError("ROOT is unavailable; set ROOT_BIN or add root to PATH.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample", choices=tuple(OUTPUT_NAMES))
    parser.add_argument("--max-events", type=int, default=-1,
                        help="-1 processes all events (default); positive values require a separate --output-dir")
    parser.add_argument("--output-dir", type=Path,
                        help="output root, relative to the current working directory (default: repository output/inputs)")
    parser.add_argument("--config", type=Path, default=BASE / "config/input_sources.json")
    parser.add_argument("--dry-run", action="store_true", help="print planned paths and ROOT arguments without writing")
    args = parser.parse_args()
    if args.max_events == 0 or args.max_events < -1:
        parser.error("--max-events must be -1 or a positive integer")
    output_dir = (args.output_dir or PRODUCTION).expanduser().resolve()
    production = PRODUCTION.resolve()
    if args.max_events != -1 and (
            args.output_dir is None or output_dir == production or production in output_dir.parents):
        parser.error("partial runs require --output-dir outside production output/inputs")
    config = json.loads(args.config.expanduser().read_text())
    source = Path(config[args.sample]).expanduser()
    if not source.is_absolute():
        source = args.config.expanduser().resolve().parent / source
    source = source.resolve()
    output = output_dir / OUTPUT_NAMES[args.sample]
    build = BASE / "output/build/skim" / args.sample
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    log = BASE / "output/logs/skim" / f"{args.sample}_{timestamp}.log"
    # json.dumps provides C++ string literals; subprocess receives an argv list,
    # so file names never become shell code.
    quote = lambda value: json.dumps(str(value))
    if args.sample == "nonprompt":
        macro = BASE / "01_skim/make_np_hijing_gen.C"
        arguments = f"{quote(source)},{quote(output)},{args.max_events}"
    else:
        macro = BASE / "01_skim/skim_to_roodataset.C"
        mc = "false,0" if args.sample == "data" else "true,1"
        arguments = (f'{args.max_events},{mc},false,"_OO25_HLT_OxyL1SingleMuOpen_v1",'
                     f"false,false,false,false,0,200,2.6,3.5,true,0,{quote(source)},{quote(output)}")
    command = [root_executable(), "-n", "-l", "-b", "-e",
               f"gSystem->SetBuildDir({quote(build)}, true);", "-q", f"{macro}+({arguments})"]
    plan = {"sample": args.sample, "max_events": args.max_events,
            "input": str(source), "input_exists": source.is_file(),
            "output": str(output), "output_exists": output.exists(),
            "cwd": str(BASE), "log": str(log), "argv": command}
    if args.dry_run:
        print(json.dumps(plan, indent=2))
        return 0
    if not source.is_file():
        raise RuntimeError(f"Input is unavailable: {source}; edit {args.config} or use --config.")
    if output.exists():
        raise RuntimeError(f"Refusing to overwrite {output}; choose a new --output-dir.")
    for directory in (output.parent, build, log.parent):
        directory.mkdir(parents=True, exist_ok=True)
    print(f"Skimming {args.sample}; log: {log}", flush=True)
    with log.open("x") as stream:
        stream.write(json.dumps(plan, indent=2) + "\n")
        stream.flush()
        process = subprocess.run(command, cwd=BASE, stdout=stream, stderr=subprocess.STDOUT)
    if process.returncode != 0 or not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"Skim failed (ROOT exit {process.returncode}); inspect {log}.")
    print(f"Created {output}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, KeyError, ValueError, RuntimeError) as error:
        print(f"skim: {error}", file=sys.stderr)
        sys.exit(1)
