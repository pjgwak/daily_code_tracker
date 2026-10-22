#!/usr/bin/env python3
"""Check the published selection against preserved source and artifact hashes."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

BASE = Path(__file__).resolve().parents[1]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify(project=BASE, code_only=False):
    project = Path(project).resolve()
    config = json.loads((project / "config/selection.json").read_text())
    sources = json.loads((project / "config/source_manifest.json").read_text())
    failures = []
    counts = {"source": 0, "runtime": 0, "rule": 0, "input": 0, "prefit": 0, "final": 0}

    def check(relative, expected, category):
        path = project / relative
        if Path(relative).is_absolute() or not path.resolve().is_relative_to(project):
            failures.append(f"Path escapes project: {relative}")
        elif not path.is_file():
            failures.append(f"Missing {category}: {relative}")
        elif sha256(path) != expected:
            failures.append(f"Changed {category}: {relative}")
        else:
            counts[category] += 1

    bins = config["bins"]
    if len(bins) != 20 or len({row["id"] for row in bins}) != 20:
        failures.append("Expected exactly 20 distinct selected bins")
    for relative, record in sources["files"].items():
        check(relative, record["sha256"], "source")
    for name, digest in config.get("runtime_files_sha256", {}).items():
        check(config["runtime_files"][name], digest, "runtime")
    for row in bins:
        check(row["rule_file"], row["rule_sha256"], "rule")
        expected = row["expected"]
        distance = (expected["B_raw"] - expected["nominal"]) / expected["nominal_stat_error"]
        if abs(distance - expected["sigma_signed"]) > 1e-12 or abs(distance) > config["threshold_sigma"]:
            failures.append(f'Invalid frozen agreement reference: {row["id"]}')
        if not code_only:
            check(row["output_root"], expected["root_sha256"], "final")
            if len(row["prefits"]) != 5:
                failures.append(f'Expected five consumed prefits: {row["id"]}')
            for prefit in row["prefits"].values():
                check(prefit["path"], prefit["sha256"], "prefit")
    if not code_only:
        for record in config["accepted_inputs"]:
            check(record["path"], record["sha256"], "input")
    return {
        "selection": config["selection"],
        "scope": "code_only" if code_only else "preserved_selection_integrity",
        "passed": not failures,
        "verified": counts,
        "failures": failures,
        "interpretation": "Hash verification of the selected artifacts; not an additional fit, model-stability or coverage test.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--code-only", action="store_true", help="check tracked code and rules without local ROOT datasets")
    parser.add_argument("--report", type=Path, help="write the verification report to this path")
    args = parser.parse_args()
    report = verify(code_only=args.code_only)
    if args.report:
        path = args.report if args.report.is_absolute() else BASE / args.report
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, KeyError, ValueError) as error:
        print(f"Verification failed: {error}", file=sys.stderr)
        raise SystemExit(1)
