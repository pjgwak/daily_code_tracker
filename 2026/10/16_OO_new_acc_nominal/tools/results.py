#!/usr/bin/env python3
"""Build the agreement CSV/plot from current fit ROOTs and frozen nominal values.

Examples (from any working directory):
  python tools/results.py --root /path/to/root
  python tools/results.py --allow-new-fits

The default requires the selected ROOT hashes. --allow-new-fits accepts rerun
ROOTs while preserving ROOT fit-quality and no-common-scale checks. Nominal
values always come from config/selection.json; no archive is consulted.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

PROJECT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def project_path(value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        raise ValueError(f"Config paths must be project-relative: {value}")
    resolved = (PROJECT / path).resolve()
    if not resolved.is_relative_to(PROJECT):
        raise ValueError(f"Config path escapes the project: {value}")
    return resolved


def output_path(value: str) -> Path:
    path = Path(value)
    return (PROJECT / path).resolve() if not path.is_absolute() else path.resolve()


def root_executable(value: str | None) -> str:
    selected = value or os.environ.get("ROOT_EXECUTABLE") or os.environ.get("ROOT_BIN") or "root"
    executable = shutil.which(selected)
    fallback = Path.home() / "build/root-6.36.06/bin/root"
    if not executable and selected == "root" and fallback.is_file():
        executable = str(fallback)
    if not executable:
        raise ValueError("ROOT is unavailable. Put root on PATH, set ROOT_EXECUTABLE, or pass --root /path/to/root.")
    return executable


def load_selection(path: Path, allow_new: bool) -> tuple[dict, list[dict]]:
    selection = json.loads(path.read_text())
    bins = selection.get("bins", [])
    if selection.get("schema_version") != 1 or len(bins) != 20:
        raise ValueError("Expected selection schema version 1 with exactly 20 bins")
    if selection.get("threshold_sigma") != 0.414:
        raise ValueError("The preserved agreement plot requires threshold_sigma = 0.414")
    identities = set()
    counts = {"fwd": 0, "mid": 0}
    records = []
    for item in bins:
        identity = (item["region"], float(item["pt_low"]), float(item["pt_high"]))
        if identity in identities or identity[0] not in counts or identity[2] <= identity[1]:
            raise ValueError(f"Duplicate or invalid bin: {item['id']}")
        identities.add(identity)
        counts[identity[0]] += 1
        expected = item["expected"]
        for key in ("B_raw", "B_raw_error", "nominal", "nominal_stat_error"):
            if not math.isfinite(float(expected[key])):
                raise ValueError(f"Nonfinite {key} in {item['id']}")
        if expected["nominal"] <= 0 or expected["nominal_stat_error"] <= 0:
            raise ValueError(f"Invalid nominal reference in {item['id']}")
        source = project_path(item["output_root"])
        if not source.is_file():
            raise ValueError(f"Missing current ROOT for {item['id']}: {source}")
        actual_hash = sha256(source)
        matches = actual_hash == expected["root_sha256"]
        if not matches and not allow_new:
            raise ValueError(f"ROOT hash changed for {item['id']}. Use --allow-new-fits to collect intentional reruns.")
        records.append({"bin": item, "source": source, "sha256": actual_hash, "frozen_match": matches})
    if counts != {"fwd": 12, "mid": 8} or len({r["bin"]["id"] for r in records}) != 20:
        raise ValueError("Expected 12 forward and 8 mid bins with unique ids")
    records.sort(key=lambda r: (r["bin"]["region"] != "fwd", r["bin"]["pt_low"]))
    return selection, records


def figure_inventory(records: list[dict], destination: Path) -> list[tuple[Path, Path]]:
    copies = []
    for record in records:
        item = record["bin"]
        source = project_path(item.get("figure_dir", f"output/fit/{item['id']}/figs"))
        files = sorted(p for p in source.rglob("*") if p.is_file() and p.suffix.lower() in {".pdf", ".png"})
        names = [p.name for p in files]
        if len(names) != len(set(names)):
            raise ValueError(f"Duplicate figure names under {source}")
        if not any(n.startswith("mass_fit_") and n.endswith(".pdf") for n in names) or not any(n.startswith("lifetime_fit_") and n.endswith(".pdf") for n in names):
            raise ValueError(f"Missing mass/lifetime PDFs for {item['id']}: {source}; use --skip-figures for plot-only checks")
        for path in files:
            target = destination / item["id"] / path.name
            if target.resolve() == path.resolve():
                raise ValueError(f"Figure output would overwrite its fit input: {path}")
            copies.append((path, target))
    return copies


def collect_table(path: Path, records: list[dict]) -> list[dict]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != len(records):
        raise ValueError("ROOT did not produce all 20 agreement rows")
    lookup = {(r["bin"]["region"], float(r["bin"]["pt_low"]), float(r["bin"]["pt_high"])): r for r in records}
    seen = set()
    output = []
    for row in rows:
        key = (row["region"], float(row["pt_low"]), float(row["pt_high"]))
        if key not in lookup or key in seen:
            raise ValueError("Unexpected or repeated agreement row from ROOT")
        seen.add(key)
        record = lookup[key]
        item, expected = record["bin"], record["bin"]["expected"]
        if Path(row["source_root"]).resolve() != record["source"]:
            raise ValueError(f"ROOT input mismatch for {item['id']}")
        for column, reference in (("B_nominal", "nominal"), ("nominal_stat_error", "nominal_stat_error")):
            if float(row[column]) != expected[reference]:
                raise ValueError(f"Nominal reference changed for {item['id']}")
        if record["frozen_match"]:
            for column, reference in (("B_raw", "B_raw"), ("B_raw_fit_error", "B_raw_error")):
                if not math.isclose(float(row[column]), expected[reference], rel_tol=0, abs_tol=1e-12):
                    raise ValueError(f"Frozen {column} mismatch for {item['id']}")
        row["source_root"] = item["output_root"]
        output.append({"id": item["id"], "rule": item["rule"], **row, "source_root_sha256": record["sha256"]})
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/selection.json")
    parser.add_argument("--output-dir", default="results", help="Agreement CSV, PDF, PNG and summary destination")
    parser.add_argument("--figures-dir", default="figs", help="Per-bin figure destination")
    parser.add_argument("--work-dir", default="output/results", help="Transient build files and ROOT log directory")
    parser.add_argument("--root", help="ROOT executable (otherwise ROOT_EXECUTABLE or PATH)")
    parser.add_argument("--allow-new-fits", action="store_true", help="Accept changed current ROOT hashes after an intentional rerun")
    parser.add_argument("--skip-figures", action="store_true", help="Generate agreement outputs without copying fit figures")
    args = parser.parse_args()
    try:
        executable = root_executable(args.root)
        renderer = shutil.which("pdftoppm")
        if not renderer:
            raise ValueError("pdftoppm is unavailable; install Poppler utilities to render the agreement PNG.")
        config = output_path(args.config)
        selection, records = load_selection(config, args.allow_new_fits)
        destination = output_path(args.output_dir)
        figure_destination = output_path(args.figures_dir)
        work = output_path(args.work_dir)
        copies = [] if args.skip_figures else figure_inventory(records, figure_destination)
        work.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="agreement-", dir=work) as temporary:
            staging = Path(temporary)
            table = staging / "inputs.tsv"
            with table.open("w") as stream:
                for record in records:
                    item, expected = record["bin"], record["bin"]["expected"]
                    values = [int(item["region"] == "mid"), item["pt_low"], item["pt_high"], expected["nominal"], expected["nominal_stat_error"]]
                    stream.write("\t".join(map(str, values)) + "\t" + json.dumps(str(record["source"])) + "\n")
            macro = PROJECT / "03_results/draw_agreement.C"
            expression = f"{macro}({json.dumps(str(table))},{json.dumps(str(staging))})"
            completed = subprocess.run([executable, "-n", "-l", "-b", "-q", expression], cwd=PROJECT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            log = work / "results.log"
            log.write_text(completed.stdout)
            if completed.returncode:
                raise ValueError(f"ROOT plotting failed (exit {completed.returncode}); see {log}")
            rendered = subprocess.run([renderer, "-singlefile", "-png", "-r", "150", str(staging / "nominal_agreement.pdf"), str(staging / "nominal_agreement")], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if rendered.returncode:
                with log.open("a") as stream:
                    stream.write(rendered.stdout)
                raise ValueError(f"PNG rendering failed; see {log}")
            for name in ("agreement.csv", "nominal_agreement.pdf", "nominal_agreement.png", "agreement_summary.json"):
                if not (staging / name).is_file() or not (staging / name).stat().st_size:
                    raise ValueError(f"ROOT did not produce {name}; see {log}")
            rows = collect_table(staging / "agreement.csv", records)
            with (staging / "agreement.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            summary = json.loads((staging / "agreement_summary.json").read_text())
            summary.update(selection=selection["selection"], bins=len(rows), count_within_sigma0414=sum(abs(float(row["signed_nominal_stat_units"])) <= .414 for row in rows), frozen_root_hashes_matched=sum(r["frozen_match"] for r in records), figures_copied=len(copies))
            (staging / "agreement_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
            destination.mkdir(parents=True, exist_ok=True)
            for source, target in copies:
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
            for name in ("agreement.csv", "nominal_agreement.pdf", "nominal_agreement.png", "agreement_summary.json"):
                shutil.copy2(staging / name, destination / name)
        print(f"Collected {len(rows)} current ROOTs; {summary['count_within_sigma0414']}/20 within 0.414 nominal-stat-error units.")
        print(f"Agreement outputs: {destination}")
        if copies:
            print(f"Copied {len(copies)} per-bin figures: {figure_destination}")
        return 0
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"results: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
