#!/usr/bin/env python3
"""Rerun selected final fits from their frozen five-stage calibration inputs.

The fit model and bin settings come only from config/selection.json. Calibration
macros are in 02_fit/calibration; this command does not regenerate those prefits.
Existing final ROOT files are protected; use --output-base for a fresh rerun.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import uuid


def project_root(project=None):
    return Path(project).resolve() if project is not None else Path(__file__).resolve().parents[1]


def load_selection(project=None):
    return json.loads((project_root(project) / 'config/selection.json').read_text())


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _checked(path, digest=None):
    if not path.is_file():
        raise FileNotFoundError(path)
    if digest is not None and sha256(path) != digest:
        raise ValueError(f'Frozen input/source hash mismatch: {path}')
    return path


def _cpp_string(value):
    # JSON's quoted strings are also valid C++ literals for these filesystem paths.
    return json.dumps(str(value))


def run_bin(bin_id, *, project=None, output_base=None, dry_run=False):
    """Run one selected bin; return its run record (or a side-effect-free plan).

    output_base holds one directory per bin. The default is output/fit, whose
    seeded results cannot be overwritten. Use a fresh base for comparisons.
    """
    base = project_root(project)
    selection = load_selection(base)
    matches = [b for b in selection['bins'] if b['id'] == bin_id]
    if len(matches) != 1:
        raise ValueError(f'Unknown bin {bin_id!r}; available: ' + ', '.join(b['id'] for b in selection['bins']))
    item = matches[0]
    source = _checked(base / item['macro'])
    rule_file = _checked(base / item['rule_file'])
    prefit_base = base / 'output/prefits' / bin_id
    for info in item['prefits'].values():
        _checked(base / info['path'])
    data_info = next(v for v in selection['accepted_inputs'] if v['role'] == 'data')
    data = _checked(base / data_info['path'])
    style = _checked(base / selection['runtime_files']['style'])
    pr_windows = _checked(base / selection['runtime_files']['pr_calibration_windows'])
    destination_base = Path(output_base) if output_base is not None else Path('output/fit')
    if not destination_base.is_absolute():
        destination_base = base / destination_base
    destination_base = destination_base.resolve()
    destination = destination_base / bin_id
    relative_result = Path(item['output_root']).relative_to(Path('output/fit') / bin_id)
    result = destination / relative_result
    root = os.environ.get('ROOT_EXECUTABLE') or os.environ.get('ROOT_BIN') or shutil.which('root')
    fallback = Path.home() / 'build/root-6.36.06/bin/root'
    if root is None and fallback.is_file():
        root = str(fallback)
    if root is None:
        raise RuntimeError('ROOT is not on PATH; source thisroot.sh or set ROOT_EXECUTABLE.')
    root = str(Path(root).resolve())
    # Never inherit interactive overrides that could alter the selected fit.
    environment = {k: v for k, v in os.environ.items() if not k.startswith('OO_')}
    settings = {
        'OO_DATA_ROOT': str(data),
        'OO_ROOTLOGON': str(style),
        'OO_PR_CALIBRATION_WINDOWS': str(pr_windows),
        'OO_MASS_PREFIT_DIR': str(prefit_base),
        'OO_PR_PREFIT_DIR': str(prefit_base),
        'OO_NP_PREFIT_DIR': str(prefit_base),
        'OO_BKG_PREFIT_DIR': str(prefit_base),
        **item.get('environment', {}),
    }
    if item['rule_mode'] == 'knobs':
        settings['OO_RULE_DIR'] = str(base / item['rule_dir'])
    environment.update(settings)
    environment['ROOT_HIST'] = '0'
    numbers = ','.join(format(item[k], '.17g') for k in ('pt_low', 'pt_high', 'y_low', 'y_high'))
    call = f'{item["function"]}({numbers},false,false,false,{_cpp_string(destination)})'
    plan = {
        'bin': bin_id,
        'macro': str(source),
        'function_call': call,
        'output_root': str(result),
        'environment': settings,
        'prefits': {role: str(base / info['path']) for role, info in item['prefits'].items()},
        'existing_output': result.exists(),
        'fit_scope': selection['fit_scope'],
        'dry_run': bool(dry_run),
    }
    if dry_run:
        return plan
    if result.exists():
        raise FileExistsError(f'Refusing to overwrite {result}; choose a fresh --output-base.')
    # Validate the exact selected inputs and relocated source before launching ROOT.
    _checked(rule_file, item['rule_sha256'])
    for info in item['prefits'].values():
        _checked(base / info['path'], info['sha256'])
    _checked(data, data_info['sha256'])
    source_manifest = json.loads((base / 'config/source_manifest.json').read_text())
    for path, info in source_manifest['files'].items():
        _checked(base / path, info['sha256'])
    _checked(pr_windows, selection['runtime_files_sha256']['pr_calibration_windows'])
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '_' + uuid.uuid4().hex[:8]
    run_dir = base / 'output/runs' / bin_id / run_id
    work = run_dir / 'work'
    work.mkdir(parents=True)
    build = run_dir / 'build'
    build.mkdir()
    destination.mkdir(parents=True, exist_ok=True)
    # Reserve the destination while the process runs. Existing selected ROOTs above
    # are never touched, and two concurrent launches cannot write the same bin.
    lock = destination / '.fit.lock'
    try:
        lock_fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise RuntimeError(f'Fit destination is already locked: {lock}') from exc
    os.close(lock_fd)
    # A competing run may have finished during the input hash checks above.
    if result.exists():
        lock.unlink()
        raise FileExistsError(f'Refusing to overwrite {result}; choose a fresh --output-base.')
    if result.exists():
        lock.unlink()
        raise FileExistsError(f'Refusing to overwrite {result}; choose a fresh --output-base.')
    (work / 'roots').symlink_to(prefit_base / 'roots', target_is_directory=True)
    (work / 'inputs').symlink_to(base / 'output/inputs', target_is_directory=True)
    (work / 'include').symlink_to(base / 'include', target_is_directory=True)
    bootstrap = run_dir / 'run_fit.C'
    bootstrap.write_text(
        '#include <TROOT.h>\n#include <TSystem.h>\n#include <iostream>\n#include <exception>\n{\n'
        f'  gSystem->SetBuildDir({_cpp_string(build)}, true);\n'
        '  Int_t error = 0;\n'
        f'  gROOT->LoadMacro({_cpp_string(str(source) + "+")}, &error);\n'
        '  if (error) gSystem->Exit(10);\n'
        '  try {\n'
        f'    gROOT->ProcessLine({_cpp_string(call + ";")}, &error);\n'
        '  } catch (const std::exception &e) { std::cerr << e.what() << std::endl; gSystem->Exit(11); }\n'
        '  gSystem->Exit(error ? 12 : 0);\n}\n'
    )
    command = [root, '-n', '-l', '-b', '-q', str(bootstrap)]
    validation_call = f'{base / "02_fit/validate_saved_fit.C"}("fit2d",{numbers},0.001,2.0,0.0,{_cpp_string(destination)})'
    validation_command = [root, '-n', '-l', '-b', '-q', validation_call]
    record = {**plan, 'run_id': run_id, 'log': str(run_dir / 'fit.log'), 'validation_log': str(run_dir / 'validation.log'),
              'command': command, 'validation_command': validation_command, 'cwd': str(work),
              'source_sha256': sha256(source), 'rule_sha256': sha256(rule_file), 'state': 'running'}
    record_path = run_dir / 'run.json'
    record_path.write_text(json.dumps(record, indent=2) + '\n')
    try:
        with (run_dir / 'fit.log').open('w') as stream:
            completed = subprocess.run(command, cwd=work, env=environment, stdout=stream, stderr=subprocess.STDOUT, check=False)
        record['fit_returncode'] = completed.returncode
        if completed.returncode != 0:
            raise RuntimeError(f'ROOT fit failed ({completed.returncode}); see {run_dir / "fit.log"}')
        if not result.is_file():
            raise RuntimeError(f'ROOT produced no final file; see {run_dir / "fit.log"}')
        with (run_dir / 'validation.log').open('w') as stream:
            completed = subprocess.run(validation_command, cwd=work, env=environment, stdout=stream, stderr=subprocess.STDOUT, check=False)
        record['validation_returncode'] = completed.returncode
        if completed.returncode != 0:
            raise RuntimeError(f'Saved-fit validation failed; see {run_dir / "validation.log"}')
        record['output_sha256'] = sha256(result)
        record['state'] = 'complete'
    except BaseException as exc:
        record['state'] = 'failed'
        record['error'] = str(exc)
        raise
    finally:
        record['finished_utc'] = datetime.now(timezone.utc).isoformat()
        record_path.write_text(json.dumps(record, indent=2) + '\n')
        lock.unlink(missing_ok=True)
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--bin', action='append', dest='bins', help='selected bin id, repeatable')
    group.add_argument('--all', action='store_true', help='run all selected bins sequentially')
    parser.add_argument('--project', type=Path)
    parser.add_argument('--output-base', type=Path, help='fresh destination containing one directory per bin')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args(argv)
    ids = [b['id'] for b in load_selection(args.project)['bins']] if args.all else args.bins
    try:
        for bin_id in ids:
            print(json.dumps(run_bin(bin_id, project=args.project, output_base=args.output_base, dry_run=args.dry_run), indent=2), flush=True)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f'fit: {exc}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
