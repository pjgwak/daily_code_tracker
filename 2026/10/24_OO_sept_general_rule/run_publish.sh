#!/usr/bin/env bash
# Draw saved production results into figs_publish/; nothing is refitted.
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$script_dir/run_ctau_bins.sh" --publish "$@"
