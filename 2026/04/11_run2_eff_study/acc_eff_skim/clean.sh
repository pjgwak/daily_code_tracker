#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

patterns=(
  "*_C.so"
  "*_C.d"
  "*_C_ACLiC_dict_rdict.pcm"
)

found_any=0

for pattern in "${patterns[@]}"; do
  while IFS= read -r -d '' file; do
    found_any=1
    rm -f "$file"
    printf 'removed %s\n' "${file#$script_dir/}"
  done < <(find "$script_dir" -maxdepth 1 -type f -name "$pattern" -print0)
done

if [[ "$found_any" -eq 0 ]]; then
  echo "No ROOT build artifacts found."
fi
