#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}"

shopt -s nullglob

compiled_files=(
  *_C.so
  *_C.d
  *_C_ACLiC_dict_rdict.pcm
  *_C_ACLiC_dict.cxx
  *_C_ACLiC_dict.h
  *_C_ACLiC_dict.o
  *_C_ACLiC_dict.rootmap
  *_C_ACLiC_linkdef.h
)

if (( ${#compiled_files[@]} == 0 )); then
  echo "[INFO] no compiled ROOT ACLiC files found in: ${script_dir}"
  exit 0
fi

echo "[INFO] removing compiled ROOT ACLiC files in: ${script_dir}"
printf '  %s\n' "${compiled_files[@]}"
rm -f -- "${compiled_files[@]}"
echo "[INFO] removed ${#compiled_files[@]} file(s)"
