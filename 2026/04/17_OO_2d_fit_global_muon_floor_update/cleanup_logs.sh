#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
target_dir="${script_dir}/logs"

if [[ ! -e "${target_dir}" ]]; then
  echo "[INFO] no logs directory: ${target_dir}"
  exit 0
fi

if [[ ! -d "${target_dir}" ]]; then
  echo "[ERROR] target exists but is not a directory: ${target_dir}" >&2
  exit 1
fi

case "${target_dir}" in
  "${script_dir}/logs")
    ;;
  *)
    echo "[ERROR] refusing to remove unexpected path: ${target_dir}" >&2
    exit 1
    ;;
esac

echo "[INFO] removing logs directory: ${target_dir}"
rm -rf -- "${target_dir}"
echo "[INFO] removed: ${target_dir}"
