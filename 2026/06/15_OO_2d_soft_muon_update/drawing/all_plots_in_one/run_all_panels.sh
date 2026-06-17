#!/usr/bin/env bash
#######################
# Need LaTex compiler to run this script.
#######################
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
FIGS_DIR="${FIGS_DIR:-${SCRIPT_DIR}/../../figs_publish}"

for dir in "${FIGS_DIR}"/*/*; do
  [ -d "${dir}" ] || continue

  count="$(find "${dir}" -maxdepth 1 -type f -name '*.pdf' | wc -l | tr -d ' ')"
  [ "${count}" -gt 0 ] || continue

  rapidity="$(basename "$(dirname "${dir}")")"
  kind="$(basename "${dir}")"

  if [ "${kind}" = "fit2d" ]; then
    for pattern in mass_fit*.pdf lifetime_fit*.pdf; do
      split_count="$(find "${dir}" -maxdepth 1 -type f -name "${pattern}" | wc -l | tr -d ' ')"
      [ "${split_count}" -gt 0 ] || continue

      split_kind="${pattern%%_fit*.pdf}"
      output="${SCRIPT_DIR}/${rapidity}_${kind}_${split_kind}.pdf"

      echo "making $(basename "${output}") from ${split_count} PDFs"
      "${SCRIPT_DIR}/make_pdf_panel.sh" "${dir}" \
        --pattern "${pattern}" \
        --cols 3 \
        --row-gap-mm 0 \
        -o "${output}"
    done
  else
    output="${SCRIPT_DIR}/${rapidity}_${kind}.pdf"

    echo "making $(basename "${output}") from ${count} PDFs"
    "${SCRIPT_DIR}/make_pdf_panel.sh" "${dir}" \
      --cols 3 \
      --row-gap-mm 0 \
      -o "${output}"
  fi
done
