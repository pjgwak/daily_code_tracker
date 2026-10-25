#!/usr/bin/env bash
#######################
# Uses LaTeX when available; otherwise Poppler/Ghostscript vector composition.
#######################
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
FIGS_DIR="${FIGS_DIR:-${SCRIPT_DIR}/../../figs_publish}"

if command -v pdflatex >/dev/null 2>&1 && command -v pdfcrop >/dev/null 2>&1; then
  panel_command=(bash "${SCRIPT_DIR}/make_pdf_panel.sh")
else
  panel_command=(python3 "${SCRIPT_DIR}/make_pdf_panel.py")
fi

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
      "${panel_command[@]}" "${dir}" \
        --pattern "${pattern}" \
        --cols 3 \
        --row-gap-mm 0 \
        -o "${output}"
    done
  else
    output="${SCRIPT_DIR}/${rapidity}_${kind}.pdf"

    echo "making $(basename "${output}") from ${count} PDFs"
    "${panel_command[@]}" "${dir}" \
      --cols 3 \
      --row-gap-mm 0 \
      -o "${output}"
  fi
done
