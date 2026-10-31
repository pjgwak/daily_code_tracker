#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  ./make_pdf_panel.sh INPUT_DIR [options]

Options:
  -o, --output FILE          Output PDF. Default: pdf_panel.pdf
  --pattern GLOB            Input PDF glob. Default: *.pdf
  --recursive               Search recursively.
  --first N                 Use only first N PDFs.
  --cols N                  Number of columns. Default: 3
  --page-width-mm X         Output page width. Default: 190
  --max-page-height-mm X    Maximum output page height. Default: 260
  --margin-mm X             Page margin. Default: 4
  --gap-mm X                Column gap. Default: 3
  --row-gap-mm X            Row gap. Default: 0
  --no-crop                 Do not crop input PDFs first.
  --keep-tex                Save generated .tex next to output PDF.
USAGE
}

die() {
  echo "error: $*" >&2
  exit 2
}

need_tool() {
  command -v "$1" >/dev/null 2>&1 || die "$1 not found"
}

abspath() {
  local path="$1"
  local dir base
  dir="$(dirname "$path")"
  base="$(basename "$path")"
  if [ -d "$dir" ]; then
    printf '%s/%s\n' "$(cd "$dir" && pwd -P)" "$base"
  else
    printf '%s/%s\n' "$(cd "$(dirname "$dir")" && pwd -P)" "$(basename "$dir")/$base"
  fi
}

calc() {
  awk "BEGIN { printf \"%.6f\", $* }"
}

tex_path_arg() {
  local path="$1"
  path="${path//\\/\/}"
  printf '\\detokenize{%s}' "$path"
}

bbox_aspect() {
  local pdf="$1"
  gs -q -dBATCH -dNOPAUSE -sDEVICE=bbox "$pdf" 2>&1 |
    awk '
      /%%HiResBoundingBox:/ { x0=$2; y0=$3; x1=$4; y1=$5 }
      END {
        w = x1 - x0
        h = y1 - y0
        if (w <= 0 || h <= 0) exit 1
        printf "%.10f", w / h
      }
    '
}

input_dir="."
output="pdf_panel.pdf"
pattern="*.pdf"
recursive=0
first=""
cols=3
page_width_mm=190
max_page_height_mm=260
margin_mm=4
gap_mm=3
row_gap_mm=0
crop_inputs=1
keep_tex=0

if [ "$#" -gt 0 ] && [[ "$1" != -* ]]; then
  input_dir="$1"
  shift
fi

while [ "$#" -gt 0 ]; do
  case "$1" in
    -o|--output)
      output="${2:?missing output path}"
      shift 2
      ;;
    --pattern)
      pattern="${2:?missing pattern}"
      shift 2
      ;;
    --recursive)
      recursive=1
      shift
      ;;
    --first)
      first="${2:?missing first value}"
      shift 2
      ;;
    --cols)
      cols="${2:?missing cols value}"
      shift 2
      ;;
    --page-width-mm)
      page_width_mm="${2:?missing page width}"
      shift 2
      ;;
    --max-page-height-mm)
      max_page_height_mm="${2:?missing max page height}"
      shift 2
      ;;
    --margin-mm)
      margin_mm="${2:?missing margin}"
      shift 2
      ;;
    --gap-mm)
      gap_mm="${2:?missing gap}"
      shift 2
      ;;
    --row-gap-mm)
      row_gap_mm="${2:?missing row gap}"
      shift 2
      ;;
    --no-crop)
      crop_inputs=0
      shift
      ;;
    --keep-tex)
      keep_tex=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown option: $1"
      ;;
  esac
done

[ -d "$input_dir" ] || die "input directory does not exist: $input_dir"
[ "$cols" -ge 1 ] || die "--cols must be >= 1"

need_tool gs
need_tool pdflatex
if [ "$crop_inputs" -eq 1 ]; then
  need_tool pdfcrop
fi

output_abs="$(abspath "$output")"
tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/pdf_panel.XXXXXX")"
trap 'rm -rf "$tmpdir"' EXIT
crop_dir="$tmpdir"
if [ "$keep_tex" -eq 1 ] && [ "$crop_inputs" -eq 1 ]; then
  crop_dir="${output_abs%.pdf}_inputs"
  mkdir -p "$crop_dir"
fi

pdfs=()
find_args=("$input_dir")
if [ "$recursive" -eq 0 ]; then
  find_args+=(-maxdepth 1)
fi
while IFS= read -r pdf; do
  pdf_abs="$(abspath "$pdf")"
  [ "$pdf_abs" = "$output_abs" ] && continue
  pdfs+=("$pdf_abs")
  if [ -n "$first" ] && [ "${#pdfs[@]}" -ge "$first" ]; then
    break
  fi
done < <(find "${find_args[@]}" -type f -name "$pattern" -print | sort -V)

[ "${#pdfs[@]}" -gt 0 ] || die "no PDFs found in $input_dir with pattern $pattern"

items=()
aspects=()
idx=1
for pdf in "${pdfs[@]}"; do
  item_pdf="$pdf"
  if [ "$crop_inputs" -eq 1 ]; then
    stem="$(basename "$pdf" .pdf)"
    stem="$(printf '%s' "$stem" | sed 's/[^A-Za-z0-9_.-]/_/g')"
    item_pdf="$crop_dir/$(printf '%03d_%s.pdf' "$idx" "$stem")"
    pdfcrop --margins 0 "$pdf" "$item_pdf" >/dev/null
  fi
  aspect="$(bbox_aspect "$item_pdf")"
  items+=("$item_pdf")
  aspects+=("$aspect")
  idx=$((idx + 1))
done

cell_w="$(calc "($page_width_mm - 2 * $margin_mm - $gap_mm * ($cols - 1)) / $cols")"
rows=$(( (${#items[@]} + cols - 1) / cols ))

row_heights=()
total_rows_h=0
for ((r = 0; r < rows; r++)); do
  row_h=0
  for ((c = 0; c < cols; c++)); do
    i=$((r * cols + c))
    [ "$i" -lt "${#items[@]}" ] || continue
    h="$(calc "$cell_w / ${aspects[$i]}")"
    row_h="$(awk -v a="$row_h" -v b="$h" 'BEGIN { printf "%.6f", (a > b ? a : b) }')"
  done
  row_heights+=("$row_h")
  total_rows_h="$(calc "$total_rows_h + $row_h")"
done

natural_page_h="$(calc "2 * $margin_mm + $total_rows_h + $row_gap_mm * ($rows - 1) + 0.8")"
page_height_mm="$(awk -v h="$natural_page_h" -v max="$max_page_height_mm" 'BEGIN { printf "%.6f", (h < max ? h : max) }')"
available_h="$(calc "$page_height_mm - 2 * $margin_mm - $row_gap_mm * ($rows - 1) - 0.8")"
scale="$(awk -v total="$total_rows_h" -v avail="$available_h" 'BEGIN { if (total > avail && total > 0) printf "%.10f", avail / total; else printf "1.0000000000" }')"

if awk -v s="$scale" 'BEGIN { exit !(s < 0.999999) }'; then
  for ((r = 0; r < rows; r++)); do
    row_heights[$r]="$(calc "${row_heights[$r]} * $scale")"
  done
fi

tex="$tmpdir/panel.tex"
{
  printf '\\documentclass{article}\n'
  printf '\\usepackage[paperwidth=%.6fmm,paperheight=%.6fmm,margin=%.6fmm]{geometry}\n' \
    "$page_width_mm" "$page_height_mm" "$margin_mm"
  printf '\\usepackage{graphicx}\n'
  printf '\\setkeys{Gin}{draft=false}\n'
  printf '\\pagestyle{empty}\n'
  printf '\\setlength{\\parindent}{0pt}\n'
  printf '\\setlength{\\parskip}{0pt}\n'
  printf '\\setlength{\\topskip}{0pt}\n'
  printf '\\begin{document}\n'
  printf '\\noindent%%\n'

  for ((r = 0; r < rows; r++)); do
    row_h="${row_heights[$r]}"
    for ((c = 0; c < cols; c++)); do
      i=$((r * cols + c))
      if [ "$i" -lt "${#items[@]}" ]; then
        printf '\\begin{minipage}[t]{%.6fmm}\n' "$cell_w"
        printf '\\centering\n'
        printf '\\includegraphics[page=1,width=%.6fmm,height=%.6fmm,keepaspectratio]{%s}\n' \
          "$cell_w" "$row_h" "$(tex_path_arg "${items[$i]}")"
        printf '\\end{minipage}%%\n'
      else
        printf '\\hspace*{%.6fmm}%%\n' "$cell_w"
      fi

      if [ "$c" -ne $((cols - 1)) ]; then
        printf '\\hspace*{%.6fmm}%%\n' "$gap_mm"
      fi
    done
    if [ "$r" -ne $((rows - 1)) ]; then
      printf '\\par\\nointerlineskip\\vspace*{%.6fmm}%%\n' "$row_gap_mm"
    fi
  done

  printf '\\end{document}\n'
} > "$tex"

pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$tmpdir" "$tex" >/dev/null
mkdir -p "$(dirname "$output_abs")"
cp "$tmpdir/panel.pdf" "$output_abs"
if [ "$keep_tex" -eq 1 ]; then
  cp "$tex" "${output_abs%.pdf}.tex"
fi

echo "wrote $output_abs"
echo "used ${#items[@]} PDFs from $(cd "$input_dir" && pwd -P)"
