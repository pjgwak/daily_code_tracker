#!/usr/bin/env python3
"""Make one compact PDF page from several single-page plot PDFs.

Example:
    python3 make_contact_sheet_example.py

This intentionally starts with one folder and a few plots, as requested in
AGENTS.md.  It uses only stdlib Python plus the system `pdftoppm` command.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import tempfile
import zlib
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent.parent
DEFAULT_INPUT_DIR = PROJECT_DIR / "figs" / "y0.00_1.60" / "mass"
DEFAULT_OUTPUT = SCRIPT_DIR / "example_y0p00_1p60_mass_first6.pdf"


def mm_to_pt(value_mm: float) -> float:
    return value_mm * 72.0 / 25.4


def natural_plot_key(path: Path) -> tuple[float, float, str]:
    match = re.search(r"_pt([0-9]+(?:\.[0-9]+)?)_([0-9]+(?:\.[0-9]+)?)", path.name)
    if not match:
        return (math.inf, math.inf, path.name)
    return (float(match.group(1)), float(match.group(2)), path.name)


def read_ppm(path: Path) -> tuple[int, int, bytes]:
    data = path.read_bytes()
    pos = 0

    def next_token() -> bytes:
        nonlocal pos
        while pos < len(data) and data[pos] in b" \t\r\n":
            pos += 1
        if pos < len(data) and data[pos] == ord("#"):
            while pos < len(data) and data[pos] not in b"\r\n":
                pos += 1
            return next_token()
        start = pos
        while pos < len(data) and data[pos] not in b" \t\r\n":
            pos += 1
        return data[start:pos]

    magic = next_token()
    if magic != b"P6":
        raise ValueError(f"{path} is not a binary PPM file")
    width = int(next_token())
    height = int(next_token())
    max_value = int(next_token())
    if max_value != 255:
        raise ValueError(f"{path} has unsupported max color value {max_value}")
    while pos < len(data) and data[pos] in b" \t\r\n":
        pos += 1
    pixels = data[pos:]
    expected = width * height * 3
    if len(pixels) != expected:
        raise ValueError(f"{path} has {len(pixels)} bytes, expected {expected}")
    return width, height, pixels


def render_pdf_to_ppm(pdf_path: Path, output_stem: Path, dpi: int) -> Path:
    cmd = [
        "pdftoppm",
        "-f",
        "1",
        "-l",
        "1",
        "-singlefile",
        "-r",
        str(dpi),
        str(pdf_path),
        str(output_stem),
    ]
    subprocess.run(cmd, check=True)
    return output_stem.with_suffix(".ppm")


def pdf_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def write_pdf(
    output_path: Path,
    images: list[tuple[Path, int, int, bytes]],
    *,
    cols: int,
    page_width_mm: float,
    page_height_mm: float,
    margin_mm: float,
    gap_mm: float,
    title: str,
) -> None:
    page_w = mm_to_pt(page_width_mm)
    page_h = mm_to_pt(page_height_mm)
    margin = mm_to_pt(margin_mm)
    gap = mm_to_pt(gap_mm)
    title_h = mm_to_pt(8.0) if title else 0.0
    rows = math.ceil(len(images) / cols)
    cell_w = (page_w - 2.0 * margin - (cols - 1) * gap) / cols
    cell_h = (page_h - 2.0 * margin - title_h - (rows - 1) * gap) / rows

    commands = ["1 1 1 rg", f"0 0 {page_w:.3f} {page_h:.3f} re f"]
    if title:
        commands.extend(
            [
                "0 0 0 rg",
                "BT",
                "/F1 9 Tf",
                f"{margin:.3f} {page_h - margin - mm_to_pt(3.0):.3f} Td",
                f"({pdf_escape(title)}) Tj",
                "ET",
            ]
        )

    for idx, (source_path, img_w, img_h, _pixels) in enumerate(images, start=1):
        row = (idx - 1) // cols
        col = (idx - 1) % cols
        x0 = margin + col * (cell_w + gap)
        y_top = page_h - margin - title_h - row * (cell_h + gap)

        scale = min(cell_w / img_w, cell_h / img_h)
        draw_w = img_w * scale
        draw_h = img_h * scale
        x = x0 + (cell_w - draw_w) / 2.0
        y = y_top - cell_h + (cell_h - draw_h) / 2.0

        commands.extend(
            [
                "q",
                f"{draw_w:.3f} 0 0 {draw_h:.3f} {x:.3f} {y:.3f} cm",
                f"/Im{idx} Do",
                "Q",
            ]
        )

    content = ("\n".join(commands) + "\n").encode("ascii")

    objects: list[bytes] = []
    image_object_ids = []
    for idx, (_source_path, img_w, img_h, pixels) in enumerate(images, start=1):
        compressed = zlib.compress(pixels, level=6)
        image_object_ids.append(5 + idx - 1)
        objects.append(
            b"<< /Type /XObject /Subtype /Image "
            + f"/Width {img_w} /Height {img_h} ".encode("ascii")
            + b"/ColorSpace /DeviceRGB /BitsPerComponent 8 "
            + b"/Filter /FlateDecode "
            + f"/Length {len(compressed)} >>\nstream\n".encode("ascii")
            + compressed
            + b"\nendstream"
        )

    xobjects = " ".join(
        f"/Im{idx} {object_id} 0 R"
        for idx, object_id in enumerate(image_object_ids, start=1)
    )
    font_resource = "/Font << /F1 << /Type /Font /Subtype /Type1 /BaseFont /Helvetica >> >>" if title else ""
    resources = f"<< /XObject << {xobjects} >> {font_resource} >>"

    compressed_content = zlib.compress(content, level=6)
    base_objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        (
            f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {page_w:.3f} {page_h:.3f}] "
            f"/Resources {resources} /Contents 4 0 R >>"
        ).encode("ascii"),
        (
            f"<< /Filter /FlateDecode /Length {len(compressed_content)} >>\nstream\n"
        ).encode("ascii")
        + compressed_content
        + b"\nendstream",
    ]
    objects = base_objects + objects

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as handle:
        handle.write(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
        offsets = [0]
        for object_number, obj in enumerate(objects, start=1):
            offsets.append(handle.tell())
            handle.write(f"{object_number} 0 obj\n".encode("ascii"))
            handle.write(obj)
            handle.write(b"\nendobj\n")
        xref_offset = handle.tell()
        handle.write(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
        handle.write(b"0000000000 65535 f \n")
        for offset in offsets[1:]:
            handle.write(f"{offset:010d} 00000 n \n".encode("ascii"))
        handle.write(
            (
                f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\n"
                f"startxref\n{xref_offset}\n%%EOF\n"
            ).encode("ascii")
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create one below-A4 PDF page from a few plot PDFs."
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-plots", type=int, default=6)
    parser.add_argument("--cols", type=int, default=2)
    parser.add_argument("--render-dpi", type=int, default=120)
    parser.add_argument("--page-width-mm", type=float, default=190.0)
    parser.add_argument("--page-height-mm", type=float, default=260.0)
    parser.add_argument("--margin-mm", type=float, default=4.0)
    parser.add_argument("--gap-mm", type=float, default=2.0)
    parser.add_argument("--title", default="")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if shutil.which("pdftoppm") is None:
        raise SystemExit("ERROR: pdftoppm is required but was not found on PATH")
    if args.cols < 1:
        raise SystemExit("ERROR: --cols must be at least 1")
    if args.max_plots < 1:
        raise SystemExit("ERROR: --max-plots must be at least 1")

    pdfs = sorted(args.input_dir.glob("*.pdf"), key=natural_plot_key)[: args.max_plots]
    if not pdfs:
        raise SystemExit(f"ERROR: no PDF files found in {args.input_dir}")

    rendered_images = []
    with tempfile.TemporaryDirectory(prefix="contact_sheet_") as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        for idx, pdf_path in enumerate(pdfs, start=1):
            ppm_path = render_pdf_to_ppm(
                pdf_path, tmp_dir / f"plot_{idx:02d}", args.render_dpi
            )
            width, height, pixels = read_ppm(ppm_path)
            rendered_images.append((pdf_path, width, height, pixels))

    write_pdf(
        args.output,
        rendered_images,
        cols=args.cols,
        page_width_mm=args.page_width_mm,
        page_height_mm=args.page_height_mm,
        margin_mm=args.margin_mm,
        gap_mm=args.gap_mm,
        title=args.title,
    )

    print(f"Wrote {args.output}")
    print("Included:")
    for pdf_path, _width, _height, _pixels in rendered_images:
        print(f"  {pdf_path}")


if __name__ == "__main__":
    main()
