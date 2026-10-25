#!/usr/bin/env python3
"""Compose vector PDF panels using installed Poppler/Ghostscript (no LaTeX)."""
import argparse
import concurrent.futures
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_dir', type=Path)
    parser.add_argument('--pattern', default='*.pdf')
    parser.add_argument('--cols', type=int, default=3)
    parser.add_argument('--row-gap-mm', type=float, default=0.)
    parser.add_argument('-o', '--output', type=Path, required=True)
    args = parser.parse_args()
    if args.cols < 1 or args.row_gap_mm < 0:
        parser.error('cols must be positive and row gap nonnegative')
    def order(path):
        match=re.search(r'_pt([0-9.]+)_([0-9.]+)\.pdf$', path.name)
        return (float(match[1]), float(match[2])) if match else (math.inf, path.name)
    inputs=sorted((p.resolve() for p in args.input_dir.glob(args.pattern)
                   if p.is_file() and p.resolve()!=args.output.resolve()), key=order)
    if not inputs:
        parser.error('no matching PDFs')
    with tempfile.TemporaryDirectory(prefix='oo_pdf_panel_') as tmp:
        work=Path(tmp)
        def prepare(pair):
            i, source=pair
            eps=work/f'{i}.eps'
            subprocess.run(['pdftops','-f','1','-l','1','-eps',str(source),str(eps)],check=True)
            result=subprocess.run(['gs','-q','-dSAFER','-dBATCH','-dNOPAUSE','-sDEVICE=bbox',str(eps)],check=True,capture_output=True,text=True)
            boxes=re.findall(r'%%HiResBoundingBox: ([^\n]+)',result.stderr)
            if not boxes:raise RuntimeError(f'Missing drawing bounds: {source}')
            x0,y0,x1,y1=map(float,boxes[-1].split())
            if x1<=x0 or y1<=y0:raise RuntimeError(f'Empty PDF: {source}')
            return eps,(x0,y0,x1,y1)
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            items=list(pool.map(prepare,enumerate(inputs)))
        mm=72/25.4; width=190*mm; margin=4*mm; gap=3*mm; rowgap=args.row_gap_mm*mm
        cell=(width-2*margin-(args.cols-1)*gap)/args.cols
        heights=[max(cell*(b[3]-b[1])/(b[2]-b[0]) for _,b in items[i:i+args.cols])
                 for i in range(0,len(items),args.cols)]
        available=260*mm-2*margin-rowgap*(len(heights)-1)
        scale=min(1.,available/sum(heights));heights=[h*scale for h in heights]
        height=2*margin+sum(heights)+rowgap*(len(heights)-1)
        ps=work/'panel.ps'
        with ps.open('wb') as out:
            def emit(s):out.write(s.encode('ascii'))
            emit(f'%!PS-Adobe-3.0\n%%Pages: 1\n<< /PageSize [{width} {height}] >> setpagedevice\n')
            top=height-margin
            for i,(eps,b) in enumerate(items):
                row,col=divmod(i,args.cols)
                if col==0 and row:top-=heights[row-1]+rowgap
                factor=min(cell/(b[2]-b[0]),heights[row]/(b[3]-b[1]))
                x=margin+col*(cell+gap)+(cell-factor*(b[2]-b[0]))/2
                y=top-factor*(b[3]-b[1])
                emit(f'save\n/showpage {{}} def\n{x} {y} translate\n{factor} {factor} scale\n{-b[0]} {-b[1]} translate\n%%BeginDocument: {i}.eps\n')
                out.write(eps.read_bytes())
                emit('\n%%EndDocument\nrestore\n')
            emit('showpage\n%%EOF\n')
        target=work/'panel.pdf'
        subprocess.run(['gs','-q','-dSAFER','-dBATCH','-dNOPAUSE','-sDEVICE=pdfwrite',f'-sOutputFile={target}',str(ps)],check=True)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_bytes(target.read_bytes())
        # args.output.with_suffix('.inputs.json').write_text(json.dumps([str(p) for p in inputs],indent=2)+'\n')
    print(f'wrote {args.output}: {len(inputs)} vector panels')


if __name__=='__main__':
    main()
