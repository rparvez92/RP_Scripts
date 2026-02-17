#!/usr/bin/env python3
"""
make_png_book.py

Create a single PDF that contains all PNGs under:
  results/**/PNGs/*.png

Ordering of settings:
  pass4, then pass5
    within pass: pi+ then pi-
      within pion: LH2 then LD2
        within target: z from low to high
          then tie-break by setting_id alphabetical

Requires:
  ImageMagick 'convert' available in PATH.

Usage:
  cd rate_dependence_v1
  python3 tools/make_png_book.py --results results -o results/all_settings_pngs.pdf
"""

import argparse
import glob
import os
import re
import subprocess
import sys
from collections import defaultdict


def natural_key(s: str):
    """Human-friendly sorting (file2 before file10)."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def parse_pass(s: str) -> int:
    sl = s.lower()
    if "pass4" in sl:
        return 0
    if "pass5" in sl:
        return 1
    return 9


def parse_pion(s: str) -> int:
    """
    Return 0 for pi+ (first), 1 for pi- (second), else 9.
    Tries to be robust to naming like:
      pi+sidis, piplus, pip, pi-
      piminus, pim
    """
    sl = s.lower()

    # Prefer explicit tokens
    if "pi+sidis" in sl or "piplus" in sl or "pi+" in sl:
        return 0
    if "pi-sidis" in sl or "piminus" in sl or "pi-" in sl:
        return 1

    # Fallback tokens
    if "pip" in sl and "pim" not in sl:
        return 0
    if "pim" in sl:
        return 1

    return 9


def parse_target(s: str) -> int:
    # LH2 first, then LD2
    if "LH2" in s:
        return 0
    if "LD2" in s:
        return 1
    return 9


def parse_z(s: str) -> float:
    """
    Extract z value from tokens like: z0p36, z0p5, z0p50, z0p3p6
    Returns +inf if not found.
    """
    m = re.search(r"\bz\d+p\d+(?:p\d+)?\b", s)
    if not m:
        return float("inf")

    token = m.group(0)  # e.g. z0p36
    val = token[1:]     # drop 'z' => 0p36
    parts = val.split("p")

    try:
        if len(parts) == 2:
            # 0p36 -> 0.36
            return float(parts[0] + "." + parts[1])
        if len(parts) == 3:
            # 0p3p6 -> 0.36
            return float(parts[0] + "." + parts[1] + parts[2])
    except Exception:
        return float("inf")

    return float("inf")


def setting_sort_key(setting_dir: str):
    """
    sorting tuple:
      (pass_order, pion_order, target_order, z_value, setting_id)
    """
    path_norm = setting_dir.replace("\\", "/")
    pass_order = parse_pass(path_norm)
    pion_order = parse_pion(path_norm)
    target_order = parse_target(setting_dir)  # keep LH2/LD2 case-sensitive check
    z_value = parse_z(path_norm)
    setting_id = os.path.basename(setting_dir).lower()
    return (pass_order, pion_order, target_order, z_value, setting_id)


def find_convert() -> str:
    """Return path to convert or raise."""
    try:
        subprocess.check_call(["bash", "-lc", "command -v convert >/dev/null 2>&1"])
        return "convert"
    except subprocess.CalledProcessError:
        raise RuntimeError("ImageMagick 'convert' not found in PATH.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results", help="Top results dir (default: results)")
    ap.add_argument("-o", "--out", default="results/all_settings_pngs.pdf", help="Output PDF path")
    ap.add_argument("--density", type=int, default=150, help="Raster density for title pages/PDF (default: 150)")
    ap.add_argument("--title_pointsize", type=int, default=28, help="Font size for setting title pages (default: 28)")
    ap.add_argument("--quality", type=int, default=95, help="PDF quality/compression (default: 95)")
    args = ap.parse_args()

    convert = find_convert()

    # Find PNGs
    pattern = os.path.join(args.results, "**", "PNGs", "*.png")
    pngs = glob.glob(pattern, recursive=True)
    if not pngs:
        print(f"ERROR: No PNGs found: {pattern}")
        sys.exit(1)

    # Group by setting directory (parent of PNGs/)
    groups = defaultdict(list)
    for p in pngs:
        setting_dir = os.path.dirname(os.path.dirname(p))  # .../<setting_id>
        groups[setting_dir].append(p)

    # Sort setting dirs with physics order
    setting_dirs = sorted(groups.keys(), key=setting_sort_key)

    # For each setting, sort images naturally
    for sd in setting_dirs:
        groups[sd] = sorted(groups[sd], key=natural_key)

    # Build ordered page list for convert, inserting a title page per setting
    ordered_pages = []
    for sd in setting_dirs:
        setting_id = os.path.basename(sd)
        # ImageMagick "label:" page
        ordered_pages.append(f"label:Setting: {setting_id}\\n{sd}")
        ordered_pages.extend(groups[sd])

    # Ensure output directory exists
    out_dir = os.path.dirname(args.out) or "."
    os.makedirs(out_dir, exist_ok=True)

    cmd = [
        convert,
        "-density", str(args.density),
        "-background", "white",
        "-fill", "black",
        "-font", "Helvetica",
        "-pointsize", str(args.title_pointsize),
    ]
    cmd += ordered_pages
    cmd += ["-quality", str(args.quality), args.out]

    print(f"[INFO] Found PNGs: {len(pngs)}")
    print(f"[INFO] Settings:   {len(setting_dirs)}")
    print(f"[INFO] Writing:    {args.out}")
    print(f"[INFO] Running:    convert (pages={len(ordered_pages)})")

    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print("\nERROR: 'convert' failed.")
        print("If you see a 'security policy PDF' error, PDF writing may be disabled on this node.")
        print("Paste the exact error line and I’ll give the gs-based fallback that always works.")
        sys.exit(e.returncode)

    print(f"[OK] Wrote: {args.out}")


if __name__ == "__main__":
    main()
