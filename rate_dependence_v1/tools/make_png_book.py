#!/usr/bin/env python3
"""
make_png_book.py

Create a single PDF that contains:
  1) global tau summary PNGs under:
       results/tau/
  2) selected per-setting PNGs under:
       results/**/PNGs/

Global tau PNGs are placed at the beginning of the PDF in this order
(if present):
  - tau_hist_all.png
  - tau_vs_current.png
  - tau_vs_shms34.png

Default per-setting PNGs included (if present):
  - yield_vs_current.png
  - yield_vs_trigger_shms34.png
  - hms_elclean_chargeNorm_vs_current.png
  - yield_ratio_vs_trigger_shms34.png

Ordering of settings:
  pass4, then pass5
    within pass: pi+ then pi-
      within pion: LH2 then LD2
        within target: z from low to high
          then tie-break by setting_id alphabetical

For each setting, a title page is inserted before its selected PNGs.
The tau summary PNGs are not setting-specific and are added without title pages.

Requires:
  ImageMagick 'convert' available in PATH.

Usage:
  cd rate_dependence_v1
  python3 tools/make_png_book.py --results results
  python3 tools/make_png_book.py --results results --out my_book.pdf

Optional:
  --include yield_vs_current.png
  --include yield_vs_trigger_shms34.png
  --include hms_elclean_chargeNorm_vs_current.png
  --include yield_ratio_vs_trigger_shms34.png

Notes:
  - --include controls only the per-setting PNGs under results/**/PNGs/.
  - The three tau summary PNGs are always considered separately and, if found,
    are prepended to the PDF in the fixed order listed above.
  - The output PDF is always written under results/PDFs/.
  - If --out is omitted, the default filename is all_settings_pngs.pdf.
"""

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict


TAU_SUMMARY_PNGS = [
    "tau_hist_all.png",
    "tau_vs_current.png",
    "tau_vs_shms34.png",
]


DEFAULT_INCLUDES = [
    "yield_vs_current.png",
    "yield_vs_trigger_shms34.png",
    "hms_elclean_chargeNorm_vs_current.png",
    "yield_ratio_vs_trigger_shms34.png",
]


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
    """Return 'convert' if available, else raise."""
    try:
        subprocess.check_call(["bash", "-lc", "command -v convert >/dev/null 2>&1"])
        return "convert"
    except subprocess.CalledProcessError:
        raise RuntimeError("ImageMagick 'convert' not found in PATH.")


def discover_settings_with_pngs(results_top: str):
    """
    Find all settings under results_top that have a PNGs/ directory.
    Returns dict: setting_dir -> png_dir
    where setting_dir is .../<setting_id> (parent of PNGs).
    """
    settings = {}
    for root, dirs, files in os.walk(results_top):
        if os.path.basename(root) == "PNGs":
            setting_dir = os.path.dirname(root)
            settings[setting_dir] = root
    return settings


def discover_tau_summary_pngs(results_top: str):
    """
    Return existing tau summary PNGs in the fixed order defined by
    TAU_SUMMARY_PNGS.
    """
    tau_dir = os.path.join(results_top, "tau")
    found = []
    for base in TAU_SUMMARY_PNGS:
        p = os.path.join(tau_dir, base)
        if os.path.isfile(p):
            found.append(p)
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results", help="Top results dir (default: results)")
    ap.add_argument("-o", "--out", default="all_settings_pngs.pdf",
                    help="Output PDF filename only (written under <results>/PDFs/). Default: all_settings_pngs.pdf")
    ap.add_argument("--density", type=int, default=150, help="Raster density for title pages/PDF (default: 150)")
    ap.add_argument("--title_pointsize", type=int, default=28, help="Font size for setting title pages (default: 28)")
    ap.add_argument("--quality", type=int, default=95, help="PDF quality/compression (default: 95)")
    ap.add_argument("--include", action="append", default=[],
                    help="PNG basename to include (repeatable). Default includes 4 per-setting files.")
    args = ap.parse_args()

    includes = args.include if args.include else list(DEFAULT_INCLUDES)

    convert = find_convert()

    tau_pages = discover_tau_summary_pngs(args.results)
    settings = discover_settings_with_pngs(args.results)
    if not settings:
        print(f"ERROR: No PNGs directories found under: {args.results}")
        sys.exit(1)

    groups = defaultdict(list)  # setting_dir -> [png_paths]
    total_found = 0
    for setting_dir, png_dir in settings.items():
        for base in includes:
            p = os.path.join(png_dir, base)
            if os.path.isfile(p):
                groups[setting_dir].append(p)
                total_found += 1

    groups = {sd: paths for sd, paths in groups.items() if paths}

    if not groups:
        print("ERROR: No matching PNGs found.")
        print("Searched for basenames:", ", ".join(includes))
        sys.exit(1)

    setting_dirs = sorted(groups.keys(), key=setting_sort_key)

    include_rank = {name: i for i, name in enumerate(includes)}
    for sd in setting_dirs:
        groups[sd] = sorted(
            groups[sd],
            key=lambda p: (include_rank.get(os.path.basename(p), 999), natural_key(p))
        )

    ordered_pages = list(tau_pages)
    for sd in setting_dirs:
        setting_id = os.path.basename(sd)
        ordered_pages.append(f"label:Setting: {setting_id}\\n{sd}")
        ordered_pages.extend(groups[sd])

    out_name = os.path.basename(args.out)
    if not out_name.lower().endswith(".pdf"):
        out_name += ".pdf"
    out_dir = os.path.join(args.results, "PDFs")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, out_name)

    cmd = [
        convert,
        "-density", str(args.density),
        "-background", "white",
        "-fill", "black",
        "-font", "Helvetica",
        "-pointsize", str(args.title_pointsize),
    ]
    cmd += ordered_pages
    cmd += ["-quality", str(args.quality), out_path]

    print(f"[INFO] Settings with selected PNGs: {len(setting_dirs)}")
    print(f"[INFO] Tau summary pages found:      {len(tau_pages)}")
    print(f"[INFO] Selected PNG pages found:   {total_found}")
    print(f"[INFO] Writing:                   {out_path}")
    print(f"[INFO] Running:                   convert (pages={len(ordered_pages)})")
    print(f"[INFO] Including files:           {', '.join(includes)}")

    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print("\nERROR: 'convert' failed.")
        print("If you see a 'security policy PDF' error, PDF writing may be disabled on this node.")
        print("In that case we can switch to a gs-based merge workflow.")
        sys.exit(e.returncode)

    print(f"[OK] Wrote: {out_path}")


if __name__ == "__main__":
    main()
