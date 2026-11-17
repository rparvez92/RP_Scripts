#!/usr/bin/env python3
"""
Batch driver: walk manifests and run a ROOT macro per setting.

- Mirrors output into a parallel results tree
- Writes output_meta.json with corrections_enabled and manifest path

Use --dry-run to print the ROOT commands without executing them.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import List


def find_manifests(settings_root: Path) -> List[Path]:
    return sorted(settings_root.rglob("manifest.json"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--settings", required=True, help="settings root (e.g., rsidis_xs_v5/settings)")
    ap.add_argument("--results", required=True, help="results root (e.g., rsidis_xs_v5/results)")
    ap.add_argument("--root-macro", default="PlotFromManifest.C", help="ROOT macro to execute")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print commands only (no ROOT execution). Useful for verifying scope.")
    ap.add_argument("--corrections-enabled", default="comp_livetime",
                    help="Comma-separated list recorded in output_meta.json (e.g., comp_livetime,boil_corr,fan_mean)")
    args = ap.parse_args()

    settings_root = Path(args.settings).expanduser().resolve()
    results_root = Path(args.results).expanduser().resolve()
    results_root.mkdir(parents=True, exist_ok=True)

    enabled = [x.strip() for x in args.corrections_enabled.split(",") if x.strip()]

    manifests = find_manifests(settings_root)
    if not manifests:
        print(f"No manifest.json found under {settings_root}")
        return

    for mpath in manifests:
        manifest = json.loads(mpath.read_text())
        setting_path = Path(manifest["setting"]["path"])
        rel = setting_path.relative_to(settings_root)
        out_dir = results_root / rel
        out_dir.mkdir(parents=True, exist_ok=True)

        (out_dir / "output_meta.json").write_text(json.dumps({
            "manifest": str(mpath),
            "corrections_enabled": enabled,
        }, indent=2) + "\n")

        cmd = ["root", "-l", "-b", "-q", f'{args.root_macro}("{str(mpath)}","{str(out_dir)}")']
        print(" ".join(cmd))

        if not args.dry_run:
            subprocess.run(cmd, check=True)

    print(f"Processed {len(manifests)} manifests. Results in: {results_root}")


if __name__ == "__main__":
    main()
