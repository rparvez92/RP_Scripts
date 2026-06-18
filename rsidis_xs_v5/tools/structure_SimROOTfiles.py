#!/usr/bin/env python3
"""
Populate the v5 Pass0_SimROOTfiles/ tree with symlinks to local SIMC files.

TableCoinXsec.C expects SIMC files in:

  Pass0_SimROOTfiles/<setting_id>/<stem>.root
  Pass0_SimROOTfiles/<setting_id>/<stem>.hist

where <setting_id> is the path below settings/ that contains manifest.json.
On the local Mac/T7 setup, the ROOT files are flat in
/Volumes/T7/RSIDIS/Pass0_SimROOTfiles and the matching SIMC .hist files are
flat in /Volumes/T7/RSIDIS/Pass0_SimOUTfiles. This tool recreates the expected
project-local layout with symlinks, without copying the large ROOT files.

Examples, from rsidis_xs_v5/:

  # Preview only
  python3 tools/structure_SimROOTfiles.py

  # Create symlinks
  python3 tools/structure_SimROOTfiles.py --execute

  # Replace existing wrong/broken symlinks
  python3 tools/structure_SimROOTfiles.py --execute --replace
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


RUN_TYPE_TO_STEM = {
    "pi+sidis": "piplus",
    "pi-sidis": "piminus",
}


@dataclass(frozen=True)
class SettingLinkPlan:
    manifest: Path
    setting_id: Path
    stem: str
    root_src: Path
    hist_src: Path
    root_dst: Path
    hist_dst: Path
    root_exists: bool
    hist_exists: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create Pass0_SimROOTfiles/<setting_id>/ symlinks for v5 SIMC ROOT/hist files."
    )
    parser.add_argument(
        "--settings-root",
        default="settings",
        type=Path,
        help="Settings tree containing manifest.json files. Default: settings",
    )
    parser.add_argument(
        "--sim-root",
        default=Path("/Volumes/T7/RSIDIS/Pass0_SimROOTfiles"),
        type=Path,
        help="Flat directory containing SIMC .root files. Default: /Volumes/T7/RSIDIS/Pass0_SimROOTfiles",
    )
    parser.add_argument(
        "--hist-root",
        default=Path("/Volumes/T7/RSIDIS/Pass0_SimOUTfiles"),
        type=Path,
        help="Flat directory containing SIMC .hist files. Default: /Volumes/T7/RSIDIS/Pass0_SimOUTfiles",
    )
    parser.add_argument(
        "--out-root",
        default=Path("Pass0_SimROOTfiles"),
        type=Path,
        help="Project-local SIMC output tree. Default: Pass0_SimROOTfiles",
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually create symlinks. Without this flag, the tool only reports a dry run.",
    )
    parser.add_argument(
        "--replace",
        action="store_true",
        help="Replace existing symlinks/files at destination paths. Use with care.",
    )
    parser.add_argument(
        "--show-missing",
        type=int,
        default=30,
        help="Maximum missing-setting examples to print. Default: 30",
    )
    return parser.parse_args()


def load_manifest(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def setting_id_from_manifest(manifest_path: Path, settings_root: Path) -> Path:
    return manifest_path.parent.relative_to(settings_root)


def stem_from_manifest(manifest_path: Path, settings_root: Path) -> Optional[str]:
    manifest = load_manifest(manifest_path)
    keys = manifest.get("path_keys", {})
    if not isinstance(keys, dict):
        return None

    pass_key = str(keys.get("pass", ""))
    run_type = str(keys.get("run_type", ""))
    target = str(keys.get("target", ""))
    setting_leaf = str(keys.get("setting_id", ""))

    pion_token = RUN_TYPE_TO_STEM.get(run_type)
    if not pass_key or not pion_token or not target or not setting_leaf:
        return None

    return f"coin_{pass_key}_{pion_token}_{target}_{setting_leaf}"


def build_plans(
    settings_root: Path,
    sim_root: Path,
    hist_root: Path,
    out_root: Path,
) -> Tuple[List[SettingLinkPlan], List[Tuple[Path, str]]]:
    plans: List[SettingLinkPlan] = []
    skipped: List[Tuple[Path, str]] = []

    for manifest_path in sorted(settings_root.rglob("manifest.json")):
        stem = stem_from_manifest(manifest_path, settings_root)
        if not stem:
            skipped.append((manifest_path, "unsupported or incomplete path_keys"))
            continue

        setting_id = setting_id_from_manifest(manifest_path, settings_root)
        root_src = sim_root / f"{stem}.root"
        hist_src = hist_root / f"{stem}.hist"
        dst_dir = out_root / setting_id

        plans.append(
            SettingLinkPlan(
                manifest=manifest_path,
                setting_id=setting_id,
                stem=stem,
                root_src=root_src,
                hist_src=hist_src,
                root_dst=dst_dir / f"{stem}.root",
                hist_dst=dst_dir / f"{stem}.hist",
                root_exists=root_src.exists(),
                hist_exists=hist_src.exists(),
            )
        )

    return plans, skipped


def link_one(src: Path, dst: Path, replace: bool) -> str:
    if dst.exists() or dst.is_symlink():
        if not replace:
            if dst.is_symlink() and Path(os.readlink(dst)) == src:
                return "already_ok"
            return "exists"
        dst.unlink()

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.symlink_to(src)
    return "linked"


def execute(plans: Iterable[SettingLinkPlan], replace: bool) -> Dict[str, int]:
    counts = {
        "root_linked": 0,
        "hist_linked": 0,
        "already_ok": 0,
        "exists": 0,
        "missing_source": 0,
    }

    for plan in plans:
        if not (plan.root_exists and plan.hist_exists):
            counts["missing_source"] += 1
            continue

        root_status = link_one(plan.root_src, plan.root_dst, replace)
        hist_status = link_one(plan.hist_src, plan.hist_dst, replace)

        for status, key in ((root_status, "root_linked"), (hist_status, "hist_linked")):
            if status == "linked":
                counts[key] += 1
            elif status == "already_ok":
                counts["already_ok"] += 1
            elif status == "exists":
                counts["exists"] += 1

    return counts


def main() -> int:
    args = parse_args()
    settings_root = args.settings_root
    sim_root = args.sim_root
    hist_root = args.hist_root
    out_root = args.out_root

    if not settings_root.is_dir():
        raise SystemExit(f"ERROR: settings root not found: {settings_root}")
    if not sim_root.is_dir():
        raise SystemExit(f"ERROR: SIM ROOT directory not found: {sim_root}")
    if not hist_root.is_dir():
        raise SystemExit(f"ERROR: SIM hist directory not found: {hist_root}")

    plans, skipped = build_plans(settings_root, sim_root, hist_root, out_root)
    complete = [p for p in plans if p.root_exists and p.hist_exists]
    missing_root = [p for p in plans if not p.root_exists]
    missing_hist = [p for p in plans if p.root_exists and not p.hist_exists]

    mode = "EXECUTE" if args.execute else "DRY RUN"
    print(f"[{mode}] settings_root={settings_root}")
    print(f"[{mode}] sim_root={sim_root}")
    print(f"[{mode}] hist_root={hist_root}")
    print(f"[{mode}] out_root={out_root}")
    print()
    print(f"manifests considered: {len(plans) + len(skipped)}")
    print(f"supported manifests:  {len(plans)}")
    print(f"complete root+hist:   {len(complete)}")
    print(f"missing root:         {len(missing_root)}")
    print(f"missing hist only:    {len(missing_hist)}")
    print(f"skipped manifests:    {len(skipped)}")

    if missing_root or missing_hist or skipped:
        print()
        print("Examples needing attention:")
        shown = 0
        for plan in missing_root[: args.show_missing]:
            print(f"  missing root: {plan.setting_id} -> {plan.root_src.name}")
            shown += 1
        remaining = max(args.show_missing - shown, 0)
        for plan in missing_hist[:remaining]:
            print(f"  missing hist: {plan.setting_id} -> {plan.hist_src.name}")
            shown += 1
        remaining = max(args.show_missing - shown, 0)
        for manifest, reason in skipped[:remaining]:
            print(f"  skipped: {manifest} ({reason})")

    if args.execute:
        print()
        counts = execute(plans, args.replace)
        print("link results:")
        for key, value in counts.items():
            print(f"  {key}: {value}")
    else:
        print()
        print("No links were created. Re-run with --execute to populate Pass0_SimROOTfiles/.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
