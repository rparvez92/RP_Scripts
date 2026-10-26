#!/usr/bin/env python3
"""Combine corrected SIDIS data and MC yields across nominal-thpq settings."""

from __future__ import annotations

import argparse
import math
from collections import Counter, defaultdict
from pathlib import Path

import ROOT

import RP_data_mc_compare as dmc
import RP_azimuthal_coverage_study as coverage
from RP_build_thpq_combination import MAP_FIELDS, POLICY_VERSION
from RP_sidis_binning import AZIMUTHAL_BINNING_VERSION, AZIMUTHAL_PHI_EDGES, AZIMUTHAL_PT_EDGES

DEFAULT_MAP = dmc.PROJECT_DIR / "results/Pass1/Tables/RP_thpq_combination/RP_thpq_combination_map.csv"
DEFAULT_OUTPUT_DIR = dmc.PROJECT_DIR / "results/Pass1/Tables/RP_combine_thpq_yields"
FIELDS = [
    "Combination_group_id", "Combination_policy_version", "Phase", "Pass",
    "Run_type", "Target", "x", "Q2", "z", "Ebeam_GeV", "N_contributors",
    "Contributing_setting_ids", "Contributing_thpq_deg", "Pt_bin", "Pt_low",
    "Pt_high", "Pt_high_closed", "phipq_bin", "phipq_low", "phipq_high",
    "phipq_center", "Binning_version", "Y_Data", "Y_Data_err",
    *(item for reaction in dmc.REACTIONS for item in
      (f"Y_MC_{reaction}", f"Y_MC_{reaction}_err")),
    "Y_MC", "Y_MC_err", "Data_by_MC", "Data_by_MC_err",
    "Combination_status", "Combination_reason",
]


def finite(value: object) -> bool:
    try: return math.isfinite(float(value))
    except (TypeError, ValueError): return False


def combine_group(group_id, members, norm_groups, catalog, root_dir):
    first = members[0]
    data_total = coverage.make_hist(f"{group_id}_data")
    components = {reaction: coverage.make_hist(f"{group_id}_{reaction}")
                  for reaction in dmc.REACTIONS}
    structural, advisory = [], []
    for member in members:
        key = (member["Phase"], member["Pass"], member["Run_type"],
               member["Target"], member["setting_id"])
        base = (key[0], key[1], key[2], key[4])
        files = norm_groups.get(base, {})
        categories = {
            "target_e": (key[3], "Elec"), "target_p": (key[3], "Pos"),
            "dummy_e": ("DUMMY", "Elec"), "dummy_p": ("DUMMY", "Pos"),
        }
        classes = {}
        for name, category in categories.items():
            if category in files:
                classes[name], reasons = coverage.data_class(
                    files[category], root_dir, f"{group_id}_{key[4]}_{name}")
                structural.extend(f"{key[4]}:{reason}" for reason in reasons)
            else:
                classes[name] = coverage.make_hist(f"{group_id}_{key[4]}_{name}_zero")
                advisory.append(f"{key[4]}:ZERO_CORRECTION_CLASS={name}")
        setting_data = coverage.combine_data(classes, dmc.S_DUMMY[(key[0], key[3])])
        data_total.Add(setting_data)
        physical = dmc.applicable_reactions(key[2])
        for reaction in dmc.REACTIONS:
            if reaction not in physical:
                continue
            source = catalog.get((*key, reaction))
            if source is None:
                structural.append(f"{key[4]}:MISSING_MC={reaction}")
                continue
            hist, reasons = coverage.simc_class(
                source, f"{group_id}_{key[4]}_{reaction}")
            structural.extend(f"{key[4]}:{reason}" for reason in reasons)
            components[reaction].Add(hist)
    reasons = sorted(set(structural + advisory))
    status = "ERROR" if structural else "WARNING" if advisory else "OK"
    settings = "|".join(member["setting_id"] for member in members)
    thpqs = "|".join(member["thpq_nominal_deg"] for member in members)
    rows = []
    for ipt, (pt_low, pt_high) in enumerate(zip(AZIMUTHAL_PT_EDGES[:-1], AZIMUTHAL_PT_EDGES[1:]), 1):
        for iphi, (phi_low, phi_high) in enumerate(zip(AZIMUTHAL_PHI_EDGES[:-1], AZIMUTHAL_PHI_EDGES[1:]), 1):
            row = {field: "" for field in FIELDS}
            data, data_error = data_total.GetBinContent(iphi, ipt), data_total.GetBinError(iphi, ipt)
            mc_values = [components[r].GetBinContent(iphi, ipt) for r in dmc.REACTIONS]
            mc_errors = [components[r].GetBinError(iphi, ipt) for r in dmc.REACTIONS]
            mc, mc_error = sum(mc_values), math.sqrt(sum(value * value for value in mc_errors))
            ratio = data / mc if mc > 0 else math.nan
            ratio_error = math.hypot(data_error / mc, data * mc_error / (mc * mc)) if mc > 0 else math.nan
            row.update({
                "Combination_group_id": group_id,
                "Combination_policy_version": POLICY_VERSION,
                **{field: first[field] for field in
                   ("Phase", "Pass", "Run_type", "Target", "x", "Q2", "z", "Ebeam_GeV")},
                "N_contributors": len(members), "Contributing_setting_ids": settings,
                "Contributing_thpq_deg": thpqs, "Pt_bin": ipt, "Pt_low": pt_low,
                "Pt_high": pt_high, "Pt_high_closed": int(ipt == len(AZIMUTHAL_PT_EDGES)-1),
                "phipq_bin": iphi, "phipq_low": phi_low, "phipq_high": phi_high,
                "phipq_center": (phi_low + phi_high) / 2,
                "Binning_version": AZIMUTHAL_BINNING_VERSION,
                "Y_Data": data, "Y_Data_err": data_error, "Y_MC": mc,
                "Y_MC_err": mc_error, "Data_by_MC": ratio,
                "Data_by_MC_err": ratio_error, "Combination_status": status,
                "Combination_reason": ";".join(reasons),
            })
            for reaction, value, error in zip(dmc.REACTIONS, mc_values, mc_errors):
                row[f"Y_MC_{reaction}"] = value; row[f"Y_MC_{reaction}_err"] = error
            rows.append(row)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--combination-map", type=Path, default=DEFAULT_MAP)
    parser.add_argument("--norm-dir", type=Path, default=dmc.DEFAULT_NORM_DIR)
    parser.add_argument("--simc-catalog", type=Path, default=dmc.DEFAULT_SIMC_CATALOG)
    parser.add_argument("--t7-root", type=Path, default=dmc.DEFAULT_T7_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--group", action="append", default=[])
    args = parser.parse_args()
    mapped = [row for row in dmc.read_csv(args.combination_map)
              if row.get("Combination_status") == "OK"]
    grouped = defaultdict(list)
    for row in mapped: grouped[row["Combination_group_id"]].append(row)
    if args.group: grouped = {key: value for key, value in grouped.items() if key in set(args.group)}
    norm_groups = dmc.discover_norm_files(args.norm_dir)
    catalog = {
        tuple(row[field] for field in ("Phase", "Pass", "Run_type", "Target", "setting_id", "Reaction")): row
        for row in dmc.read_csv(args.simc_catalog)
    }
    rows = []
    for index, (group_id, members) in enumerate(sorted(grouped.items()), 1):
        members.sort(key=lambda row: int(row["Contributor_index"]))
        found = combine_group(group_id, members, norm_groups, catalog, args.t7_root)
        rows.extend(found)
        print(f"[{index}/{len(grouped)}] {found[0]['Combination_status']}: {group_id}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    coverage.atomic_csv(args.output_dir / "RP_combined_thpq_yields.csv", FIELDS, rows)
    problematic = [row for row in rows if row["Combination_status"] != "OK"]
    coverage.atomic_csv(args.output_dir / "RP_combined_thpq_yields_Problematic.csv", FIELDS, problematic)
    print("Rows:", len(rows), "Groups:", len(grouped),
          "Status:", Counter(row["Combination_status"] for row in rows))
    return 1 if any(row["Combination_status"] == "ERROR" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
