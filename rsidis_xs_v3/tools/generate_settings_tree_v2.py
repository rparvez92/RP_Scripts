#!/usr/bin/env python3
"""
Generate RSIDIS settings directory tree + manifests/runlists from rsidis_bigtable_pass0.csv.

Tree (human browsing UI):
  Pass -> RunType -> Target(LH2/LD2) -> z -> x -> Q2 -> SettingID

Per-setting files (deterministic truth for processing):
  - manifest.json
  - README.txt
  - runs_signal.txt
  - runs_dummy.txt
  - runs_positron.txt
  - run_metadata.csv

Rules:
- Case-insensitive matching for categorical fields (target, run_type)
- Placeholder rows with x/Q2/z == -999 are excluded
- Dummy runs are matched into each LH2/LD2 setting (no Dummy folder branch)
- Positron runs are those with HMS momentum sign > 0 (hms_p > 0)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


SCHEMA_VERSION = "1.2"


def norm_str(x: object) -> str:
    return "" if x is None else str(x).strip().upper()


def token_from_numeric_string(raw: str) -> str:
    s = raw.strip()
    if s.startswith("+"):
        s = s[1:]
    return s.replace(".", "p")


def signed_token(raw: str) -> str:
    s = raw.strip()
    if s.startswith("+"):
        s = s[1:]
    if s.startswith("-"):
        return "neg" + token_from_numeric_string(s[1:])
    return token_from_numeric_string(s)


def safe_float(raw: object) -> float:
    try:
        return float(raw)
    except Exception:
        return float("nan")


def beam_pass_key(ebeam: float, tol: float = 1e-3) -> Tuple[str, str]:
    if math.isfinite(ebeam) and abs(ebeam - 8.5831) < tol:
        return ("pass4", "E8p5831")
    if math.isfinite(ebeam) and abs(ebeam - 10.6716) < tol:
        return ("pass5", "E10p6716")
    if not math.isfinite(ebeam):
        return ("passUNK", "EUNK")
    s = f"{ebeam:.4f}".replace(".", "p")
    return ("passUNK", f"E{s}")


def choose_ps_factor(ps5: float, ps6: float, tol: float = 1e-6) -> Tuple[str, float, List[str]]:
    warnings: List[str] = []
    is5 = math.isfinite(ps5) and abs(ps5 - 1.0) < tol
    is6 = math.isfinite(ps6) and abs(ps6 - 1.0) < tol
    if is5 and not is6:
        return ("ps5", ps5, warnings)
    if is6 and not is5:
        return ("ps6", ps6, warnings)
    if is5 and is6:
        return ("ps5", ps5, warnings)
    warnings.append(f"Neither ps5 nor ps6 is ~1.0 (ps5={ps5}, ps6={ps6}). Defaulted to ps5.")
    return ("ps5", ps5, warnings)


def build_setting_id(row: pd.Series) -> str:
    # Keep exact textual values from bigtable (no extra rounding)
    hmsP = "hmsP" + signed_token(row["hms_p_raw"])
    shmsP = "shmsP" + signed_token(row["shms_p_raw"])
    hmsTh = "hmsTh" + token_from_numeric_string(row["hms_th_raw"])
    shmsTh = "shmsTh" + token_from_numeric_string(row["shms_th_raw"])
    thpq = "thpq" + token_from_numeric_string(row["thpq_raw"])
    return "_".join([hmsP, shmsP, hmsTh, shmsTh, thpq])


def lepton_species(hms_p: float) -> str:
    return "POSITRON" if (math.isfinite(hms_p) and hms_p > 0) else "ELECTRON"


def hadron_charge(shms_p: float) -> str:
    return "PI+" if (math.isfinite(shms_p) and shms_p > 0) else "PI-"


def write_runlist(path: Path, runs: List[int]) -> None:
    path.write_text("\n".join(str(r) for r in sorted(set(runs))) + ("\n" if runs else ""))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bigtable", required=True, help="Path to rsidis_bigtable_pass0.csv")
    ap.add_argument("--out", required=True, help="Output directory (e.g., rsidis_xs_v2/settings)")
    ap.add_argument("--targets", default="LH2,LD2", help="Physics targets to generate (comma-separated)")
    ap.add_argument("--run-types", default="", help="Optional run_type filter (comma-separated)")
    ap.add_argument("--include-unknown-pass", action="store_true", help="Also generate passUNK folders")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing manifest/runlists")
    ap.add_argument("--boil-col", default="boil_corr",
                    help="Column name for boil correction. If missing, stored as NaN.")
    args = ap.parse_args()

    bigtable_path = Path(args.bigtable).expanduser().resolve()
    out_root = Path(args.out).expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    targets_keep = [norm_str(t) for t in args.targets.split(",") if t.strip()]
    run_types_keep = [norm_str(t) for t in args.run_types.split(",") if t.strip()]

    df = pd.read_csv(bigtable_path, dtype=str)

    # Required columns from current pass0 bigtable
    required_cols = [
        "run", "ebeam", "target", "run_type",
        "x", "Q2", "z", "thpq",
        "hms_p", "hms_th", "shms_p", "shms_th",
        "BCM4A_Q", "h_esing_Eff", "ps5", "ps6",
        "comp_livetime", "fan_mean", "electr_deadtime",
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing required columns in bigtable: {missing}")

    boil_col = args.boil_col
    has_boil = boil_col in df.columns

    # Helpers
    df["target_norm"] = df["target"].map(norm_str)
    df["run_type_norm"] = df["run_type"].map(norm_str)

    df["ebeam_f"] = df["ebeam"].map(safe_float)
    df["x_f"] = df["x"].map(safe_float)
    df["Q2_f"] = df["Q2"].map(safe_float)
    df["z_f"] = df["z"].map(safe_float)
    df["hms_p_f"] = df["hms_p"].map(safe_float)
    df["shms_p_f"] = df["shms_p"].map(safe_float)

    for c in ["hms_p", "hms_th", "shms_p", "shms_th", "thpq", "x", "Q2", "z"]:
        df[c + "_raw"] = df[c].astype(str).str.strip()

    # Exclude placeholders
    df = df[(df["z_f"] != -999) & (df["x_f"] != -999) & (df["Q2_f"] != -999)]

    if run_types_keep:
        df = df[df["run_type_norm"].isin(run_types_keep)]

    df_phys = df[df["target_norm"].isin(targets_keep)].copy()
    df_dummy = df[df["target_norm"] == "DUMMY"].copy()

    df_phys["pass_key"] = df_phys["ebeam_f"].map(lambda e: beam_pass_key(e)[0])
    df_phys["ebeam_key"] = df_phys["ebeam_f"].map(lambda e: beam_pass_key(e)[1])

    if not args.include_unknown_pass:
        df_phys = df_phys[df_phys["pass_key"].isin(["pass4", "pass5"])]

    df_phys["z_dir"] = df_phys["z_raw"].map(lambda s: "z" + token_from_numeric_string(s))
    df_phys["x_dir"] = df_phys["x_raw"].map(lambda s: "x" + token_from_numeric_string(s))
    df_phys["q2_dir"] = df_phys["Q2_raw"].map(lambda s: "Q2_" + token_from_numeric_string(s))

    df_phys["target_dir"] = df_phys["target_norm"]
    df_phys["run_type_dir"] = df_phys["run_type_norm"]

    df_phys["setting_id"] = df_phys.apply(build_setting_id, axis=1)
    df_phys["lepton_species"] = df_phys["hms_p_f"].map(lepton_species)
    df_phys["hadron_charge"] = df_phys["shms_p_f"].map(hadron_charge)

    # Pairing keys (keep for future automation)
    # abs(hms_p) uses 3 decimals because bigtable is mostly 3 decimals
    df_phys["abs_hms_p_3"] = df_phys["hms_p_f"].abs().round(3)
    df_phys["abs_shms_p_3"] = df_phys["shms_p_f"].abs().round(3)

    # Dummy lookup
    df_dummy["run_type_norm"] = df_dummy["run_type"].map(norm_str)
    df_dummy["ebeam_f"] = df_dummy["ebeam"].map(safe_float)
    df_dummy["pass_key"] = df_dummy["ebeam_f"].map(lambda e: beam_pass_key(e)[0])
    df_dummy["ebeam_key"] = df_dummy["ebeam_f"].map(lambda e: beam_pass_key(e)[1])

    for c in ["hms_p", "hms_th", "shms_p", "shms_th", "thpq", "x", "Q2", "z"]:
        df_dummy[c + "_raw"] = df_dummy[c].astype(str).str.strip()

    df_dummy["z_dir"] = df_dummy["z_raw"].map(lambda s: "z" + token_from_numeric_string(s))
    df_dummy["x_dir"] = df_dummy["x_raw"].map(lambda s: "x" + token_from_numeric_string(s))
    df_dummy["q2_dir"] = df_dummy["Q2_raw"].map(lambda s: "Q2_" + token_from_numeric_string(s))
    df_dummy["run_type_dir"] = df_dummy["run_type_norm"]
    df_dummy["setting_id"] = df_dummy.apply(build_setting_id, axis=1)

    dummy_key_cols = ["pass_key", "ebeam_key", "run_type_dir", "z_dir", "x_dir", "q2_dir", "setting_id"]
    dummy_lookup: Dict[Tuple[str, ...], List[int]] = {}
    for k, g in df_dummy.groupby(dummy_key_cols, dropna=False):
        dummy_lookup[tuple(k)] = [int(x) for x in g["run"].tolist() if str(x).strip().isdigit()]

    # Positron partner lookup
    partner_key_cols = [
        "pass_key", "ebeam_key", "run_type_dir", "target_dir",
        "z_dir", "x_dir", "q2_dir",
        "shms_p_raw", "hms_th_raw", "shms_th_raw", "thpq_raw",
        "abs_hms_p_3",
    ]
    partner_lookup: Dict[Tuple[str, ...], Dict[str, List[int]]] = {}
    for k, g in df_phys.groupby(partner_key_cols, dropna=False):
        key = tuple(str(x) for x in k)
        e_runs = [int(r) for r in g.loc[g["hms_p_f"] < 0, "run"].tolist() if str(r).strip().isdigit()]
        p_runs = [int(r) for r in g.loc[g["hms_p_f"] > 0, "run"].tolist() if str(r).strip().isdigit()]
        partner_lookup[key] = {"electron": e_runs, "positron": p_runs}

    group_cols = ["pass_key", "ebeam_key", "run_type_dir", "target_dir", "z_dir", "x_dir", "q2_dir", "setting_id"]
    created = 0

    for k, g in df_phys.groupby(group_cols, dropna=False):
        pass_key, ebeam_key, run_type_dir, target_dir, z_dir, x_dir, q2_dir, setting_id = k

        pass_dir = f"{pass_key}_{ebeam_key}"
        setting_path = out_root / pass_dir / run_type_dir / target_dir / z_dir / x_dir / q2_dir / setting_id
        setting_path.mkdir(parents=True, exist_ok=True)

        manifest_path = setting_path / "manifest.json"
        if manifest_path.exists() and not args.overwrite:
            continue

        warnings: List[str] = []
        if not has_boil:
            warnings.append(f"Boil correction column '{boil_col}' not found in bigtable. Stored as NaN.")

        signal_runs = [int(r) for r in g["run"].tolist() if str(r).strip().isdigit()]

        dummy_key = (pass_key, ebeam_key, run_type_dir, z_dir, x_dir, q2_dir, setting_id)
        dummy_runs = dummy_lookup.get(tuple(dummy_key), [])

        # Positron runs only for electron settings
        hms_p_val = safe_float(g["hms_p_raw"].iloc[0])
        is_electron_setting = math.isfinite(hms_p_val) and hms_p_val < 0

        partner_key = (
            pass_key, ebeam_key, run_type_dir, target_dir,
            z_dir, x_dir, q2_dir,
            g["shms_p_raw"].iloc[0], g["hms_th_raw"].iloc[0], g["shms_th_raw"].iloc[0], g["thpq_raw"].iloc[0],
            f"{float(g['abs_hms_p_3'].iloc[0]):.3f}",
        )
        partner_key = tuple(str(x) for x in partner_key)
        partner = partner_lookup.get(partner_key, {"electron": [], "positron": []})
        positron_runs = partner["positron"] if is_electron_setting else []

        write_runlist(setting_path / "runs_signal.txt", signal_runs)
        write_runlist(setting_path / "runs_dummy.txt", dummy_runs)
        write_runlist(setting_path / "runs_positron.txt", positron_runs)

        # Per-run metadata
        df_indexed = df.set_index("run", drop=False)

        rows_out = []
        for category, run_list in [("signal", signal_runs), ("dummy", dummy_runs), ("positron", positron_runs)]:
            for run in run_list:
                rr = df_indexed.loc[str(run)]
                if isinstance(rr, pd.DataFrame):
                    rr = rr.iloc[0]

                ps5 = safe_float(rr.get("ps5", "nan"))
                ps6 = safe_float(rr.get("ps6", "nan"))
                ps_choice, ps_factor, ps_warnings = choose_ps_factor(ps5, ps6)
                for w in ps_warnings:
                    warnings.append(f"Run {run}: {w}")

                boil_val = safe_float(rr.get(boil_col, "nan")) if has_boil else float("nan")

                rows_out.append({
                    "run": int(run),
                    "category": category,
                    "target_raw": str(rr.get("target", "")).strip(),
                    "run_type_raw": str(rr.get("run_type", "")).strip(),
                    "BCM4A_Q": safe_float(rr.get("BCM4A_Q", "nan")),
                    "h_esing_Eff": safe_float(rr.get("h_esing_Eff", "nan")),
                    "comp_livetime": safe_float(rr.get("comp_livetime", "nan")),
                    "electr_deadtime": safe_float(rr.get("electr_deadtime", "nan")),
                    "fan_mean": safe_float(rr.get("fan_mean", "nan")),
                    "boil_corr": boil_val,
                    "ps5": ps5,
                    "ps6": ps6,
                    "ps_choice": ps_choice,
                    "ps_factor": ps_factor,
                })

        meta_path = setting_path / "run_metadata.csv"
        with meta_path.open("w", newline="") as f:
            if rows_out:
                writer = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
                writer.writeheader()
                writer.writerows(rows_out)
            else:
                f.write("run,category\n")

        selection_text = (
            f"pass={pass_key}, ebeam={ebeam_key}, run_type={run_type_dir}, target={target_dir}, "
            f"{z_dir}, {x_dir}, {q2_dir}, setting_id={setting_id}\n"
            f"Case-insensitive matching used for categorical fields. Placeholders (-999) excluded.\n"
        )

        manifest = {
            "schema_version": SCHEMA_VERSION,
            "generated_utc": pd.Timestamp.utcnow().isoformat(),
            "bigtable": {
                "path": str(bigtable_path),
                "filename": bigtable_path.name,
                "boil_corr_column": boil_col,
            },
            "setting": {
                "pass": pass_key,
                "ebeam_key": ebeam_key,
                "run_type": run_type_dir,
                "target": target_dir,
                "z": z_dir,
                "x": x_dir,
                "Q2": q2_dir,
                "setting_id": setting_id,
                "path": str(setting_path),
            },
            "pairing": {
                "lepton_species": g["lepton_species"].iloc[0],
                "hadron_charge": g["hadron_charge"].iloc[0],
                "abs_hms_p_3": float(g["abs_hms_p_3"].iloc[0]),
                "abs_shms_p_3": float(g["abs_shms_p_3"].iloc[0]),
            },
            "selection": {
                "text": selection_text,
                "counts": {
                    "signal": len(signal_runs),
                    "dummy": len(dummy_runs),
                    "positron": len(positron_runs),
                },
            },
            "runs": {
                "signal": sorted(signal_runs),
                "dummy": sorted(dummy_runs),
                "positron": sorted(positron_runs),
            },
            "files": {
                "runs_signal": "runs_signal.txt",
                "runs_dummy": "runs_dummy.txt",
                "runs_positron": "runs_positron.txt",
                "run_metadata": "run_metadata.csv",
                "readme": "README.txt",
            },
            "corrections": {
                "available": [
                    "comp_livetime",
                    "h_esing_Eff",
                    "ps_factor",
                    "fan_mean",
                    "boil_corr",
                    "electr_deadtime",
                ],
            },
            "warnings": warnings,
        }

        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

        lines: List[str] = []
        lines.append("RSIDIS setting directory (auto-generated)")
        lines.append("")
        lines.append("Selection:")
        lines.append(selection_text.rstrip())
        lines.append("")
        lines.append("Counts:")
        lines.append(f"  signal runs:   {len(signal_runs)}")
        lines.append(f"  dummy runs:    {len(dummy_runs)}")
        lines.append(f"  positron runs: {len(positron_runs)}")
        lines.append("")
        lines.append("Runlists:")
        lines.append("  runs_signal.txt   : physics runs for this target/setting")
        lines.append("  runs_dummy.txt    : matched Dummy runs (same setting)")
        lines.append("  runs_positron.txt : matched positron runs for subtraction (hms_p > 0)")
        lines.append("")
        lines.append("Metadata:")
        lines.append("  run_metadata.csv  : per-run scalers and correction columns (boil_corr may be NaN if absent)")
        lines.append("")
        if warnings:
            lines.append("Warnings:")
            for w in warnings[:80]:
                lines.append(f"  - {w}")
            if len(warnings) > 80:
                lines.append(f"  ... ({len(warnings) - 80} more)")
            lines.append("")
        (setting_path / "README.txt").write_text("\n".join(lines) + "\n")

        created += 1

    print(f"Done. Created/updated {created} setting directories under: {out_root}")


if __name__ == "__main__":
    main()
