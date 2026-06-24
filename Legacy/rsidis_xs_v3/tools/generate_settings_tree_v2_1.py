#!/usr/bin/env python3
"""
generate_settings_tree_v2_1.py

Run command example:
  python3 generate_settings_tree_v2_1.py \
    --bigtable /path/to/rsidis_bigtable_pass0.csv \
    --outdir   /path/to/rsidis_xs_v2/settings_tree \
    --targets  LH2 LD2

Purpose:
  Generate a directory tree:
    Pass -> RunType -> Target -> z -> x -> Q2 -> SettingID

  For each electron setting (hms_p < 0) on LH2 or LD2:
    - runs_signal.txt
    - runs_dummy.txt
    - runs_positron.txt
    - runs_positron_dummy.txt
    - run_metadata.csv
    - manifest.json
    - README.txt

Physics subtraction supported by these runlists:
  (Data - PositronSignal) - (Dummy - PositronDummy)

Notes:
  - Positron identification uses HMS sign: hms_p > 0.
  - Dummy target is identified by target == "dummy" (case-insensitive).
  - run_metadata.csv schema matches:
      run,category,target_raw,run_type_raw,BCM4A_Q,h_esing_Eff,comp_livetime,
      electr_deadtime,fan_mean,boil_corr,ps5,ps6,ps_choice,ps_factor
    boil_corr is filled as 1.0 if not present in the bigtable.
"""

import re
import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


def norm_str(x: str) -> str:
    return str(x).strip()


def norm_lower(x: str) -> str:
    return norm_str(x).lower()


def safe_float(x) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def pass_key_from_ebeam(ebeam: float) -> str:
    if not pd.notna(ebeam):
        return "passUnknown"
    if abs(ebeam - 8.5831) < 1e-3:
        return "pass4"
    if abs(ebeam - 10.6716) < 1e-3:
        return "pass5"
    return f"passE{format_dir_value('', ebeam, digits=4)}".replace("passE", "passE")


def format_dir_value(prefix: str, val: float, digits: int = None) -> str:
    if not pd.notna(val):
        return prefix + "nan"
    if digits is None:
        s = str(val)
    else:
        s = f"{val:.{digits}f}"
        s = s.rstrip("0").rstrip(".") if "." in s else s
    s = s.replace("-", "neg")
    s = s.replace(".", "p")
    return f"{prefix}{s}"


def signed_mom_label(prefix: str, val: float, digits: int = 3) -> str:
    if not pd.notna(val):
        return prefix + "nan"
    sign = "neg" if val < 0 else ""
    mag = abs(val)
    s = f"{mag:.{digits}f}".replace(".", "p")
    return f"{prefix}{sign}{s}"


def choose_ps(ps5: float, ps6: float, tol: float = 1e-6) -> Tuple[str, float]:
    def is_one(x: float) -> bool:
        return pd.notna(x) and abs(x - 1.0) < tol

    if is_one(ps5) and not is_one(ps6):
        return ("ps5", float(ps5))
    if is_one(ps6) and not is_one(ps5):
        return ("ps6", float(ps6))
    if is_one(ps5) and is_one(ps6):
        return ("ps5", float(ps5))
    return ("ps5", float(ps5) if pd.notna(ps5) else 1.0)


def build_setting_id(row: pd.Series) -> str:
    hms_p = safe_float(row["hms_p"])
    shms_p = safe_float(row["shms_p"])
    hms_th = safe_float(row["hms_th"])
    shms_th = safe_float(row["shms_th"])
    thpq = safe_float(row["thpq"])

    parts = [
        signed_mom_label("hmsP", hms_p, 3),
        format_dir_value("shmsP", shms_p, 3),
        format_dir_value("hmsTh", hms_th, 3),
        format_dir_value("shmsTh", shms_th, 3),
        format_dir_value("thpq", thpq, 3),
    ]
    return "_".join(parts)


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def write_runlist(path: Path, runs: List[int]) -> None:
    runs_sorted = sorted(set(int(r) for r in runs))
    path.write_text("\n".join(str(r) for r in runs_sorted) + ("\n" if runs_sorted else ""))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bigtable", required=True, help="Path to rsidis_bigtable_pass0.csv")
    ap.add_argument("--outdir", required=True, help="Output settings_tree directory")
    ap.add_argument("--targets", nargs="+", default=["LH2", "LD2"], help="Targets to generate settings for (default LH2 LD2)")
    ap.add_argument("--overwrite", action="store_true", help="Allow overwriting existing files")
    args = ap.parse_args()

    bigtable_path = Path(args.bigtable).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()

    df = pd.read_csv(bigtable_path)

    required = ["run", "ebeam", "target", "run_type", "x", "Q2", "z", "hms_p", "hms_th", "shms_p", "shms_th", "thpq",
                "BCM4A_Q", "h_esing_Eff", "comp_livetime", "electr_deadtime", "fan_mean", "ps5", "ps6"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing required columns in bigtable: {missing}")

    # Fill boil_corr if missing in bigtable
    if "boil_corr" not in df.columns:
        df["boil_corr"] = 1.0

    df["target_raw"] = df["target"].map(norm_str)
    df["target_norm"] = df["target"].map(norm_lower)

    df["run_type_raw"] = df["run_type"].map(norm_str)
    df["run_type_norm"] = df["run_type"].map(norm_lower)

    df["ebeam_f"] = df["ebeam"].map(safe_float)
    df["pass_key"] = df["ebeam_f"].map(pass_key_from_ebeam)

    df["x_f"] = df["x"].map(safe_float)
    df["q2_f"] = df["Q2"].map(safe_float)
    df["z_f"] = df["z"].map(safe_float)

    df["hms_p_f"] = df["hms_p"].map(safe_float)
    df["shms_p_f"] = df["shms_p"].map(safe_float)

    df["abs_hms_p_3"] = df["hms_p_f"].abs().round(3)

    df["z_dir"] = df["z_f"].map(lambda v: format_dir_value("z", v, digits=2))
    df["x_dir"] = df["x_f"].map(lambda v: format_dir_value("x", v, digits=2))
    df["q2_dir"] = df["q2_f"].map(lambda v: format_dir_value("q2", v, digits=2))

    df["run_type_dir"] = df["run_type_norm"].map(lambda s: re.sub(r"[^a-z0-9\+\-]+", "", s))
    df["setting_id"] = df.apply(build_setting_id, axis=1)

    targets_keep = set(t.lower() for t in args.targets)
    df_phys = df[df["target_norm"].isin(targets_keep)].copy()
    df_dummy = df[df["target_norm"] == "dummy"].copy()

    # Partner lookup for physics targets (LH2/LD2), keyed by abs(hms_p) and other kinematics, excluding hms sign.
    partner_key_cols = [
        "pass_key", "run_type_dir", "target_norm",
        "z_dir", "x_dir", "q2_dir",
        "shms_p", "hms_th", "shms_th", "thpq",
        "abs_hms_p_3",
    ]

    partner_lookup: Dict[Tuple[str, ...], Dict[str, List[int]]] = {}
    for k, g in df_phys.groupby(partner_key_cols, dropna=False):
        key = tuple(str(x) for x in k)
        e_runs = [int(r) for r in g.loc[g["hms_p_f"] < 0, "run"].tolist()]
        p_runs = [int(r) for r in g.loc[g["hms_p_f"] > 0, "run"].tolist()]
        partner_lookup[key] = {"electron": e_runs, "positron": p_runs}

    # Partner lookup for dummy target
    dummy_partner_key_cols = [
        "pass_key", "run_type_dir",
        "z_dir", "x_dir", "q2_dir",
        "shms_p", "hms_th", "shms_th", "thpq",
        "abs_hms_p_3",
    ]

    dummy_partner_lookup: Dict[Tuple[str, ...], Dict[str, List[int]]] = {}
    for k, g in df_dummy.groupby(dummy_partner_key_cols, dropna=False):
        key = tuple(str(x) for x in k)
        e_runs = [int(r) for r in g.loc[g["hms_p_f"] < 0, "run"].tolist()]
        p_runs = [int(r) for r in g.loc[g["hms_p_f"] > 0, "run"].tolist()]
        dummy_partner_lookup[key] = {"electron": e_runs, "positron": p_runs}

    ensure_dir(outdir)

    # Generate settings only for electron configurations in physics targets
    df_e_settings = df_phys[df_phys["hms_p_f"] < 0].copy()

    group_cols = ["pass_key", "run_type_dir", "target_norm", "z_dir", "x_dir", "q2_dir", "setting_id"]
    n_settings = 0

    for k, g in df_e_settings.groupby(group_cols, dropna=False):
        pass_key, run_type_dir, target_norm, z_dir, x_dir, q2_dir, setting_id = [str(x) for x in k]

        # Primary electron runs for this setting
        signal_runs = [int(r) for r in g["run"].tolist()]

        # Find positron partner runs on the same physics target
        partner_key = (
            pass_key, run_type_dir, target_norm,
            z_dir, x_dir, q2_dir,
            str(g["shms_p"].iloc[0]), str(g["hms_th"].iloc[0]), str(g["shms_th"].iloc[0]), str(g["thpq"].iloc[0]),
            str(float(g["abs_hms_p_3"].iloc[0])),
        )
        partner_key = tuple(str(x) for x in partner_key)
        partner = partner_lookup.get(partner_key, {"electron": [], "positron": []})
        positron_runs = partner["positron"]

        # Find dummy electron and dummy positron partner runs
        dummy_key = (
            pass_key, run_type_dir,
            z_dir, x_dir, q2_dir,
            str(g["shms_p"].iloc[0]), str(g["hms_th"].iloc[0]), str(g["shms_th"].iloc[0]), str(g["thpq"].iloc[0]),
            str(float(g["abs_hms_p_3"].iloc[0])),
        )
        dummy_key = tuple(str(x) for x in dummy_key)
        dummy_partner = dummy_partner_lookup.get(dummy_key, {"electron": [], "positron": []})
        dummy_runs = dummy_partner["electron"]
        positron_dummy_runs = dummy_partner["positron"]

        # Build directory path
        target_dir = target_norm.upper()
        setting_path = outdir / pass_key / run_type_dir / target_dir / z_dir / x_dir / q2_dir / setting_id
        ensure_dir(setting_path)

        # Write runlists
        write_runlist(setting_path / "runs_signal.txt", signal_runs)
        write_runlist(setting_path / "runs_dummy.txt", dummy_runs)
        write_runlist(setting_path / "runs_positron.txt", positron_runs)
        write_runlist(setting_path / "runs_positron_dummy.txt", positron_dummy_runs)

        # Build run_metadata.csv for union of all runs in this setting
        all_runs = []
        for cat, lst in [("signal", signal_runs), ("dummy", dummy_runs), ("positron", positron_runs), ("positron_dummy", positron_dummy_runs)]:
            for r in lst:
                all_runs.append((int(r), cat))

        rows = []
        df_index = df.set_index("run", drop=False)
        for run, cat in all_runs:
            if run not in df_index.index:
                continue
            row = df_index.loc[run]
            ps_choice, ps_factor = choose_ps(safe_float(row["ps5"]), safe_float(row["ps6"]))
            rows.append({
                "run": int(run),
                "category": cat,
                "target_raw": norm_str(row["target_raw"]),
                "run_type_raw": norm_str(row["run_type_raw"]),
                "BCM4A_Q": safe_float(row["BCM4A_Q"]),
                "h_esing_Eff": safe_float(row["h_esing_Eff"]),
                "comp_livetime": safe_float(row["comp_livetime"]),
                "electr_deadtime": safe_float(row["electr_deadtime"]),
                "fan_mean": safe_float(row["fan_mean"]),
                "boil_corr": safe_float(row["boil_corr"]),
                "ps5": safe_float(row["ps5"]),
                "ps6": safe_float(row["ps6"]),
                "ps_choice": ps_choice,
                "ps_factor": ps_factor,
            })

        meta_df = pd.DataFrame(rows)
        meta_df.to_csv(setting_path / "run_metadata.csv", index=False)

        # Manifest
        manifest = {
            "schema_version": 1,
            "bigtable": str(bigtable_path),
            "path_keys": {
                "pass": pass_key,
                "run_type": run_type_dir,
                "target": target_dir,
                "z": z_dir,
                "x": x_dir,
                "q2": q2_dir,
                "setting_id": setting_id,
            },
            "pairing_keys": {
                "abs_hms_p_3": float(g["abs_hms_p_3"].iloc[0]),
                "shms_p": safe_float(g["shms_p_f"].iloc[0]),
                "run_type": norm_str(g["run_type_raw"].iloc[0]),
            },
            "selection_counts": {
                "signal": len(set(signal_runs)),
                "dummy": len(set(dummy_runs)),
                "positron": len(set(positron_runs)),
                "positron_dummy": len(set(positron_dummy_runs)),
            },
            "files": {
                "runs_signal": "runs_signal.txt",
                "runs_dummy": "runs_dummy.txt",
                "runs_positron": "runs_positron.txt",
                "runs_positron_dummy": "runs_positron_dummy.txt",
                "run_metadata": "run_metadata.csv",
            }
        }
        (setting_path / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

        # README.txt per setting
        lines = []
        lines.append("Setting directory generated by generate_settings_tree_v2_1.py")
        lines.append(f"pass: {pass_key}")
        lines.append(f"run_type: {run_type_dir}")
        lines.append(f"target: {target_dir}")
        lines.append(f"z: {z_dir}  x: {x_dir}  q2: {q2_dir}")
        lines.append(f"setting_id: {setting_id}")
        lines.append("")
        lines.append("Runlists:")
        lines.append(f"  runs_signal.txt         : {len(set(signal_runs))}")
        lines.append(f"  runs_positron.txt       : {len(set(positron_runs))}")
        lines.append(f"  runs_dummy.txt          : {len(set(dummy_runs))}")
        lines.append(f"  runs_positron_dummy.txt : {len(set(positron_dummy_runs))}")
        lines.append("")
        lines.append("Subtraction formula intended:")
        lines.append("  (Data - PositronSignal) - (Dummy - PositronDummy)")
        (setting_path / "README.txt").write_text("\n".join(lines) + "\n")

        n_settings += 1

    print(f"Generated {n_settings} settings under {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
