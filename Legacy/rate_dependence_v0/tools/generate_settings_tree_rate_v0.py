#!/usr/bin/env python3
"""
generate_settings_tree_rate_v0.py  (FIXED: v5-style setting_id + target-clean partner runlists)

This version fixes two issues you reported:

1) setting_id sign:
   - Settings are seeded from ELECTRON runs (HMS p < 0) on physics targets (LH2/LD2),
     and the setting_id is built using the SIGNED HMS momentum, so it becomes:
       hmsPneg1p531_...
     (exactly like rsidis_xs_v5)

2) positron runlist contamination:
   - runs_positron.txt is now restricted to the SAME physics target as the seed setting
     (so LH2 settings only pick LH2 positron runs, etc.).
   - Dummy partners come only from target == Dummy (case-insensitive).

Schema mirrors rsidis_xs_v5:
  settings/<pass>/<run_type>/<target>/<z>/<x>/<Q2>/<setting_id>/
    runs_signal.txt
    runs_positron.txt
    runs_dummy.txt
    runs_positron_dummy.txt
    run_metadata.csv
    manifest.json
    README.txt

Rate-dependence columns written to run_metadata.csv:
  BCM2_I, BCM2_Q, comp_livetime, h_esing_Eff, boil_corr, ps5/ps6 chooser, normyield fields.

Run-type:
  bigtable run_type like 'PI+SIDIS'/'PI-SIDIS' is normalized to 'pi+sidis'/'pi-sidis'
  --run-types coin expands to both pi+sidis and pi-sidis

Usage:
  python3 tools/generate_settings_tree_rate_v0.py \
    --bigtable bigtable/rsidis_bigtable_pass0.csv \
    --outdir settings \
    --targets LH2 LD2 \
    --run-types coin
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


def norm_str(x) -> str:
    return str(x).strip()


def norm_lower(x) -> str:
    return norm_str(x).lower()


def safe_float(x) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def write_runlist(path: Path, runs: List[int]) -> None:
    runs_u = sorted(set(int(r) for r in runs))
    path.write_text("\n".join(str(r) for r in runs_u) + ("\n" if runs_u else ""))


def fmt_dir(prefix: str, val: float, digits: int = 2) -> str:
    if not pd.notna(val):
        return f"{prefix}nan"
    s = f"{val:.{digits}f}"
    s = s.rstrip("0").rstrip(".") if "." in s else s
    s = s.replace("-", "neg").replace(".", "p")
    return f"{prefix}{s}"


def fmt_mom_signed(prefix: str, val: float, digits: int = 3) -> str:
    """Like v5: include 'neg' when val<0, otherwise no sign tag."""
    if not pd.notna(val):
        return f"{prefix}nan"
    sign = "neg" if val < 0 else ""
    mag = abs(val)
    s = f"{mag:.{digits}f}".replace(".", "p")
    return f"{prefix}{sign}{s}"


def pass_key_from_ebeam(ebeam: float) -> str:
    if not pd.notna(ebeam):
        return "passUnknown"
    if abs(ebeam - 8.5831) < 1e-3:
        return "pass4"
    if abs(ebeam - 10.6716) < 1e-3:
        return "pass5"
    return "passUnknown"


def normalize_run_type(s: str) -> str:
    return norm_lower(s).replace(" ", "")


def choose_ps(ps5: float, ps6: float, tol: float = 1e-6) -> Tuple[str, float]:
    def is_one(x: float) -> bool:
        return pd.notna(x) and abs(x - 1.0) < tol

    if is_one(ps5) and not is_one(ps6):
        return ("ps5", float(ps5))
    if is_one(ps6) and not is_one(ps5):
        return ("ps6", float(ps6))
    if is_one(ps5) and is_one(ps6):
        return ("ps5", float(ps5))
    return ("ps5", float(ps5) if pd.notna(ps5) and ps5 > 0 else 1.0)


def build_setting_id_v5_style(row: pd.Series) -> str:
    """Setting id should be *electron signed* HMS momentum (neg), like rsidis_xs_v5."""
    hms_p  = safe_float(row.get("hms_p", float("nan")))   # signed
    shms_p = safe_float(row.get("shms_p", float("nan")))
    hms_th = safe_float(row.get("hms_th", float("nan")))
    shms_th= safe_float(row.get("shms_th", float("nan")))
    thpq   = safe_float(row.get("thpq", float("nan")))

    parts = [
        fmt_mom_signed("hmsP", hms_p, digits=3),
        fmt_dir("shmsP", shms_p, digits=3),
        fmt_dir("hmsTh", hms_th, digits=3),
        fmt_dir("shmsTh", shms_th, digits=3),
        fmt_dir("thpq", thpq, digits=3),
    ]
    return "_".join(parts)


def make_partner_key(row: pd.Series) -> str:
    """Partner key uses abs(HMS p) so electron/positron match, plus other kinematics."""
    hms_p_abs = abs(safe_float(row.get("hms_p", float("nan"))))
    shms_p    = safe_float(row.get("shms_p", float("nan")))
    hms_th    = safe_float(row.get("hms_th", float("nan")))
    shms_th   = safe_float(row.get("shms_th", float("nan")))
    thpq      = safe_float(row.get("thpq", float("nan")))
    z         = safe_float(row.get("z", float("nan")))
    x         = safe_float(row.get("x", float("nan")))
    Q2        = safe_float(row.get("Q2", row.get("q2", float("nan"))))
    rt        = normalize_run_type(row.get("run_type", ""))
    return f"{rt}|{z:.4f}|{x:.4f}|{Q2:.4f}|{hms_p_abs:.4f}|{shms_p:.4f}|{hms_th:.4f}|{shms_th:.4f}|{thpq:.4f}"


def build_map(sub: pd.DataFrame) -> Dict[str, List[int]]:
    m: Dict[str, List[int]] = {}
    for _, r in sub.iterrows():
        k = make_partner_key(r)
        m.setdefault(k, []).append(int(r["run"]))
    for k in list(m.keys()):
        m[k] = sorted(set(m[k]))
    return m


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bigtable", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--targets", nargs="+", default=["LH2", "LD2"], help="Physics targets to seed settings")
    ap.add_argument("--run-types", nargs="+", default=["coin"], help="Run types (use 'coin' for pi± sidis)")
    args = ap.parse_args()

    bigtable_path = Path(args.bigtable).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    ensure_dir(outdir)

    df = pd.read_csv(bigtable_path)

    required = ["run","target","run_type","hms_p","shms_p","hms_th","shms_th","thpq","z","x"]
    for req in required:
        if req not in df.columns:
            raise RuntimeError(f"bigtable missing required column: {req}")
    if "Q2" not in df.columns and "q2" not in df.columns:
        raise RuntimeError("bigtable missing required column: Q2 (or q2)")

    df["run"] = df["run"].astype(int)

    df["target_raw"] = df["target"].map(norm_str)
    df["target_norm"] = df["target_raw"].map(norm_lower)
    df["is_dummy"] = df["target_norm"].eq("dummy")

    df["run_type_raw"] = df["run_type"].map(norm_str)
    df["run_type_norm"] = df["run_type_raw"].map(normalize_run_type)

    # expand run types
    rt_in = [normalize_run_type(x) for x in args.run_types]
    run_types_keep = set()
    for x in rt_in:
        if x == "coin":
            run_types_keep |= {"pi+sidis","pi-sidis"}
        else:
            run_types_keep.add(x)

    targets_keep = set(t.lower() for t in args.targets)

    # pass key
    if "ebeam" in df.columns:
        df["pass_key"] = df["ebeam"].map(safe_float).map(pass_key_from_ebeam)
    else:
        df["pass_key"] = "passUnknown"

    # dirs
    df["z_dir"] = df["z"].map(safe_float).map(lambda v: fmt_dir("z", v, digits=2))
    df["x_dir"] = df["x"].map(safe_float).map(lambda v: fmt_dir("x", v, digits=2))
    if "Q2" in df.columns:
        df["Q2_dir"] = df["Q2"].map(safe_float).map(lambda v: fmt_dir("Q2", v, digits=2))
    else:
        df["Q2_dir"] = df["q2"].map(safe_float).map(lambda v: fmt_dir("Q2", v, digits=2))

    # filter run types now
    df = df[df["run_type_norm"].isin(run_types_keep)].copy()
    df["hms_p_f"] = df["hms_p"].map(safe_float)

    # partitions
    df_phys = df[~df["is_dummy"]].copy()
    df_dum  = df[df["is_dummy"]].copy()

    # Seed settings ONLY from electron runs on requested physics targets
    df_seed = df_phys[(df_phys["target_norm"].isin(targets_keep)) & (df_phys["hms_p_f"] < 0)].copy()

    # IMPORTANT FIX: setting_id must be v5-style (signed HMS momentum)
    df_seed["setting_id"] = df_seed.apply(build_setting_id_v5_style, axis=1)
    # Also keep setting_id on full df for manifest lookup convenience (not strictly necessary)
    df["setting_id_seed_style"] = df.apply(build_setting_id_v5_style, axis=1)

    # Partner maps:
    #  - positron map is built PER physics target, to avoid mixing LH2/LD2/C/Cu etc.
    #  - dummy maps are global dummy-only (target Dummy), split by electron/positron.
    df_phys_pos = df_phys[df_phys["hms_p_f"] > 0].copy()

    pos_map_by_tgt: Dict[str, Dict[str, List[int]]] = {}
    for tgt in targets_keep:
        sub = df_phys_pos[df_phys_pos["target_norm"] == tgt].copy()
        pos_map_by_tgt[tgt] = build_map(sub)

    df_dum["hms_p_f"] = df_dum["hms_p"].map(safe_float)
    dum_map    = build_map(df_dum[df_dum["hms_p_f"] < 0].copy())
    posdum_map = build_map(df_dum[df_dum["hms_p_f"] > 0].copy())

    df_index = df.set_index("run", drop=False)

    # Grouping columns mirror v5 structure
    group_cols = ["pass_key","run_type_norm","target_norm","z_dir","x_dir","Q2_dir","setting_id"]
    n_settings = 0

    for k, g in df_seed.groupby(group_cols, dropna=False):
        pass_key, rt, tgt_norm, z_dir, x_dir, Q2_dir, setting_id = [str(x) for x in k]
        tgt_dir = tgt_norm.upper()
        setting_path = outdir / pass_key / rt / tgt_dir / z_dir / x_dir / Q2_dir / setting_id
        ensure_dir(setting_path)

        runs_signal = sorted(set(int(r) for r in g["run"].tolist()))

        # partner lookups restricted to same physics target
        pos_map = pos_map_by_tgt.get(tgt_norm, {})
        runs_pos, runs_dum, runs_posdum = [], [], []
        for _, row in g.iterrows():
            key = make_partner_key(row)
            runs_pos    += pos_map.get(key, [])
            runs_dum    += dum_map.get(key, [])
            runs_posdum += posdum_map.get(key, [])
        runs_pos = sorted(set(runs_pos))
        runs_dum = sorted(set(runs_dum))
        runs_posdum = sorted(set(runs_posdum))

        write_runlist(setting_path/"runs_signal.txt", runs_signal)
        write_runlist(setting_path/"runs_positron.txt", runs_pos)
        write_runlist(setting_path/"runs_dummy.txt", runs_dum)
        write_runlist(setting_path/"runs_positron_dummy.txt", runs_posdum)

        # metadata for all referenced runs
        all_runs = [("signal",r) for r in runs_signal] +                    [("positron",r) for r in runs_pos] +                    [("dummy",r) for r in runs_dum] +                    [("positron_dummy",r) for r in runs_posdum]

        rows = []
        for cat, run in all_runs:
            row = df_index.loc[int(run)]
            ps_choice, ps_factor = choose_ps(safe_float(row.get("ps5", float("nan"))),
                                             safe_float(row.get("ps6", float("nan"))))
            rows.append({
                "run": int(run),
                "category": cat,
                "target_raw": norm_str(row.get("target","")),
                "run_type_raw": norm_str(row.get("run_type","")),
                "BCM2_I": safe_float(row.get("BCM2_I", float("nan"))),
                "BCM2_Q": safe_float(row.get("BCM2_Q", float("nan"))),
                "h_esing_Eff": safe_float(row.get("h_esing_Eff", float("nan"))),
                "comp_livetime": safe_float(row.get("comp_livetime", float("nan"))),
                "boil_corr": safe_float(row.get("boil_corr", 1.0)),
                "ps5": safe_float(row.get("ps5", float("nan"))),
                "ps6": safe_float(row.get("ps6", float("nan"))),
                "ps_choice": ps_choice,
                "ps_factor": float(ps_factor),
                "normyield": safe_float(row.get("normyield", float("nan"))),
                "normyield_err": safe_float(row.get("normyield_err", float("nan"))),
            })

        pd.DataFrame(rows).to_csv(setting_path/"run_metadata.csv", index=False)

        manifest = {
            "schema_version": 1,
            "bigtable": str(bigtable_path),
            "path_keys": {
                "pass": pass_key,
                "run_type": rt,
                "target": tgt_dir,
                "z": z_dir,
                "x": x_dir,
                "Q2": Q2_dir,
                "setting_id": setting_id,
            },
            "files": {
                "runs_signal": "runs_signal.txt",
                "runs_positron": "runs_positron.txt",
                "runs_dummy": "runs_dummy.txt",
                "runs_positron_dummy": "runs_positron_dummy.txt",
                "run_metadata": "run_metadata.csv",
            },
            "rootfiles": {
                "type": "skimmed",
                "tree_name": "T",
                "base_dir": "./Skimmed_ROOTfiles",
                "pattern": "skimmed_coin_replay_production_{run}_-1.root",
            }
        }
        (setting_path/"manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

        (setting_path/"README.txt").write_text(
            "\n".join([
                "Setting directory generated by generate_settings_tree_rate_v0.py (v5-style fixes applied)",
                f"pass: {pass_key}",
                f"run_type: {rt}",
                f"target: {tgt_dir}",
                f"z: {z_dir}  x: {x_dir}  Q2: {Q2_dir}",
                f"setting_id: {setting_id}",
                "",
                f"runs_signal.txt : {len(runs_signal)}",
                f"runs_positron.txt : {len(runs_pos)}",
                f"runs_dummy.txt : {len(runs_dum)}",
                f"runs_positron_dummy.txt : {len(runs_posdum)}",
                "",
                "Notes:",
                "  - Settings are seeded from electron runs on physics targets only (v5 behavior).",
                "  - setting_id uses signed HMS momentum (neg for signal), matching rsidis_xs_v5.",
                "  - positron partners are restricted to the same physics target as the seed setting.",
                "  - dummy partners come only from target == Dummy.",
            ]) + "\n"
        )

        n_settings += 1

    print(f"Generated {n_settings} settings under {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
