#!/usr/bin/env python3
"""
Plot SIDIS fit parameters from:
  results/<group_id>/tables/fit_parameters.csv

This version is robust against a common failure mode where the header (and/or rows)
get written without '\n', causing the whole file to become one giant line.

Outputs:
  results/<group_id>/PNGs/fitparam_plots/*.png

Run from rsidis_xs_v5/:
  python3 macros/plot_fit_parameters.py --group grp_LH2_zOv_x0p25_Q23p3_tpq2p0
"""

import argparse
import re
from io import StringIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for farm/VDI

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


EXPECTED_COLS = [
    "mode","group_id","curve_label","setting_id",
    "pt_bin","z_bin","pt_lo","pt_hi","pt_center","z_lo","z_hi","z_center",
    "n_points","M0","M0_err","A","A_err","B","B_err","chi2","ndf","prob"
]


def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def _looks_broken(raw: str) -> bool:
    # If the file has only 0–1 newlines, it's almost certainly concatenated.
    if raw.count("\n") <= 1:
        return True
    # Or: first line contains header AND also contains "single," etc later => header+row glued.
    first = raw.splitlines()[0] if raw.splitlines() else raw
    if "mode,group_id" in first and re.search(r"(single|group|overlay),", first):
        return True
    return False


def repair_fit_csv_text(raw: str) -> str:
    """
    Repairs:
      1) Missing newline after header
      2) Missing newlines between rows (rows start with mode = single/group/overlay)
    """
    s = raw

    # Normalize newlines
    s = s.replace("\r\n", "\n").replace("\r", "\n")

    # Ensure a newline after the header (header ends with 'prob')
    # Find the first occurrence of the exact header prefix.
    hdr_pos = s.find("mode,group_id")
    if hdr_pos != -1:
        prob_pos = s.find("prob", hdr_pos)
        if prob_pos != -1:
            hdr_end = prob_pos + len("prob")
            if hdr_end < len(s) and s[hdr_end] != "\n":
                s = s[:hdr_end] + "\n" + s[hdr_end:]

    # Now ensure each row begins on a new line.
    # Rows begin with mode field: 'single,' or 'group,' or 'overlay,'
    # Insert newline BEFORE these tokens if not already preceded by '\n' and not at start.
    s = re.sub(r"(?<!\n)(?=(?:single|group|overlay),)", "\n", s)

    return s


def read_fit_parameters_csv(path: Path) -> pd.DataFrame:
    raw = path.read_text(errors="replace")

    if _looks_broken(raw):
        raw2 = repair_fit_csv_text(raw)
        df = pd.read_csv(StringIO(raw2))
    else:
        df = pd.read_csv(path)

    # If columns are still wrong, do a last-resort repair attempt
    if len(df.columns) < 5 or ("mode" not in df.columns and "curve_label" not in df.columns):
        raw2 = repair_fit_csv_text(raw)
        df = pd.read_csv(StringIO(raw2))

    # Normalize / ensure expected columns exist (some older files may omit z_center)
    for c in EXPECTED_COLS:
        if c not in df.columns:
            df[c] = np.nan

    # Coerce numeric columns
    num_cols = [
        "pt_bin","z_bin","pt_lo","pt_hi","pt_center","z_lo","z_hi","z_center",
        "n_points","M0","M0_err","A","A_err","B","B_err","chi2","ndf","prob"
    ]
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # If z_center missing, derive from curve_label like z0p36
    if df["z_center"].isna().all():
        def z_from_label(s):
            s = str(s).strip()
            if s.startswith("z"):
                s = s[1:]
            try:
                return float(s.replace("p", "."))
            except Exception:
                return np.nan
        df["z_center"] = df["curve_label"].map(z_from_label)

    # Drop totally broken rows
    df = df.dropna(subset=["curve_label","pt_bin","pt_center","A","B","M0"], how="any")

    return df


def plot_vs_pt(df, ycol, yerr, outpath, title, ylabel):
    fig, ax = plt.subplots(figsize=(7, 5))

    for lab, g in df.groupby("curve_label", sort=True):
        g = g.sort_values("pt_center")
        x = g["pt_center"].to_numpy()
        y = g[ycol].to_numpy()
        ye = g[yerr].to_numpy() if yerr in g.columns else None
        ax.errorbar(x, y, yerr=ye, fmt="o", capsize=2, label=str(lab))

    ax.set_xlabel(r"$p_T$ (GeV)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def plot_prob_vs_pt(df, outpath):
    if "prob" not in df.columns:
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    for lab, g in df.groupby("curve_label", sort=True):
        g = g.sort_values("pt_center")
        ax.plot(g["pt_center"].to_numpy(), g["prob"].to_numpy(), marker="o", linestyle="-", label=str(lab))
    ax.set_xlabel(r"$p_T$ (GeV)")
    ax.set_ylabel("Fit probability")
    ax.set_title("Fit quality: probability vs $p_T$")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def plot_vs_z_grid(df, ycol, yerr, outpath, suptitle, ylabel):
    pt_bins = sorted(df["pt_bin"].dropna().unique().tolist())
    if len(pt_bins) == 0:
        return

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharey=True)
    axes = axes.ravel()

    for i, ptb in enumerate(pt_bins[:4]):
        ax = axes[i]
        gd = df[df["pt_bin"] == ptb].copy()

        gd = gd.sort_values("z_center")
        # One point per curve_label (z setting) in this pt bin:
        for lab, g in gd.groupby("curve_label", sort=True):
            x = g["z_center"].to_numpy()
            y = g[ycol].to_numpy()
            ye = g[yerr].to_numpy() if yerr in g.columns else None
            ax.errorbar(x, y, yerr=ye, fmt="o", capsize=2, label=str(lab))

        # pt edges if available
        if "pt_lo" in gd.columns and "pt_hi" in gd.columns and len(gd) > 0:
            ax.set_title(fr"$p_T \in$ [{gd['pt_lo'].iloc[0]:.2f}, {gd['pt_hi'].iloc[0]:.2f}] GeV")
        else:
            ax.set_title(f"pt_bin={int(ptb)}")

        ax.set_xlabel("z")
        ax.grid(True, alpha=0.3)
        if i == 0:
            ax.legend(frameon=False)

    # hide unused pads
    for j in range(min(4, len(pt_bins)), 4):
        axes[j].axis("off")

    axes[0].set_ylabel(ylabel)
    axes[2].set_ylabel(ylabel)

    fig.suptitle(suptitle)
    fig.tight_layout(rect=[0, 0.0, 1, 0.95])
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    ap.add_argument("--results-dir", default="results")
    args = ap.parse_args()

    results_dir = Path(args.results_dir).resolve()
    csv_path = results_dir / args.group / "tables" / "fit_parameters.csv"
    if not csv_path.exists():
        raise SystemExit(f"ERROR: not found: {csv_path}")

    df = read_fit_parameters_csv(csv_path)

    outdir = results_dir / args.group / "PNGs" / "fitparam_plots"
    ensure_dir(outdir)

    plot_vs_pt(df, "A", "A_err", outdir / "A_vs_pT.png", "Fit parameter A vs $p_T$", "A")
    plot_vs_pt(df, "B", "B_err", outdir / "B_vs_pT.png", "Fit parameter B vs $p_T$", "B")
    plot_vs_pt(df, "M0", "M0_err", outdir / "M0_vs_pT.png", "Fit parameter $M_0$ vs $p_T$", r"$M_0$")
    plot_prob_vs_pt(df, outdir / "prob_vs_pT.png")

    plot_vs_z_grid(df, "A", "A_err", outdir / "A_vs_z_grid.png", "A vs z (pads: $p_T$ bins)", "A")
    plot_vs_z_grid(df, "B", "B_err", outdir / "B_vs_z_grid.png", "B vs z (pads: $p_T$ bins)", "B")

    print(f"Wrote plots to: {outdir}")
    # Small debug summary
    print("Rows read:", len(df))
    print("curve_label:", sorted(df["curve_label"].unique().tolist()))
    print("pt_bin:", sorted(df["pt_bin"].unique().tolist()))


if __name__ == "__main__":
    main()
