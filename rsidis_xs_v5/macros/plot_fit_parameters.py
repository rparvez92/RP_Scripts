#!/usr/bin/env python3
"""
Plot SIDIS fit parameters produced by PlotCoinXsec group-overlay mode.

Reads:
  results/<group_id>/tables/fit_parameters.csv

Writes PNGs to:
  results/<group_id>/PNGs/fitparam_plots/*.png

Usage (run from rsidis_xs_v5/):
  python3 macros/plot_fit_parameters.py --group grp_LH2_zOv_x0p25_Q23p3_tpq2p0

Optional filters (recommended once you go to many phi bins):
  --min-points 10 --min-prob 0.01
"""

import argparse
import math
from pathlib import Path

# Force a non-interactive backend (prevents EGL/GLX issues on farm/VDI)
import matplotlib
matplotlib.use("Agg")

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def _safe_float(x):
    try:
        return float(x)
    except Exception:
        return float("nan")


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def plot_vs_pt(df: pd.DataFrame, ycol: str, yerr: str, outpath: Path, title: str, ylabel: str,
               min_points: int = 0, min_prob: float = 0.0) -> None:
    """One figure: y vs pT with separate curves by curve_label."""
    d = df.copy()

    if "n_points" in d.columns:
        d = d[d["n_points"].astype(int) >= int(min_points)]
    if "prob" in d.columns:
        d = d[d["prob"].astype(float) >= float(min_prob)]

    fig, ax = plt.subplots(figsize=(7, 5))

    for lab, g in d.groupby("curve_label", sort=True):
        g = g.sort_values("pt_center")
        x = g["pt_center"].astype(float).to_numpy()
        y = g[ycol].astype(float).to_numpy()
        ye = g[yerr].astype(float).to_numpy() if yerr in g.columns else None
        ax.errorbar(x, y, yerr=ye, fmt='o', capsize=2, label=str(lab))

    ax.set_xlabel(r"$p_T$ (GeV)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)

    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def plot_prob_vs_pt(df: pd.DataFrame, outpath: Path, min_points: int = 0) -> None:
    """One figure: fit probability vs pT with separate curves by curve_label."""
    d = df.copy()
    if "n_points" in d.columns:
        d = d[d["n_points"].astype(int) >= int(min_points)]

    fig, ax = plt.subplots(figsize=(7, 5))

    for lab, g in d.groupby("curve_label", sort=True):
        g = g.sort_values("pt_center")
        x = g["pt_center"].astype(float).to_numpy()
        y = g["prob"].astype(float).to_numpy()
        ax.plot(x, y, marker='o', linestyle='-', label=str(lab))

    ax.set_xlabel(r"$p_T$ (GeV)")
    ax.set_ylabel("Fit probability")
    ax.set_title("Fit quality: probability vs $p_T$")
    ax.set_yscale("log")  # probability spans orders of magnitude
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False)

    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def _flatten_axes(axes):
    """Return a flat list of matplotlib Axes from plt.subplots output."""
    if isinstance(axes, np.ndarray):
        return list(axes.ravel())
    # single Axes
    return [axes]


def plot_AorB_vs_z_grid(df: pd.DataFrame, ycol: str, yerr: str, outpath: Path,
                        title: str, ylabel: str, min_points: int = 0, min_prob: float = 0.0) -> None:
    """
    One figure: 2x2 grid (pt_bin pads), y vs z for each pad.
    Points are the different curve_label entries (z settings).
    """
    d = df.copy()

    if "n_points" in d.columns:
        d = d[d["n_points"].astype(int) >= int(min_points)]
    if "prob" in d.columns:
        d = d[d["prob"].astype(float) >= float(min_prob)]

    pt_bins = sorted(d["pt_bin"].unique().tolist())
    n = len(pt_bins)
    if n == 0:
        raise RuntimeError("No rows left after filtering for z-grid plot.")

    # Always use 2x2 for 4 pT bins; if fewer, keep layout but hide unused pads
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 7), sharey=True)
    axes_flat = _flatten_axes(axes)

    for i, ptb in enumerate(pt_bins):
        ax = axes_flat[i]
        gd = d[d["pt_bin"] == ptb].copy()

        for lab, g in gd.groupby("curve_label", sort=True):
            g = g.sort_values("z_center")
            x = g["z_center"].astype(float).to_numpy()
            y = g[ycol].astype(float).to_numpy()
            ye = g[yerr].astype(float).to_numpy() if yerr in g.columns else None
            ax.errorbar(x, y, yerr=ye, fmt='o', capsize=2, label=str(lab))

        # title with pt edges if available
        if "pt_lo" in gd.columns and "pt_hi" in gd.columns and len(gd) > 0:
            pt_lo = float(gd["pt_lo"].iloc[0])
            pt_hi = float(gd["pt_hi"].iloc[0])
            ax.set_title(f"$p_T \\in$ [{pt_lo:.2f}, {pt_hi:.2f}] GeV")
        else:
            ax.set_title(f"pt_bin={ptb}")

        ax.set_xlabel("z")
        ax.grid(True, alpha=0.3)

        # Put legend on first pad only
        if i == 0:
            ax.legend(frameon=False)

    # Hide unused pads
    for j in range(len(pt_bins), len(axes_flat)):
        axes_flat[j].axis("off")

    # y label on left column pads
    axes_flat[0].set_ylabel(ylabel)
    if len(axes_flat) > 2:
        axes_flat[2].set_ylabel(ylabel)

    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.0, 1, 0.96])
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True, help="group_id directory name under results/, e.g. grp_LH2_zOv_x0p25_Q23p3_tpq2p0")
    ap.add_argument("--results-dir", default="results", help="base results directory (default: results)")
    ap.add_argument("--min-points", type=int, default=0, help="filter: require n_points >= this value (default: 0)")
    ap.add_argument("--min-prob", type=float, default=0.0, help="filter: require prob >= this value (default: 0.0)")
    args = ap.parse_args()

    group_id = args.group
    results_dir = Path(args.results_dir).resolve()
    fit_csv = results_dir / group_id / "tables" / "fit_parameters.csv"
    if not fit_csv.exists():
        raise SystemExit(f"ERROR: fit_parameters.csv not found: {fit_csv}")

    df = pd.read_csv(fit_csv)

    # Compute z_center if missing but curve_label looks like z0p36 etc.
    if "z_center" not in df.columns:
        def z_from_label(s):
            s = str(s).strip()
            if s.startswith("z"):
                s = s[1:]
            return _safe_float(s.replace("p", "."))
        df["z_center"] = df["curve_label"].map(z_from_label)

    outdir = results_dir / group_id / "PNGs" / "fitparam_plots"
    _ensure_dir(outdir)

    # Core physics plots
    plot_vs_pt(df, "A", "A_err", outdir / "A_vs_pT.png",
               title="Fit parameter A vs $p_T$", ylabel="A",
               min_points=args.min_points, min_prob=args.min_prob)

    plot_vs_pt(df, "B", "B_err", outdir / "B_vs_pT.png",
               title="Fit parameter B vs $p_T$", ylabel="B",
               min_points=args.min_points, min_prob=args.min_prob)

    plot_vs_pt(df, "M0", "M0_err", outdir / "M0_vs_pT.png",
               title="Fit parameter $M_0$ vs $p_T$", ylabel=r"$M_0$",
               min_points=args.min_points, min_prob=args.min_prob)

    # Fit quality
    if "prob" in df.columns:
        plot_prob_vs_pt(df, outdir / "prob_vs_pT.png", min_points=args.min_points)

    # z-dependence in 2x2 grid by pT bin
    plot_AorB_vs_z_grid(df, "A", "A_err", outdir / "A_vs_z_grid.png",
                        title="A vs z (pads: $p_T$ bins)", ylabel="A",
                        min_points=args.min_points, min_prob=args.min_prob)

    plot_AorB_vs_z_grid(df, "B", "B_err", outdir / "B_vs_z_grid.png",
                        title="B vs z (pads: $p_T$ bins)", ylabel="B",
                        min_points=args.min_points, min_prob=args.min_prob)

    print(f"Wrote plots to: {outdir}")
    print("Tip (recommended with many phi bins): try --min-points 10 --min-prob 0.01")


if __name__ == "__main__":
    main()
