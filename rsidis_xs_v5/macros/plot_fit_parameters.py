#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

# Force non-interactive backend (avoids EGL/DRI crashes on farm nodes)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


NUM_COLS = [
    "pt_lo","pt_hi","pt_center","z_lo","z_hi","z_center",
    "n_points","M0","M0_err","A","A_err","B","B_err","chi2","ndf","prob","fit_tier"
]

KEY_COLS = ["group_id","curve_label","setting_id","pt_bin","z_bin"]


def read_fit_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    # Backward-compat: if fit_tier is missing, infer a best-effort tier
    if "fit_tier" not in df.columns:
        # Tier inference: if B is finite -> 2, else if A finite -> 1, else 0
        def infer_tier(row):
            a = row.get("A", np.nan)
            b = row.get("B", np.nan)
            if pd.notna(b) and b != -999:
                return 2
            if pd.notna(a) and a != -999:
                return 1
            return 0
        df["fit_tier"] = df.apply(infer_tier, axis=1)

    # Coerce numerics safely
    for c in NUM_COLS:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Treat -999 as invalid marker
    df.replace(-999, np.nan, inplace=True)

    # Derived columns
    df["chi2ndf"] = np.where(df["ndf"] > 0, df["chi2"] / df["ndf"], np.nan)

    return df


def dedup(df: pd.DataFrame) -> pd.DataFrame:
    """If repeated runs appended duplicate rows, keep the last occurrence."""
    if all(c in df.columns for c in KEY_COLS):
        return df.drop_duplicates(subset=KEY_COLS, keep="last").copy()
    return df


def apply_quality_filters(df: pd.DataFrame, min_prob=None, max_chi2ndf=None) -> pd.DataFrame:
    out = df.copy()
    if min_prob is not None:
        out = out[(out["prob"].notna()) & (out["prob"] >= min_prob)]
    if max_chi2ndf is not None:
        out = out[(out["chi2ndf"].notna()) & (out["chi2ndf"] <= max_chi2ndf)]
    return out


def errorbar_by_curve(ax, subdf: pd.DataFrame, xcol, ycol, yerrcol, labelcol="curve_label"):
    """Group by curve_label and draw default-colored errorbars."""
    for lab, g in subdf.groupby(labelcol):
        x = g[xcol].to_numpy()
        y = g[ycol].to_numpy()
        ye = g[yerrcol].to_numpy() if yerrcol in g.columns else None
        ax.errorbar(x, y, yerr=ye, fmt="o", capsize=2, label=str(lab))


def plot_vs_pt(df: pd.DataFrame, y, yerr, outpath: Path, title: str):
    fig, ax = plt.subplots()
    errorbar_by_curve(ax, df, "pt_center", y, yerr)
    ax.set_xlabel(r"$p_T$ (GeV)")
    ax.set_ylabel(y)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def plot_vs_z_grid(df: pd.DataFrame, y, yerr, outpath: Path, title: str):
    # pads = pt_bin
    pt_bins = sorted(df["pt_bin"].dropna().unique().tolist())
    n = len(pt_bins)
    if n == 0:
        return

    ncols = 2 if n > 1 else 1
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 6), squeeze=False)
    axes = axes.flatten()

    for i, ptb in enumerate(pt_bins):
        ax = axes[i]
        g = df[df["pt_bin"] == ptb]
        errorbar_by_curve(ax, g, "z_center", y, yerr)
        # nice title from pt range if available
        if "pt_lo" in g.columns and "pt_hi" in g.columns and g["pt_lo"].notna().any():
            ptlo = float(g["pt_lo"].dropna().iloc[0])
            pthi = float(g["pt_hi"].dropna().iloc[0])
            ax.set_title(rf"$p_T \in [{ptlo:.2f}, {pthi:.2f}]$ GeV")
        else:
            ax.set_title(f"pt_bin={ptb}")
        ax.set_xlabel("z")
        ax.set_ylabel(y)
        ax.grid(True, alpha=0.3)

    # turn off unused pads
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    # legend only once (top-left)
    axes[0].legend()

    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def plot_prob_vs_pt(df: pd.DataFrame, outpath: Path):
    # log scale is helpful
    fig, ax = plt.subplots()
    # only valid prob
    gdf = df[(df["prob"].notna()) & (df["ndf"] > 0)]
    for lab, g in gdf.groupby("curve_label"):
        ax.plot(g["pt_center"].to_numpy(), g["prob"].to_numpy(), marker="o", label=str(lab))
    ax.set_yscale("log")
    ax.set_xlabel(r"$p_T$ (GeV)")
    ax.set_ylabel("Fit probability")
    ax.set_title(r"Fit quality: probability vs $p_T$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True, help="group id, e.g. grp_LH2_zOv_x0p25_Q23p3_tpq2p0")
    ap.add_argument("--results", default="results", help="top-level results directory (default: results)")
    ap.add_argument("--min-prob", type=float, default=None, help="optional: require prob >= X")
    ap.add_argument("--max-chi2ndf", type=float, default=None, help="optional: require chi2/ndf <= Y")
    ap.add_argument("--dedup", action="store_true", help="drop duplicate rows by (group,curve,pt_bin,z_bin)")
    args = ap.parse_args()

    fitcsv = Path(args.results) / args.group / "tables" / "fit_parameters.csv"
    outdir = Path(args.results) / args.group / "PNGs" / "fitparam_plots"
    outdir.mkdir(parents=True, exist_ok=True)

    df = read_fit_csv(fitcsv)
    if args.dedup:
        df = dedup(df)

    # Print a quick tier summary before any cuts
    print("Rows total:", len(df))
    print("Tier counts:\n", df["fit_tier"].value_counts(dropna=False).sort_index())

    # Optional quality filters (applied AFTER tier selection per plot, below)
    q_min_prob = args.min_prob
    q_max_chi2ndf = args.max_chi2ndf

    # --- M0 plots (tiers 0/1/2) ---
    d0 = df[df["M0"].notna()].copy()
    d0 = apply_quality_filters(d0, q_min_prob, q_max_chi2ndf)
    plot_vs_pt(d0, "M0", "M0_err", outdir / "M0_vs_pT.png", r"Fit parameter $M_0$ vs $p_T$")

    # --- A plots (tiers 1/2 only) ---
    dA = df[(df["fit_tier"] >= 1) & (df["A"].notna())].copy()
    dA = apply_quality_filters(dA, q_min_prob, q_max_chi2ndf)
    plot_vs_pt(dA, "A", "A_err", outdir / "A_vs_pT.png", r"Fit parameter $A$ vs $p_T$")
    plot_vs_z_grid(dA, "A", "A_err", outdir / "A_vs_z_grid.png", r"$A$ vs $z$ (pads: $p_T$ bins)")

    # --- B plots (tier 2 only) ---
    dB = df[(df["fit_tier"] >= 2) & (df["B"].notna())].copy()
    dB = apply_quality_filters(dB, q_min_prob, q_max_chi2ndf)
    plot_vs_pt(dB, "B", "B_err", outdir / "B_vs_pT.png", r"Fit parameter $B$ vs $p_T$")
    plot_vs_z_grid(dB, "B", "B_err", outdir / "B_vs_z_grid.png", r"$B$ vs $z$ (pads: $p_T$ bins)")

    # --- prob vs pT (only where ndf>0 & prob valid) ---
    dp = df[(df["ndf"] > 0) & (df["prob"].notna())].copy()
    dp = apply_quality_filters(dp, q_min_prob, q_max_chi2ndf)
    plot_prob_vs_pt(dp, outdir / "prob_vs_pT.png")

    # Skipped summary
    def nvalid(x): return int(np.sum(pd.notna(x)))
    print("\nValid counts after tier logic (before your optional filters were applied inside plots):")
    print("M0 valid:", nvalid(df["M0"]))
    print("A  valid (tier>=1):", nvalid(df.loc[df["fit_tier"]>=1, "A"]))
    print("B  valid (tier>=2):", nvalid(df.loc[df["fit_tier"]>=2, "B"]))
    print("prob valid (ndf>0):", len(df[(df["ndf"]>0) & (df["prob"].notna())]))

    print("\nWrote plots to:", outdir)


if __name__ == "__main__":
    main()
