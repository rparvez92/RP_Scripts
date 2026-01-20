#!/usr/bin/env python3
"""
plot_xsec_phipq_overlay.py

Reads xsec_phipq_overlay.csv and plots dσ vs φ_pq with:
  - pt bins as subplots (pads)
  - curve_label (z bins) as overlaid curves per pad
  - acceptance-aware fit tiers + robust fitting + optional sigma clipping

Path conventions for rsidis_xs_v5:
  - Script lives in:   rsidis_xs_v5/macros/
  - Run from:          rsidis_xs_v5/
  - Default CSV:       rsidis_xs_v5/results/<group_id>/tables/xsec_phipq_overlay.csv
  - Default output:    rsidis_xs_v5/results/<group_id>/PNGs/<outname>.png

If group_id inside CSV doesn't match directory naming (e.g. zOv vs z0v),
script attempts to resolve a best directory match. If it cannot, it writes
output into the current working directory as fallback (as requested).
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import difflib


# ----------------------------
# Repo/path helpers
# ----------------------------

def repo_root_from_script() -> Path:
    # .../rsidis_xs_v5/macros/plot_xsec_phipq_overlay.py -> root is parents[1]
    return Path(__file__).resolve().parents[1]


def normalize_key(s: str) -> str:
    # normalize for fuzzy comparisons
    return re.sub(r"[^a-z0-9]+", "", s.lower())


def generate_group_variants(gid: str) -> List[str]:
    """
    Generate plausible variants to handle zOv vs z0v and similar issues.
    We keep it conservative (only a few obvious swaps).
    """
    variants = {gid}

    # common confusions: O <-> 0, especially in "...z0v..." segment
    variants.add(gid.replace("Ov", "0v"))
    variants.add(gid.replace("0v", "Ov"))
    variants.add(gid.replace("O", "0"))
    variants.add(gid.replace("0", "O"))

    # also try lower/upper case versions (directory is usually exact, but just in case)
    variants.add(gid.lower())
    variants.add(gid.upper())

    # keep stable ordering
    out = []
    seen = set()
    for v in variants:
        if v not in seen:
            out.append(v)
            seen.add(v)
    return out


def resolve_results_group_dir(results_dir: Path, gid_from_csv: str) -> Optional[Path]:
    """
    Try to map group_id in CSV to an existing directory under results/.
    Returns Path to group directory if found, else None.
    """
    if not results_dir.is_dir():
        return None

    # 1) exact / variant match
    for v in generate_group_variants(gid_from_csv):
        cand = results_dir / v
        if cand.is_dir():
            return cand

    # 2) fuzzy match among existing group dirs
    subdirs = [p for p in results_dir.iterdir() if p.is_dir()]
    names = [p.name for p in subdirs]
    if not names:
        return None

    target = normalize_key(gid_from_csv)
    choices = [normalize_key(n) for n in names]

    # find closest
    close = difflib.get_close_matches(target, choices, n=1, cutoff=0.75)
    if close:
        idx = choices.index(close[0])
        return subdirs[idx]

    return None


def default_csv_path(root: Path, group_id: str) -> Path:
    return root / "results" / group_id / "tables" / "xsec_phipq_overlay.csv"


# ----------------------------
# Plot helpers
# ----------------------------

def parse_z_from_label(label: str) -> float:
    """
    'z0p36' -> 0.36
    """
    m = re.search(r"z(\d+)p(\d+)", label)
    if not m:
        return float("nan")
    a = int(m.group(1))
    b = int(m.group(2))
    return a + b / (10 ** len(m.group(2)))


def nice_pt_title(pt_lo: float, pt_hi: float) -> str:
    return rf"$P_T \in [{pt_lo:.2f},\, {pt_hi:.2f}]$ GeV"


def default_color_cycle(n: int) -> List[str]:
    base = ["k", "r", "b", "g", "m", "c", "y"]
    if n <= len(base):
        return base[:n]
    extra = [f"C{i}" for i in range(n - len(base))]
    return base + extra


def make_grid(n: int) -> Tuple[int, int]:
    ncols = int(math.ceil(math.sqrt(n)))
    nrows = int(math.ceil(n / ncols))
    return nrows, ncols


def safe_rel_err(y: np.ndarray, yerr: np.ndarray) -> np.ndarray:
    denom = np.where(np.abs(y) > 0, np.abs(y), np.nan)
    return yerr / denom


# ----------------------------
# Fit models & fitting
# ----------------------------

def model_tier(phi: np.ndarray, p: np.ndarray, tier: int) -> np.ndarray:
    # tier 0: a0
    # tier 1: a0 + a1 cos(phi)
    # tier 2: a0 + a1 cos(phi) + a2 cos(2phi)
    a0 = p[0]
    if tier == 0:
        return a0 + 0.0 * phi
    a1 = p[1]
    if tier == 1:
        return a0 + a1 * np.cos(phi)
    a2 = p[2]
    return a0 + a1 * np.cos(phi) + a2 * np.cos(2.0 * phi)


def initial_guess(phi: np.ndarray, y: np.ndarray, tier: int) -> np.ndarray:
    a0 = float(np.nanmedian(y))
    if tier == 0:
        return np.array([a0], dtype=float)
    c1 = np.cos(phi)
    a1 = float(np.nan_to_num(np.dot(y - a0, c1) / (np.dot(c1, c1) + 1e-12)))
    if tier == 1:
        return np.array([a0, a1], dtype=float)
    c2 = np.cos(2.0 * phi)
    a2 = float(np.nan_to_num(np.dot(y - a0, c2) / (np.dot(c2, c2) + 1e-12)))
    return np.array([a0, a1, a2], dtype=float)


@dataclass
class FitResult:
    tier: int
    params: np.ndarray
    chi2: float
    ndf: int
    mask_used: np.ndarray


def fit_curve(
    phi: np.ndarray,
    y: np.ndarray,
    yerr: np.ndarray,
    tier: int,
    robust: bool,
    sigma_clip: bool,
    pull_cut: float,
    max_clip_iter: int,
) -> Optional[FitResult]:
    base = np.isfinite(phi) & np.isfinite(y) & np.isfinite(yerr) & (yerr > 0)
    if base.sum() < (tier + 2):
        return None

    mask = base.copy()
    loss = "huber" if robust else "linear"
    f_scale = 1.0

    def resid(p: np.ndarray, m: np.ndarray) -> np.ndarray:
        return (y[m] - model_tier(phi[m], p, tier)) / yerr[m]

    p0 = initial_guess(phi[mask], y[mask], tier)

    for _ in range(max_clip_iter if sigma_clip else 1):
        if mask.sum() < (tier + 2):
            return None

        res = least_squares(
            lambda p: resid(p, mask),
            x0=p0,
            loss=loss,
            f_scale=f_scale,
            max_nfev=5000,
        )
        pfit = res.x
        pulls = resid(pfit, mask)

        if sigma_clip:
            bad = np.abs(pulls) > pull_cut
            if not np.any(bad):
                p0 = pfit
                break
            idx_mask = np.where(mask)[0]
            mask[idx_mask[bad]] = False
            p0 = pfit
        else:
            p0 = pfit
            break

    final_pulls = resid(p0, mask)
    chi2 = float(np.sum(final_pulls**2))
    ndf = int(mask.sum() - len(p0))
    return FitResult(tier=tier, params=p0, chi2=chi2, ndf=ndf, mask_used=mask)


# ----------------------------
# Acceptance-aware selection
# ----------------------------

@dataclass
class TierCuts:
    frac_t0: float
    frac_t1: float
    frac_t2: float
    max_rel_stat: float


def eligible_mask_for_tier(
    sub: pd.DataFrame,
    tier: int,
    prefer_valid: bool,
    cuts: TierCuts,
    allow_negative_xsec: bool,
) -> np.ndarray:
    x = sub["xsec"].to_numpy(float)
    xe = sub["xsec_stat"].to_numpy(float)

    if "yield_sim" in sub.columns:
        ysim = sub["yield_sim"].to_numpy(float)
    else:
        ysim = np.full_like(x, np.nan)

    valid = sub["valid"].to_numpy(int) if "valid" in sub.columns else np.ones_like(x, dtype=int)

    m = np.isfinite(x) & np.isfinite(xe) & (xe > 0)

    if not allow_negative_xsec:
        m &= (x > 0)

    rel = safe_rel_err(x, xe)
    m &= np.isfinite(rel) & (rel <= cuts.max_rel_stat)

    # acceptance proxy by yield_sim fraction (per curve & pt pad)
    if np.any(np.isfinite(ysim)):
        ymax = np.nanmax(ysim)
        if np.isfinite(ymax) and ymax > 0:
            frac = ysim / ymax
            thr = {0: cuts.frac_t0, 1: cuts.frac_t1, 2: cuts.frac_t2}.get(tier, cuts.frac_t2)
            m &= np.isfinite(frac) & (frac >= thr)

    if prefer_valid:
        m &= (valid == 1)

    return m


def auto_choose_fit_tier(npts_t2: int, npts_t1: int, npts_t0: int) -> int:
    """
    Conservative auto-tier selection:
      - use tier2 only if plenty of good points
      - else tier1 if enough
      - else tier0
    """
    if npts_t2 >= 8:
        return 2
    if npts_t1 >= 5:
        return 1
    if npts_t0 >= 3:
        return 0
    return 0


# ----------------------------
# Main
# ----------------------------

def main() -> int:
    root = repo_root_from_script()
    results_dir = root / "results"

    ap = argparse.ArgumentParser()

    # Optional positional CSV: if omitted, we build it from --group-id
    ap.add_argument(
        "csv",
        nargs="?",
        default=None,
        help="CSV path (relative to rsidis_xs_v5) OR absolute path. "
             "If omitted, uses results/<group_id>/tables/xsec_phipq_overlay.csv",
    )
    ap.add_argument("--group-id", default=None,
                    help="Group ID (directory under results/). If csv is omitted, this is required.")

    ap.add_argument("-o", "--out", default="xsec_vs_phipq_overlay.png",
                    help="Output PNG filename (name only or path). "
                         "If name only, goes to results/<group>/PNGs/ unless group mismatch fallback.")
    ap.add_argument("--prefer-valid", action="store_true",
                    help="Require valid==1 for fit eligibility (acceptance-aware)")
    ap.add_argument("--fit-tier", default="auto", choices=["auto", "0", "1", "2"],
                    help="Fit model tier: auto/0/1/2")

    ap.add_argument("--frac-t0", type=float, default=0.20)
    ap.add_argument("--frac-t1", type=float, default=0.10)
    ap.add_argument("--frac-t2", type=float, default=0.05)
    ap.add_argument("--max-rel-stat", type=float, default=0.80)

    ap.add_argument("--allow-negative-xsec", action="store_true")
    ap.add_argument("--robust", action="store_true")
    ap.add_argument("--sigma-clip", action="store_true")
    ap.add_argument("--pull-cut", type=float, default=3.5)
    ap.add_argument("--max-clip-iter", type=int, default=6)

    ap.add_argument("--inflate-low-acc", action="store_true",
                    help="Inflate errors for low-yield_sim points in fit to reduce leverage (conservative).")

    ap.add_argument("--title", default=None)
    ap.add_argument("--ylog", action="store_true")
    ap.add_argument("--xlim", nargs=2, type=float, default=None)

    # optional view helper (off by default)
    ap.add_argument("--smart-ylim", action="store_true",
                    help="Set pad ylim using 99th percentile of (y+err) of eligible points to avoid huge bars dominating view.")

    args = ap.parse_args()

    # Determine CSV path (relative to root if needed)
    if args.csv is None:
        if args.group_id is None:
            print("ERROR: Provide either CSV path or --group-id.", file=sys.stderr)
            return 2
        csv_path = default_csv_path(root, args.group_id)
    else:
        p = Path(args.csv)
        csv_path = p if p.is_absolute() else (Path.cwd() / p)

    if not csv_path.exists():
        print(f"ERROR: CSV not found: {csv_path}", file=sys.stderr)
        return 2

    df = pd.read_csv(csv_path)

    # Determine group_id (prefer from CSV if present)
    gid_csv = None
    if "group_id" in df.columns:
        # could be repeated; take the first non-null
        s = df["group_id"].dropna()
        if len(s) > 0:
            gid_csv = str(s.iloc[0])

    group_id_effective = args.group_id or gid_csv

    # Resolve output directory
    out_arg = Path(args.out)
    out_is_path = (out_arg.parent != Path("."))  # user provided directory components
    out_name_only = out_arg.name

    resolved_group_dir = None
    if group_id_effective is not None:
        resolved_group_dir = resolve_results_group_dir(results_dir, group_id_effective)

    if out_is_path:
        # user explicitly gave a path -> respect it (relative to CWD)
        out_path = out_arg if out_arg.is_absolute() else (Path.cwd() / out_arg)
        out_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        # name only -> place into results/<group>/PNGs if we can resolve, else fallback to CWD
        if resolved_group_dir is not None:
            png_dir = resolved_group_dir / "PNGs"
            png_dir.mkdir(parents=True, exist_ok=True)
            out_path = png_dir / out_name_only
        else:
            print(
                f"WARNING: Could not resolve results/<group_id> directory for group_id='{group_id_effective}'.\n"
                f"         This can happen if CSV uses zOv but directory uses z0v.\n"
                f"         Falling back to writing output in current directory: {Path.cwd()}",
                file=sys.stderr,
            )
            out_path = Path.cwd() / out_name_only

    # Basic sanity checks for required columns
    required = ["pt_bin", "pt_lo", "pt_hi", "phi_center", "curve_label", "xsec", "xsec_stat"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"ERROR: Missing required columns in CSV: {missing}", file=sys.stderr)
        return 2

    # pt bin ordering
    pt_bins = (
        df[["pt_bin", "pt_lo", "pt_hi"]]
        .drop_duplicates()
        .sort_values(["pt_lo", "pt_hi"])
        .reset_index(drop=True)
    )
    pt_bin_list = pt_bins["pt_bin"].tolist()

    # curve ordering by z-value when possible
    z_labels = sorted(df["curve_label"].unique().tolist(),
                      key=lambda s: (np.nan_to_num(parse_z_from_label(s), nan=1e9), s))
    colors = default_color_cycle(len(z_labels))
    color_map = {zl: colors[i] for i, zl in enumerate(z_labels)}

    nrows, ncols = make_grid(len(pt_bin_list))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.4 * ncols, 4.8 * nrows), squeeze=False)

    cuts = TierCuts(frac_t0=args.frac_t0, frac_t1=args.frac_t1, frac_t2=args.frac_t2, max_rel_stat=args.max_rel_stat)

    for i, ptb in enumerate(pt_bin_list):
        r, c = divmod(i, ncols)
        ax = axes[r][c]

        sub_pt = df[df["pt_bin"] == ptb].copy()
        if sub_pt.empty:
            ax.set_visible(False)
            continue

        pt_lo = float(sub_pt["pt_lo"].iloc[0])
        pt_hi = float(sub_pt["pt_hi"].iloc[0])
        ax.set_title(nice_pt_title(pt_lo, pt_hi))

        # for smart ylim (optional)
        eligible_y_plus = []

        for zl in z_labels:
            sub = sub_pt[sub_pt["curve_label"] == zl].copy()
            if sub.empty:
                continue

            sub = sub.sort_values("phi_center")
            phi = sub["phi_center"].to_numpy(float)
            y = sub["xsec"].to_numpy(float)
            ye = sub["xsec_stat"].to_numpy(float)

            valid_flag = sub["valid"].to_numpy(int) if "valid" in sub.columns else np.ones_like(y, dtype=int)
            m_invalid = (valid_flag == 0)

            # plot invalid/faint points for visibility
            if np.any(m_invalid):
                ax.errorbar(phi[m_invalid], y[m_invalid], yerr=ye[m_invalid],
                            fmt="o", ms=4, lw=1, color=color_map[zl], alpha=0.20, capsize=2)

            # eligibility masks for tiers
            m0 = eligible_mask_for_tier(sub, tier=0, prefer_valid=args.prefer_valid, cuts=cuts,
                                        allow_negative_xsec=args.allow_negative_xsec)
            m1 = eligible_mask_for_tier(sub, tier=1, prefer_valid=args.prefer_valid, cuts=cuts,
                                        allow_negative_xsec=args.allow_negative_xsec)
            m2 = eligible_mask_for_tier(sub, tier=2, prefer_valid=args.prefer_valid, cuts=cuts,
                                        allow_negative_xsec=args.allow_negative_xsec)

            if args.fit_tier == "auto":
                tier = auto_choose_fit_tier(int(m2.sum()), int(m1.sum()), int(m0.sum()))
            else:
                tier = int(args.fit_tier)

            mfit = {0: m0, 1: m1, 2: m2}[tier]

            # always plot eligible points in solid style
            if np.any(mfit):
                ax.errorbar(phi[mfit], y[mfit], yerr=ye[mfit],
                            fmt="o", ms=5, lw=1.2, color=color_map[zl], alpha=0.95, capsize=2, label=zl)
                if args.smart_ylim:
                    eligible_y_plus.extend(list((y[mfit] + ye[mfit])[np.isfinite(y[mfit] + ye[mfit])]))

            # Effective errors for fit (optional conservative inflation by acceptance)
            if args.inflate_low_acc and "yield_sim" in sub.columns:
                ysim = sub["yield_sim"].to_numpy(float)
                ymax = np.nanmax(ysim) if np.any(np.isfinite(ysim)) else np.nan
                if np.isfinite(ymax) and ymax > 0:
                    frac = np.clip(ysim / ymax, 1e-3, 1.0)
                    ye_eff = ye / np.sqrt(frac)
                else:
                    ye_eff = ye.copy()
            else:
                ye_eff = ye.copy()

            # Fit on eligible points only
            fitres = None
            if np.any(mfit) and mfit.sum() >= (tier + 2):
                fitres = fit_curve(
                    phi=phi[mfit],
                    y=y[mfit],
                    yerr=ye_eff[mfit],
                    tier=tier,
                    robust=args.robust,
                    sigma_clip=args.sigma_clip,
                    pull_cut=args.pull_cut,
                    max_clip_iter=args.max_clip_iter,
                )

            if fitres is not None:
                xx = np.linspace(0.0, 2.0 * math.pi, 400)
                yy = model_tier(xx, fitres.params, fitres.tier)
                ax.plot(xx, yy, "-", color=color_map[zl], lw=2.0, alpha=0.95)

        ax.set_xlabel(r"$\phi_{pq}$ (rad)")
        ax.set_ylabel(r"$d\sigma$")

        if args.ylog:
            ax.set_yscale("log")

        if args.xlim is not None:
            ax.set_xlim(args.xlim[0], args.xlim[1])
        else:
            ax.set_xlim(0.0, 2.0 * math.pi)

        ax.grid(True, which="both", alpha=0.35)
        ax.legend(frameon=False, fontsize=11)

        # optional: avoid giant error bars dominating view
        if args.smart_ylim and len(eligible_y_plus) > 0 and not args.ylog:
            y99 = float(np.nanpercentile(np.array(eligible_y_plus), 99))
            if np.isfinite(y99) and y99 > 0:
                ax.set_ylim(bottom=0.0, top=1.25 * y99)

    # hide unused pads
    for j in range(len(pt_bin_list), nrows * ncols):
        rr, cc = divmod(j, ncols)
        axes[rr][cc].set_visible(False)

    if args.title:
        fig.suptitle(args.title, fontsize=16)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    print(f"Wrote: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
