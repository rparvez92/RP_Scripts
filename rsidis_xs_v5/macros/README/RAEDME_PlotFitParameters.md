# README — `PlotFitParameters.C` (R‑SIDIS `rsidis_xs_v5`)
**Macro:** `macros/PlotFitParameters.C`  
**README filename:** `README_PlotFitParameters.md`  
**Input:** `fit_parameters_<mode>_tierOn.csv` (or `_tierOff.csv`) produced by `PlotCoinXsecFromTable_wTier.C`  
**Output:** a directory of PNG summary plots under `results/<id>/PNGs/FitParamPlots_*`

---

## 0) What this macro does (in plain language)
`PlotFitParameters.C` is the **summary-visualization step** for the Fourier-fit results. It reads the fit-parameter table (CSV) written by `PlotCoinXsecFromTable_wTier.C` and generates:

- **\(M_0\), \(A\), \(B\) vs \(z\)** for each \(p_T\) bin (curves = \(p_T\) bins)
- **\(M_0\), \(A\), \(B\) vs \(p_T\)** for each \(z\) bin (curves = \(z\) bins)
- **\(\chi^2/\mathrm{ndf}\) vs \(z\)** and **\(\chi^2/\mathrm{ndf}\) vs \(p_T\)** (both in log-y scale)

This macro is intentionally “read-only physics”: it does **not** recompute yields, cross sections, or fits. It only **filters**, **groups**, and **plots**.

---

## 1) Where these fit parameters come from (context)
Upstream, `PlotCoinXsecFromTable_wTier.C` fits the azimuthal dependence of the extracted coincidence cross section in each \((p_T, z)\) bin using a tiered model:

- Tier 0: constant \(M_0\)
- Tier 1: \(M_0(1 + A\cos\phi)\)
- Tier 2: \(M_0(1 + A\cos\phi + B\cos 2\phi)\)

This produces a per-bin summary:
- \(M_0\) (scale)
- \(A\) (cos φ modulation amplitude)
- \(B\) (cos 2φ modulation amplitude)
plus diagnostics:
- `chi2`, `ndf`, `chi2_ndf`, `prob`
- `fit_tier`, `fit_status`, `cov_status`, `fit_flag_bits`
- `valid_fit_default` and component validity flags

**This macro plots those summaries across bins** to let you quickly see trends and stability.

---

## 2) How to run

From `rsidis_xs_v5/`:

```bash
# Group-mode fit table:
root -l -b -q 'macros/PlotFitParameters.C("results/<group_id>/tables/fit_parameters_group_tierOn.csv")'

# Single-setting fit table:
root -l -b -q 'macros/PlotFitParameters.C("results/<setting_id>/tables/fit_parameters_single_tierOn.csv")'

# For tierOff tables:
root -l -b -q 'macros/PlotFitParameters.C("results/<group_id>/tables/fit_parameters_group_tierOff.csv")'
```

If you call it with an empty path, it prints the usage string and returns.

---

## 3) Outputs and where they go

### 3.1 Output directory naming (tierOn vs tierOff)
The output directory is inferred from the input CSV path by `InferOutDir()`:

- if CSV path contains `tieron` → `.../PNGs/FitParamPlots_tierOn/`
- if CSV path contains `tieroff` → `.../PNGs/FitParamPlots_tierOff/`
- otherwise → `.../PNGs/FitParamPlots/` (fallback)

This is designed so you can run the macro on both tierOn and tierOff tables without overwriting.

### 3.2 Output files written
Inside the output directory, the macro writes:

- `M0_vs_z.png`
- `A_vs_z.png`
- `B_vs_z.png`
- `M0_vs_pt.png`
- `A_vs_pt.png`
- `B_vs_pt.png`
- `chi2ndf_vs_z.png`
- `chi2ndf_vs_pt.png`

At the end it prints:
```
Wrote FitParamPlots to: <outDir>
```

---

## 4) Physics meaning of the plotted quantities

### 4.1 The fit model being summarized
For each \((p_T, z)\) bin, the upstream fit approximates:
\[
\sigma(\phi_{pq}) \approx M_0\left[1 + A\cos(\phi_{pq}) + B\cos(2\phi_{pq})\right]
\]
(when tier 2 is used; tiers 0/1 drop terms).

**Interpretation:**
- **\(M_0\)**: overall scale (roughly the \(\phi\)-averaged cross section level in that bin)
- **\(A\)**: normalized \(\cos\phi\) modulation amplitude
- **\(B\)**: normalized \(\cos 2\phi\) modulation amplitude

These are not “fundamental constants”; they are compact summaries of the measured azimuthal dependence within the analysis definition (cuts, binning, random subtraction, positron/dummy subtraction, model dependence, etc.).

### 4.2 Why plot vs \(z\) and vs \(p_T\)?
- **Trends vs \(p_T\)** often reveal transverse-momentum physics and acceptance/edge issues.
- **Trends vs \(z\)** often reveal fragmentation/kinematic dependence and can separate regions where the modulation behaves differently.
- Seeing the same parameter plotted in both projections helps distinguish:
  - a real smooth physics trend
  - from an artifact that only appears in one projection due to missing bins.

### 4.3 Why include \(\chi^2/\mathrm{ndf}\) in log scale?
\(\chi^2/\mathrm{ndf}\) is a fit-quality diagnostic:
- values near 1 suggest reasonable agreement given the error model
- very large values suggest mismodeling, underestimated errors, or unstable points

Log scale is used because \(\chi^2/\mathrm{ndf}\) can span orders of magnitude when fits go bad.

**Important:** a “good” \(\chi^2/\mathrm{ndf}\) is necessary but not sufficient; covariance issues, parameter-at-limit flags, and coverage tests are encoded upstream in `fit_flag_bits` and `valid_fit_default`.

---

## 5) Filters/guards and their physics consequences

### 5.1 The central filter: `valid_fit_default`
The macro, by default, only plots rows satisfying:
- `valid_fit_default == 1`

This behavior is controlled by:

```cpp
PlotConfig cfg;
cfg.require_valid_fit_default = true;
```

and implemented in:
- `PassRowPhysics(const Row&, const PlotConfig&)`

**Why this filter exists:**  
It ensures that these summary plots represent **physics-usable** fits, not numerical failures or untrustworthy tier-forced results.

**Could this lose good physics?**  
Yes, in principle:
- A real modulation could exist in a bin where the fit was flagged invalid due to poor \(\phi\) coverage, large uncertainties, or dominance by one point.
- In that case, the correct interpretation is *not* “physics absent”; it is “physics not constrained with current binning/statistics/cuts”.

**Diagnostic option:**  
For debugging, you can set:
```cpp
cfg.require_valid_fit_default = false;
```
This will plot everything the CSV contains, including invalid fits. That’s useful to understand *why* bins are failing, but should not be treated as final physics.

### 5.2 Expected binning and “missing bin” notes
This macro contains hard-coded “expected” bins used only to generate a small reminder note in the plot:

```cpp
static const std::vector<double> kExpectedZCenters = {0.36, 0.50, 0.67, 0.90};
static const std::vector<std::pair<double,double>> kExpectedPtBins = {
  {0.0, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.4}
};
```

If an expected z or pT bin has no plotted points, a note appears in the top margin, e.g.:
- `No data for z={0.36,0.67}`

**This is not a physics cut.**  
It is purely a QA hint:
- maybe a bin truly has no usable fit,
- or the fit table schema/bins changed and the “expected” arrays should be updated.

### 5.3 Symmetric axis for A and B
For `A` and `B`, the y-axis range is forced to be symmetric around 0 using:
- `ComputeYRangeSymmetric0(...)`

This avoids misleading visuals where all values are positive and the plot auto-range hides the negative side (or vice versa).

The minimum half-range is controlled by:
```cpp
cfg.symmetric_AB_min_abs = 0.25;
```

---

## 6) Code structure (schematic map)

### 6.1 High-level call flow
```
PlotFitParameters(csvPath)
 ├─ ReadRows(csvPath)                 → vector<Row>
 ├─ outDir = InferOutDir(csvPath)     → "results/.../PNGs/FitParamPlots_*"
 ├─ exPt = ExemplarByPtBin(rows)      → legend label helpers
 ├─ exZ  = ExemplarByZBin(rows)
 ├─ BuildGraphsVsZ(...) + Plot1D...   → M0/A/B vs z, chi2ndf vs z
 ├─ BuildGraphsVsPt(...) + Plot1D...  → M0/A/B vs pt, chi2ndf vs pt
 └─ print outDir
```

### 6.2 “Graph building” concept
A `GraphSet` is:
```cpp
struct GraphSet { std::map<int, TGraphErrors*> g; };
```

- For **vs z** plots: each graph corresponds to one `pt_bin`
- For **vs pt** plots: each graph corresponds to one `z_bin`

---

## 7) Function-by-function reference

### 7.1 Constants and simple helpers
- **`kMissing`**  
  Sentinel used in the CSVs (`-999`) to represent “undefined/uncomputed”. This macro skips points where a needed value equals `kMissing`.

- **`Trim(const std::string&)`**  
  Strips leading/trailing whitespace (used for parsing robustness).

- **`SplitCSV(const std::string&)`**  
  Minimal CSV splitter with quoted-field support.

- **`IsFinite(double)`**  
  Wrapper for `std::isfinite` checks.

- **`NearlyEqual(a,b,eps)`**  
  Used for matching expected bin values when writing “missing bin” notes.

- **`InferBinIndex(center, centers, tol)`**  
  For a given `center`, finds the nearest expected center within tolerance. Used when bin indices are missing or inconsistent.

- **`PtLabel(lo,hi)` / `ZLabel(z)`**  
  Formatting helpers for note strings.

### 7.2 Data container
- **`struct Row`**  
  A parsed row from the `fit_parameters_*.csv`. Includes:
  - kinematics: `pt_bin`, `z_bin`, `pt_lo/hi/center`, `z_lo/hi/center`
  - fit parameters: `M0,A,B` and their errors
  - diagnostics: `chi2_ndf`, `prob`, `fit_tier`, statuses, flags
  - validity flags: `valid_fit_default`, `valid_M0_default`, `valid_A_default`, `valid_B_default`
  - text: `note`

### 7.3 Row key helpers
- **`ZKeyFromRow(const Row&)`**  
  Uses `z_bin` if present; otherwise infers from `z_center` using `kExpectedZCenters`.

- **`PtKeyFromRow(const Row&)`**  
  Uses `pt_bin` if present; otherwise infers from `pt_center` using `kExpectedPtBins`.

- **`ZLabelFromIndex(int)` / `PtBinLabelFromIndex(int)`**  
  Fallback legend labels if exemplar rows are missing.

### 7.4 Robust parsing utilities
- **`SafeAtoi(string, int&)`**  
  Converts string to int without throwing; leaves defaults if conversion fails.

- **`SafeAtof(string, double&)`**  
  Converts string to double safely.

- **`ReadRows(csvPath)`**  
  Reads the CSV header, builds a name→index map, and fills `std::vector<Row>`.  
  Key design: column lookup by name makes it resilient to column order changes.

### 7.5 Output directory logic
- **`InferOutDir(csvPath)`**  
  Converts `.../tables/fit_parameters_*.csv` to `.../PNGs/FitParamPlots_*`.  
  Detects tierOn/off tags from the filename.

### 7.6 Plot selection (“physics rows”)
- **`struct PlotConfig`**  
  Plot styling + selection behavior:
  - `require_valid_fit_default` (main physics guard)
  - margins, font sizes, marker sizes
  - symmetric range handling for A/B

- **`PassRowPhysics(const Row&, cfg)`**  
  Returns true only if the row is accepted by the current plot policy.

### 7.7 Styling and annotations
- **`ApplyCanvasStyle()`**  
  Turns off ROOT stat box, removes title, enables axis ticks on both sides.

- **`DrawHeader(cfg, title)`**  
  Places a left-aligned title in the top margin using `TLatex`.

- **`DrawMissingNote(cfg, msg)`**  
  Places “missing bins” note under the title, still inside the top margin.

- **`MakeTopLegend(cfg, ...)`**  
  Places a legend in the top margin on the right (transparent background).

### 7.8 Missing-bin notes
- **`MissingZNote(presentZ)`**  
  Compares present z centers in the plotted data to `kExpectedZCenters` and returns a string listing absent z values.

- **`MissingPtNote(presentPtBins)`**  
  Same idea for pt bins vs `kExpectedPtBins`.

These notes are informational only.

### 7.9 Graph management and construction
- **`struct GraphSet`**  
  Holds a map from key → `TGraphErrors*`.

- **`DeleteGraphs(GraphSet&)`**  
  Deletes owned graphs (prevents memory leaks between plots).

- **`StyleGraph(TGraphErrors*, idx, markerSize)`**  
  Applies a color/marker palette and sets marker size/line width.

- **`CollectAllY(GraphSet&, y, ey)`**  
  Collects all y-values (and errors) across curves for auto range computation.

- **`ComputeYRangeAuto(y, ey, padFrac)`**  
  Computes y-range from min(y−ey) to max(y+ey) with padding.

- **`ComputeYRangeSymmetric0(y, ey, minAbs)`**  
  Computes a symmetric y-range around zero (used for A/B).

- **`BuildGraphsVsZ(rows, cfg, yField, eyField, ...)`**  
  Builds per-pt-bin curves of y vs z:
  - x = `z_center`
  - key = `pt_bin`
  - y chosen by `yField` string (M0, A, B, chi2_ndf, etc.)
  - error chosen by `eyField` string (M0_err, A_err, B_err)
  Sorts points in x per curve.

- **`BuildGraphsVsPt(rows, cfg, yField, eyField, ...)`**  
  Builds per-z-bin curves of y vs pt:
  - x = `pt_center`
  - key = `ZKeyFromRow(r)` (robust even if z_bin missing)
  Sorts points in x per curve.

### 7.10 Plotting routine
- **`Plot1D_MultiCurves(outPng, cfg, gs, title, xTitle, yTitle, logy, missingNote, legendIsPtBin, exemplarRowByKey)`**  
  The single “plot engine” used by all outputs:
  - creates canvas with custom margins (large top margin)
  - chooses y-range (auto or symmetric 0 for A/B)
  - for logy: enforces positive min range using smallest positive y
  - draws frame axes
  - draws title + missing note in top margin
  - draws all curves with legend (labels from exemplar rows when available)
  - saves to PNG

### 7.11 Entry point
- **`PlotFitParameters(const char* csvPathC)`**  
  Orchestrates everything:
  - reads CSV
  - builds exemplar rows (for clean legend labels)
  - makes the 8 canonical plots
  - writes them to the inferred output directory

---

## 8) Common failure modes and how to interpret them

### 8.1 “No points for plot …”
This means that after `PassRowPhysics` filtering (usually `valid_fit_default==1`), no rows were available for that plot type. Causes:
- too strict upstream validity (many bins failing)
- wrong CSV path
- schema mismatch (columns missing so values stay `kMissing`)

Debug steps:
- temporarily set `cfg.require_valid_fit_default = false` to see whether rows exist but are invalid
- check your input CSV actually has non-`-999` values for that field

### 8.2 Missing-bin notes showing many bins
This usually means:
- you changed binning upstream (so update `kExpectedZCenters` / `kExpectedPtBins`)
- or many bins are failing validity cuts

Do not interpret missing-note text as a physics conclusion; it is a QA hint.

### 8.3 \(\chi^2/\mathrm{ndf}\) plot range looks weird
Because the plot is log-y, the range is set using the smallest positive y. If your curves include zeros or negative values (should not happen for chi2/ndf), they are ignored for minPos selection.

If you see spikes:
- they reflect bins where the fit is still marked valid but has relatively larger chi2/ndf; check upstream `fit_flag_bits` and `prob`.

---

## 9) FAQ

### Q: Why only `valid_fit_default==1` rows by default?
Because this macro is meant to produce “physics-ready” summary plots.  
Invalid fits are still valuable (debugging), but they can dominate the visual narrative and mislead casual readers.

### Q: Could invalid fits still contain real physics?
Yes—but “invalid” usually means the fit is not constrained (coverage, errors, domination), not that the physics is absent. Treat those bins as “needs better statistics/better binning”, not as “zero modulation”.

### Q: Why are A and B plotted with symmetric axes?
Because the sign matters physically/interpretively. A non-symmetric auto-range can visually exaggerate small positive-only fluctuations or hide small negative deviations.

### Q: Why doesn’t this macro plot `dominant_wfrac`?
The fit CSV contains it, but `Row` does not parse that column in this macro (commented in the file). Add it to `Row` + `ReadRows` + `BuildGraphs*` if you want a dominant-weight diagnostic plot.

---

## 10) Practical customization notes

### 10.1 If your binning changes
Update the “expected bin” arrays near the top:
- `kExpectedZCenters`
- `kExpectedPtBins`

Otherwise the missing-bin note becomes misleading.

### 10.2 If you want to include invalid fits for diagnostics
Set:
```cpp
cfg.require_valid_fit_default = false;
```

Optionally, you can also create two runs:
- “physics” plots (valid only)
- “debug” plots (all rows)

### 10.3 Add new summary plots (easy)
To add a new plot:
1. Ensure the quantity is parsed into `Row` in `ReadRows`.
2. Build graphs with `BuildGraphsVsZ` or `BuildGraphsVsPt`.
3. Call `Plot1D_MultiCurves`.

Examples that are already supported in `BuildGraphs*`:
- `prob`
- `fit_tier`
- `npts_fit`
- `A_sig`, `B_sig` (vs z only currently)

---

## 11) Final takeaway
- `PlotFitParameters.C` is the quick “dashboard” for your modulation parameters.
- It is conservative by default (plots only physics-usable fits).
- It provides two complementary projections (vs z and vs pt) plus fit-quality diagnostics.
- It does not decide physics; it helps you **see** whether the physics summary is stable, smooth, and trustworthy.

---
