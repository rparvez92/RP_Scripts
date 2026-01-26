# README — `PlotCoinXsecFromTable_wTier.C` (R‑SIDIS `rsidis_xs_v5`)
**Macro:** `macros/PlotCoinXsecFromTable_wTier.C`  
**README filename:** `README_PlotCoinXsecFromTable_wTier.md`  
**Consumes:** CSV table(s) produced by `TableCoinXsec.C`  
**Produces:** multi‑pad PNG plots + tiered Fourier fit tables + a “parsed points” audit table

---

## 0) What this macro does (mental model)
This macro takes the **bin-by-bin coincidence cross sections** already computed by `TableCoinXsec.C` and answers two practical questions:

1. **Visualization:**  
   “What does \(d\sigma\) vs \(\phi_{pq}\) look like in each \(p_T\) bin, and how do curves (by \(z\)) compare?”

2. **Quantification (physics fit):**  
   “Can I summarize the azimuthal modulation in each \(p_T\) bin using a stable Fourier‑like parameterization?”
   \[
   \sigma(\phi_{pq}) \approx M_0\left[1 + A\cos(\phi_{pq}) + B\cos(2\phi_{pq})\right]
   \]
   where \(M_0\) is an overall scale, \(A\) is the \(\cos\phi\) amplitude, and \(B\) is the \(\cos2\phi\) amplitude.

Crucially, it implements a **tiered fit strategy** (“tiering”) so the fit does **not** become numerically ill-conditioned when \(\phi\) coverage is incomplete or statistics are poor.

---

## 1) Inputs and outputs

### 1.1 Inputs
This macro reads one of the following, depending on run mode:

- **Single setting mode** (manifest path passed):
  - `results/<setting_id>/tables/xsec_phipq_z_pt_overlayed_single.csv`

- **Group overlay mode** (group `.list` passed):
  - `results/<group_id>/tables/xsec_phipq_z_pt_overlayed_group.csv`

These CSVs are written by `TableCoinXsec.C`. They already contain:
- `xsec`, `xsec_err`
- `phi_center`, `pt_bin`, `z_bin` (single mode) or `curve_label` (group mode)
- `valid_default` and `flag_bits` quality metadata

### 1.2 Outputs
Given `tierOn` or `tierOff`, a `tierTag` is created:

- `tierOn` → `tierTag = "tierOn"`
- `tierOff` → `tierTag = "tierOff"`

The macro writes:

1. **PNG plot**
   - Single:
     - `results/<setting_id>/PNGs/xsec_phipq_z_pt_overlayed_single_<tierTag>.png`
   - Group:
     - `results/<group_id>/PNGs/xsec_phipq_z_pt_overlayed_group_<tierTag>.png`

2. **Fit parameter table**
   - Single:
     - `results/<setting_id>/tables/fit_parameters_single_<tierTag>.csv`
   - Group:
     - `results/<group_id>/tables/fit_parameters_group_<tierTag>.csv`

3. **Parsed points audit table** (highly recommended for debugging)
   - Single:
     - `results/<setting_id>/tables/parsed_dataPoints_single_<tierTag>.csv`
   - Group:
     - `results/<group_id>/tables/parsed_dataPoints_group_<tierTag>.csv`

The “parsed points” CSV contains **every row from the input table**, plus:
- whether it was plotted
- why it was dropped from the plot (if dropped)
- whether it was used in the fit
- why it was excluded from the fit (if excluded)
- fit-predicted value at that \(\phi\), residual, and pull (if a fit exists)

This file is your best friend when someone asks:  
“Why is this point missing?” or “Why didn’t the fit include that bin?”

---

## 2) How to run

From `rsidis_xs_v5/`:

```bash
# Single mode (manifest)
root -l -b -q 'macros/PlotCoinXsecFromTable_wTier.C("settings/.../manifest.txt","results","settings")'

# Group mode (.list of curves/settings)
root -l -b -q 'macros/PlotCoinXsecFromTable_wTier.C("groups/.../grp_*.list","results","settings")'

# Explicit example with options:
#   tierOn=true, savePNGs=true, saveFitCsv=true
root -l -b -q 'macros/PlotCoinXsecFromTable_wTier.C("groups/.../grp_pass4_....list","results","settings",true,true,true)'

# Compare against "no tiering" reference (forces tier2 everywhere)
root -l -b -q 'macros/PlotCoinXsecFromTable_wTier.C("groups/.../grp_pass4_....list","results","settings",false,true,true)'
```

---

## 3) Single mode vs Group mode (important!)
This macro supports two “semantic” modes because `TableCoinXsec.C` writes slightly different meaning into the `curve_label` field.

### 3.1 Single mode (manifest)
- `z` is binned in the table (`z_lo`, `z_hi`, `z_bin` are meaningful).
- Curves represent **z-bins**.
- Legend entries are displayed as `z ∈ [z_lo, z_hi]`.

### 3.2 Group mode (group `.list`)
- The group CSV typically integrates over `z_bin` (or uses only one z-bin).
- Each curve corresponds to a **setting**; the `curve_label` usually encodes `z0p36` etc.
- Legend is formatted cleanly as `z = 0.36` (via parsing of the label).

---

## 4) The physics/mathematics being performed

### 4.1 What is \(\phi_{pq}\) and why fit cosines?
In semi-inclusive DIS (and in SIDIS-like coincidence analyses), the unpolarized cross section often exhibits azimuthal structure that can be expanded in harmonics of \(\phi\):
- \(\cos(\phi)\) term: commonly associated with kinematic effects and twist‑3 / transverse‑momentum correlations (context-dependent).
- \(\cos(2\phi)\) term: can be related to transverse momentum effects and spin‑momentum correlations in certain formalisms (context-dependent).

This macro does **not** claim a specific mechanism; it provides a *compact, analysis-friendly summary* of the measured \(\phi\)-shape using the minimal harmonic basis:
\[
\sigma(\phi) = M_0(1 + A\cos\phi + B\cos2\phi)
\]
so you can later:
- compare modulations across \(p_T\) bins,
- see whether \(A,B\) are consistent with expectations,
- and test stability vs selection, subtraction, and MC assumptions.

### 4.2 Why tiering exists (numerical conditioning and physics honesty)
A three-parameter fit requires adequate \(\phi\) coverage and enough points. In practice:
- bins near acceptance edges may have missing sectors,
- some \(p_T\) bins have sparse statistics,
- a single point with a tiny error bar can dominate the fit.

If you force a 3‑parameter model onto poor coverage, you get:
- unstable covariance matrices,
- large or “at limit” parameters,
- huge \(\chi^2/\mathrm{ndf}\),
- or physically nonsensical oscillations.

**Tiering** prevents that by fitting only what the data can support.

---

## 5) Tiering: what “tier 0/1/2” mean

### 5.1 Tier definitions
Tiering chooses which functional form to fit:

- **Tier 0**: constant
  \[
  \sigma(\phi) = M_0
  \]
- **Tier 1**: constant + \(\cos\phi\)
  \[
  \sigma(\phi) = M_0\left(1 + A\cos\phi\right)
  \]
- **Tier 2**: constant + \(\cos\phi\) + \(\cos2\phi\)
  \[
  \sigma(\phi) = M_0\left(1 + A\cos\phi + B\cos2\phi\right)
  \]

### 5.2 How the tier is chosen (`ChooseTier`)
Tier selection depends on:
- how many points survive the **fit-only filter** (`npts_fit`)
- the \(\phi\) span covered (`phi_span`)
- a crude **uniformity** measure: how many \(\phi\) sectors are populated (`phi_occ`)

Rules implemented:

1. Start from a tier based on \(\phi\) span:
   - span < `tier_span0_max` (default 1 rad) → tier 0
   - span < `tier_span1_max` (default \(\pi\)) → tier 1
   - else → tier 2

2. Downgrade if too few fit points:
   - tier 2 needs ≥4 points ideally; if <4 → tier 1 or 0
   - tier 1 needs ≥3 points; else → tier 0

3. Downgrade if coverage is not uniform enough:
   - tier 2 requires `phi_occ ≥ phi_occ_min_tier2` (default 4/8 sectors)
   - tier 1 requires `phi_occ ≥ phi_occ_min_tier1` (default 3/8 sectors)

### 5.3 tierOn vs tierOff
- **tierOn = true**  
  Uses tiering and **tries tiers from chosen tier down to 0** until it finds a “valid” fit.

- **tierOn = false**  
  Forces **tier 2 only** (always tries the 3-parameter model).  
  This is useful as a diagnostic reference: you can see what the “full” model would do even when it is not trustworthy.

---

## 6) Filters and guards: why they exist and what they risk

There are two distinct stages of filtering:

### 6.1 Plot-point filter (DEFAULT) — `PassPointDefault`
This controls what appears on the plot.

A point is plotted only if:
- \(\phi\), `xsec`, `xsec_err` are finite and not missing,
- `xsec_err > 0`,
- if enabled: `valid_default == 1` (from `TableCoinXsec.C`)
- if enabled: `xsec > 0` (negative xsec suppressed by default)

Config switches:
- `cfg.use_valid_default` (default true)
- `cfg.allow_negative_xsec` (default false)

**Physics tradeoff:**  
`valid_default` is a conservative pre-filter set by the table maker. If it is too strict, you could hide real structure in low-stat bins. You can loosen this by:
- turning `use_valid_default = false` (not exposed as an argument right now; edit `PlotConfig`)
- or allowing negative xsec for diagnostic plots.

### 6.2 Fit-only filter (stricter) — `BuildFitVectors`
Even if a point is plotted, it may be excluded from the fit.

Fit-only selection requires (per point):
- finite values, `y > 0`, `e > 0`
- `relerr = e/y ≤ cfg.fit_relerr_max` (default 0.60)

**Physics tradeoff:**  
Large-error points carry little statistical weight and can destabilize the fit tier logic (coverage tests) while contributing essentially nothing to parameter constraints. Excluding them makes fits more stable, but:
- if many bins have large uncertainties, you might downgrade to tier 0/1 more often, potentially hiding a real \(\cos2\phi\) modulation.

The solution is not to “force tier 2”; it is to:
- increase statistics,
- adjust binning,
- or revisit upstream cuts/background subtraction.

---

## 7) Robust fitting logic (how we prevent one point from owning the fit)

### 7.1 Error floor (weight regularization)
After selecting fit candidates, the macro applies:
- median of positive errors: `med = MedianPositive(ef)`
- floor: `floor = fit_sigma_floor_frac * med` (default 0.25×median)
- each `ef[i] = max(ef[i], floor)`

This prevents extremely small error bars (often due to low-stat quirks, or overly optimistic error estimates in pathological bins) from giving a single point absurd weight.

### 7.2 “Dominant point” detection
Dominant weight fraction:
\[
f_{\max} = \frac{\max_i(1/\sigma_i^2)}{\sum_i (1/\sigma_i^2)}
\]
implemented in `DominantWeightFrac`.

### 7.3 Robust single-point drop
If after the floor a single point still dominates:
- if `dominant_wfrac > robust_drop_dom_wfrac_gt` (default 0.60)
- and there are at least 4 fit points

then the macro **drops the most influential point** (smallest \(\sigma\)) and refits once.
This is recorded by:
- `ROBUST_DROPPED_PT` flag
- per-point reason `FIT_ROBUST_DROPPED` in the parsed table

### 7.4 Hard reject if still dominated
Even after robust drop, the fit is marked invalid if:
- `dominant_wfrac > max_weight_frac_reject` (default 0.80)

This is a “trust guard”: a fit dominated by one point is typically not physics—it's numerics.

---

## 8) Fit quality evaluation: `ComputeFitValidity`

After each fit attempt, the macro sets `fit_flag_bits` and derives:

- `valid_fit_default`
- `valid_M0_default`
- `valid_A_default`
- `valid_B_default`

### 8.1 FitFlagBits (meaning)
Flags include (bitmask):

- `FIT_STATUS_BAD`  
  Minuit/fit call returned non-zero status.

- `COV_FAILED`, `COV_NOT_POSDEF`  
  Covariance matrix absent or not positive-definite (parameter errors unreliable).

- `NDF_NONPOSITIVE`  
  Fit has no degrees of freedom.

- `CHI2NDF_HUGE` (`chi2/ndf > chi2ndf_max`, default 10)  
  Poor fit or underestimated errors.

- `PROB_TINY` (`prob < prob_min`, default 1e-6)  
  Another “trust” criterion.

- `PHI_UNIFORM_FAIL`  
  Coverage too non-uniform for chosen tier.

- `PARAM_AT_LIMIT`  
  \(A\) or \(B\) at/near constraint boundary (|A|≈1 or |B|≈1). Often indicates insufficient coverage or strong correlations.

- `M0_BAD`, `M0_RELERR_HUGE`  
  \(M_0\) negative/zero or fractional uncertainty too large.

- `A_BAD`, `A_ERR_BAD`, `A_SIG_LOW`  
  \(A\) out of sanity bounds, error nonsensical, or significance too low.

- `B_BAD`, `B_ERR_BAD`, `B_SIG_LOW`  
  same for \(B\).

- `DOMINANT_WEIGHT_PT`  
  the fit is still dominated by one point.

- `NO_FIT_POINTS`  
  not enough points to attempt the fit.

- `ROBUST_DROPPED_PT`  
  robust procedure dropped one point.

### 8.2 Default validity logic
A fit is “default valid” if it avoids “catastrophic” bits plus several trust guards:
- catastrophic: status/covariance failure/nonfinite params/no points/no ndf
- plus: non-posdef covariance, huge chi2/ndf, tiny prob, phi-uniform fail, dominant weight

This produces:
- `valid_fit_default = 1` for fits worth trusting by default
- specialized validity bits for reporting \(M_0\), \(A\), \(B\)

---

## 9) Drawing policy (why some fits are not drawn)

### 9.1 tierOn mode
- Fit curve is drawn only if:
  - `valid_fit_default == 1`
  - OR `cfg.draw_invalid_fits == true` (then drawn dashed)

This prevents “garbage” fits from being visually interpreted as physics.

### 9.2 tierOff mode
- The tier 2 fit is always drawn (solid if valid, dashed if invalid).
- This is intentional: tierOff is a diagnostic reference, not a production result.

---

## 10) Code structure: schematic map

### 10.1 High-level pipeline

```
PlotCoinXsecFromTable_wTier(manifestOrGroupPath, resultsRoot, settingsRoot, tierOn, ...)
 ├─ Determine mode:
 │    ├─ if path ends with ".list" → group mode
 │    └─ else → single mode
 ├─ Build paths:
 │    ├─ input csvPath
 │    ├─ output pngPath
 │    ├─ fitPath
 │    └─ parsedPath
 ├─ ReadTableCSV(csvPath) → rows[]
 ├─ Group rows into curves (after plot filter):
 │    ├─ CurveKey = (pt_bin, z_bin or setting_id, legend, sort_z)
 │    └─ CurvePoints stores x=phi, y=xsec, ey=xsec_err + mapping to original row indices
 ├─ For each pT bin (one canvas pad):
 │    ├─ Draw frame + legend
 │    ├─ For each curve in this pad:
 │    │    ├─ Draw points
 │    │    ├─ FitCurve(curve, tierOn) → FitOut + TF1*
 │    │    ├─ Write fit CSV row (always)
 │    │    └─ Draw TF1 depending on validity + tierOn/tierOff policy
 │    └─ done
 ├─ Write parsed points CSV (loop over ALL rows)
 └─ Save PNG
```

### 10.2 FitCurve internals (most important)
```
FitCurve(curve)
 ├─ BuildFitVectors:
 │    ├─ exclude bad values
 │    ├─ exclude relerr > fit_relerr_max
 │    ├─ apply sigma floor based on median errors
 │    └─ robust drop dominant point if needed
 ├─ compute fit-point phi span, occupancy, dominant weight fraction
 ├─ determine tierStart (ChooseTier) or force tier2 (tierOff)
 ├─ loop over tier = tierStart → ... → minTier:
 │    ├─ build TF1 for that tier
 │    ├─ fit TGraphErrors
 │    ├─ compute validity flags (ComputeFitValidity)
 │    ├─ if first valid fit found: accept and stop
 │    └─ else keep best invalid backup
 ├─ build human-readable "attempts" summary string
 └─ return FitOut + TF1* (for drawing)
```

---

## 11) Function-by-function reference (what each does)

### 11.1 Path and string utilities
- **`NormalizeSlashes`**  
  Normalizes filesystem paths (removes duplicate slashes etc.) so concatenation is robust.

- **`Dirname`, `Basename`, `EndsWith`, `StripExtension`**  
  Used to form `setting_id` / `group_id` and build output paths.

- **`MakeSettingIdFromManifestPath`**  
  Converts manifest path into a consistent `setting_id` relative to `settingsRoot`.  
  This matches the directory conventions used by `TableCoinXsec.C`.

- **`MkdirP`**  
  Creates output directories (`results/<id>/{tables,PNGs}`).

### 11.2 CSV parsing helpers
- **`SplitCSVLine`**  
  Minimal CSV splitter supporting quoted fields.

- **`ToInt`, `ToLL`, `ToD`**  
  Safe conversion helpers.

- **`ReadTableCSV`**  
  Reads the `xsec_phipq_z_pt_overlayed_*.csv` and fills a vector of `Row` structs.

### 11.3 Kinematic helpers
- **`WrapPhi`**  
  Wraps \(\phi\) to \([0,2\pi)\). This prevents discontinuity issues near 0/2π.

- **`ParseZSetting`**  
  Extracts a numeric z from strings like `z0p36` → `0.36` (group legend formatting).

### 11.4 Plot filtering
- **`PassPointDefault`**  
  Implements plot-point filter using `valid_default` and sign conventions.

- **`PassPointDefaultReason`**  
  Same as above but returns an explicit reason token (used in parsed points CSV).

### 11.5 Fit robustness helpers
- **`MedianPositive`**  
  Computes median of positive errors for the sigma-floor logic.

- **`PhiOccupancy`**  
  Counts how many angular sectors contain at least one fit point.

- **`DominantWeightFrac`**  
  Computes the maximum weight fraction \(f_{\max}\) to detect domination.

- **`BuildFitVectors`**  
  Applies fit-only filter, sigma floor, and robust point drop.  
  Returns per-point inclusion/exclusion tokens.

### 11.6 Tiering and validity
- **`ChooseTier`**  
  Picks tier (0/1/2) based on span, points, and occupancy.

- **`ComputeFitValidity`**  
  Converts fit output + configuration thresholds into:
  - `fit_flag_bits`
  - `valid_fit_default`, and component validity flags.

- **`AttemptReasonToken`**  
  Converts failure reasons into short tokens for the attempts history string.

### 11.7 Fit + reporting
- **`FitCurve`**  
  Main fit engine that tries tiers, falls back, writes attempt history, returns `FitOut` and a `TF1`.

- **`WriteFitHeader`, `WriteFitRow`**  
  Fit parameter CSV writer.

- **`EvalSigmaModel`**  
  Evaluates the fitted function at a given \(\phi\) for tier 0/1/2.

- **`WriteParsedHeader`, `WriteParsedRow`**  
  Writes the per-table-row audit CSV including fit residual/pull when possible.

### 11.8 Entry point
- **`PlotCoinXsecFromTable_wTier`**  
  Orchestrates mode selection, reading, grouping, plotting, fitting, and writing outputs.

---

## 12) Interpreting the output files

### 12.1 Fit CSV (`fit_parameters_*_<tierTag>.csv`)
Key columns:
- `npts` vs `npts_fit`: how many were plotted vs used in fit
- `phi_span`, `phi_occ`, `dominant_wfrac`: coverage / dominance diagnostics
- `fit_tier`: final tier chosen
- `M0, A, B` and errors; plus significance `A_sig`, `B_sig`
- `chi2_ndf`, `prob`: fit quality
- `fit_flag_bits`: bitmask
- `valid_fit_default` and component validities
- `attempts`: a compact fit history (e.g. `t2:bad(COV); t1:ok_or_best`)

### 12.2 Parsed points CSV (`parsed_dataPoints_*_<tierTag>.csv`)
Use this when someone asks:
- “Why didn’t this point show up?”
- “Why did the fit ignore this bin?”
- “Which point got robust-dropped?”
- “How big are residuals and pulls?”

Useful columns:
- `is_plotted`, `plot_drop_reason`
- `is_used_in_fit`, `fit_excl_reason`
- `y_fit`, `residual`, `pull`
- `attempts`, `note`

---

## 13) FAQ — questions a new analyzer will ask

### Q1: What does “tierOn” buy me compared to “tierOff”?
TierOn is a stability strategy:
- it avoids pretending you measured a \(\cos2\phi\) modulation when your coverage can’t constrain it.
TierOff is a stress test:
- it shows you what the “full” model tries to do everywhere, which is useful for diagnosing pathologies.

### Q2: Could tiering hide real physics (e.g., true \(\cos2\phi\))?
Yes, in a very specific sense:
- If your data truly has a strong \(\cos2\phi\) term but your \(\phi\) coverage or statistics are too weak, tiering may downgrade to tier1 or tier0.
But that’s not “hiding” physics—it's refusing to over-interpret insufficient information.

The right fix is:
- improve statistics,
- adjust binning,
- or improve acceptance coverage / selection consistency.

### Q3: Why are negative cross sections not plotted by default?
A negative extracted cross section usually means:
- background subtraction dominates, or
- the bin is statistically insignificant.
It can still be useful for debugging (so it remains in the tables), but plotting it as “physics” is misleading.

### Q4: Why can a plotted point be excluded from the fit?
Because the fit aims to be **stable** and interpretable. High relative error points:
- contribute little information,
- but can confuse tier selection and coverage metrics.

### Q5: What is the meaning of \(M_0\), \(A\), and \(B\)?
- \(M_0\): overall scale of \(\sigma(\phi)\) in that bin (roughly the \(\phi\)-averaged level).
- \(A\): amplitude of \(\cos\phi\) modulation.
- \(B\): amplitude of \(\cos2\phi\) modulation.
They are dimensionless modulation parameters in this chosen normalized form.

---

## 14) Practical debugging workflow (recommended)

1. Open the PNG:
   - are there missing points near 0 and \(2\pi\)?
   - do some pads have only a few points?

2. Check `fit_parameters_*_<tierTag>.csv`:
   - look at `phi_occ`, `dominant_wfrac`, `fit_tier`, `valid_fit_default`

3. If anything is surprising, open `parsed_dataPoints_*_<tierTag>.csv`:
   - filter by `pt_bin` and `plot_label`
   - inspect `plot_drop_reason` and `fit_excl_reason`
   - look at `pull` distribution (big pulls suggest underestimated errors or wrong model)

---

## 15) Summary (what to tell someone new)
- `TableCoinXsec.C` computes the physics cross section bins and provides conservative `valid_default` flags.
- `PlotCoinXsecFromTable_wTier.C` **visualizes** those bins and performs a **tiered Fourier fit** to extract stable modulation parameters.
- The macro is designed to be **honest about data limitations**:
  - it avoids overfitting when \(\phi\) coverage or statistics are weak,
  - it flags fits that are numerically or statistically untrustworthy,
  - and it exports an audit table explaining every plotted/fit decision.

---
