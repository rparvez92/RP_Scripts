# rate_dependence_v1 — RSIDIS Rate Dependence Diagnostic

This directory contains the **v1** workflow for a **rate-dependence diagnostic** study in the Hall‑C RSIDIS coincidence analysis.

The goal is to check whether the **charge-normalized coincidence yield** (computed with **my coincidence-time random subtraction**) depends on **beam current** and/or correlates with **trigger rates**.

This is a *diagnostic* workflow: we intentionally do **minimal filtering** so we can spot problematic runs/outliers.

The current workflow is based on **Pass0p1** inputs:
- `bigtable/rsidis_bigtable_pass0p1.csv`
- `Pass0p1_ROOTfiles/`
- `Pass0p1_REPORTfiles/`
- `Skimmed_ROOTfiles/` symlinked to the Pass0p1 skim files

---

## Summary of What v1 Produces

For each kinematic **setting** (one `manifest.json`), v1 generates:

1. **Yield vs Current** (signal only)
   - Uses skimmed ROOT files + coincidence-time random subtraction (`CoincidenceRandomSubtraction.h`)
   - Applies charge / livetime / tracking eff / prescale / boil_corr corrections (via bigtable metadata)
   - Fits a **constant** (`pol0`) and reports `C`, `chi2/ndf`, `Prob`.

2. **Yield vs SHMS 3/4 Trigger Rate** (signal only)
   - Uses the **yield already computed** by (1), then parses report files to get SHMS 3/4 trigger rate in kHz
   - Fits a **constant** (`pol0`) and reports fit quality.

3. **HMS_hEL_CLEAN Charge-Normalized Counts vs Current** (signal only)
   - Parses report files to get `HMS_hEL_CLEAN` total counts
   - Charge-normalizes by BCM2_Q from the Pass0p1 bigtable
   - Fits a **constant** (`pol0`) as another current-stability diagnostic.

4. **Yield Ratio vs SHMS 3/4 Trigger Rate**
   - Combines yield and HMS EL-clean diagnostics into a normalized ratio.
   - Fits the normalized ratio vs trigger rate and converts the slope to a tau estimate.

5. **Coincidence Blocking by Run**
   - Uses Pass0p1 replay ROOT files to measure a simple coincidence-blocking ratio from raw coincidence time.

6. **SHMS Hodo 3-of-4 Livetime by Run**
   - Parses Pass0p1 report files to get `P1X`, `P1Y`, `P2X`, and `P2Y` plane rates.
   - Computes the Dave Mack combinatoric/incoherent 3-of-4 hodoscope trigger livetime using a configurable DPR, default `48 ns`.
   - Writes a run-keyed CSV used by `YieldVsCurrent.C`.

7. **Outliers by Current**
   - Scans existing `yield_vs_current_signal.csv` files and writes a combined list of flagged runs.
   - Flags both pull outliers and runs excluded from the fit.

Most per-setting products are saved as **PNG**, **ROOT canvases**, **CSV tables**, and **logs**. Summary macros write combined CSV/PNG products under `results/`.

---

## How v1 Differs from v0

### v0 (previous)
- Supported plotting **multiple categories** per setting: `signal`, `positron`, `dummy`, `positron_dummy`
- Supported comparisons with “bigtable normyield” vs “myCTime” in some modes
- Used more complex overlay logic and multiple output naming conventions during development

### v1 (current)
- **Signal-only** plots (focus on diagnosing primary yield stability vs current)
- Adds **trigger-based diagnostics** using report files:
  - Yield vs SHMS 3/4 trigger rate (kHz)
  - HMS_hEL_CLEAN counts / BCM2_Q vs current with a **constant** stability test
- Standardizes output layout:
  - `PNGs/`, `tables/`, `canvases/`, `logs/`
- Keeps the workflow modular with separate macros:
  - `YieldVsCurrent.C`
  - `YieldVsTrigger.C`
  - `TriggerVsCurrent.C`
  - `YieldRatioVsTrigger.C`
  - `CoincidenceBlockingByRun.C`
  - `OutlierByCurrent.C`
  - `TauHistogram.C`
  - `TauVsRateAndCurrent.C`
- Does **not** draw a top x-axis (we intentionally decided against multi-axis plots for simplicity)

---

## Directory Layout

```
rate_dependence_v1/
  bigtable/
    rsidis_bigtable_pass0p1.csv # master metadata table (current, charge, efficiencies, boil_corr, etc.)

  settings/
    pass4/pi+sidis/LH2/z.../x.../Q.../<setting_id>/
      manifest.json
      runs_signal.txt
      runs_dummy.txt
      runs_positron.txt
      runs_positron_dummy.txt
      run_metadata.csv          # per-run metadata extracted from bigtable

  include/
    CoincidenceRandomSubtraction.h
    (other headers if needed)

  macros/
    YieldVsCurrent.C            # compute yield (myCTime) from ROOT + random subtraction; plot vs current; const fit
    YieldVsTrigger.C            # read YieldVsCurrent CSV + parse report files; plot yield vs SHMS 3/4 rate; const fit
    TriggerVsCurrent.C          # HMS_hEL_CLEAN counts / BCM2_Q vs current; const fit
    YieldRatioVsTrigger.C       # normalized yield/EL-clean ratio vs SHMS 3/4 rate; tau extraction
    CoincidenceBlockingByRun.C  # coincidence-blocking ratio by run from Pass0p1 replay ROOT files
    Hodo3of4LivetimeByRun.C     # SHMS hodo 3-of-4 LT by run from report P1X/P1Y/P2X/P2Y rates
    OutlierByCurrent.C          # combined current-yield outlier/problem-run CSV
    TauHistogram.C              # histogram tau values from yield-ratio outputs
    TauVsRateAndCurrent.C       # tau summaries vs SHMS 3/4 rate and current

  results/
    pass4/pi+sidis/LH2/z.../x.../Q.../<setting_id>/
      PNGs/
      tables/
      canvases/
      logs/

  Skimmed_ROOTfiles    -> Pass0p1 skim files                     (symlink)
  Pass0p1_ROOTfiles    -> Pass0p1 replay ROOT files              (symlink)
  Pass0p1_REPORTfiles  -> Pass0p1 report files                   (symlink; see below)
```

---

## Inputs and Conventions

### Generate `settings/`
The settings tree is generated from the Pass0p1 bigtable. From `rate_dependence_v1/`:

```bash
python3 tools/generate_settings_tree_rate_v1.py \
  --bigtable bigtable/rsidis_bigtable_pass0p1.csv \
  --outdir settings \
  --targets LH2 LD2 \
  --run-types coin
```

This creates `settings/<pass>/<run_type>/<target>/<z>/<x>/<Q2>/<setting_id>/` directories with run lists, `run_metadata.csv`, and `manifest.json`.

### ROOT files (skimmed)
We use **skimmed coin ROOT files** via the symlink:

- `Skimmed_ROOTfiles/skimmed_coin_replay_production_<RUN>_-1.root`
- Tree name: `T`
- Branch naming uses underscore format, e.g. `CTime_ePiCoinTime_ROC2`.

### Report files (trigger rates)
We use report files via symlink `Pass0p1_REPORTfiles`:

- `Pass0p1_REPORTfiles/COIN/PRODUCTION/replay_coin_production_<RUN>_-1.report`

The macros parse lines such as:
- SHMS 3/4 trigger rate:
  - `SHMS 3/4 Trigger Rate         : 864.802 kHz`
- HMS EL_CLEAN counts:
  - `HMS_hEL_CLEAN : 590009    [ 0.327 kHz ]`
- SHMS hodoscope plane rates:
  - `P1X    : 2250064938.000000 [ 1263.920 kHz ] AND between + and - sides of P1X`
  - `P1Y    : ... [ ... kHz ] ...`
  - `P2X    : ... [ ... kHz ] ...`
  - `P2Y    : ... [ ... kHz ] ...`

`Hodo3of4LivetimeByRun.C` parses the `P1X`, `P1Y`, `P2X`, and `P2Y` lines by label, not by line number.

### Pass0p1 replay ROOT files
`CoincidenceBlockingByRun.C` reads full replay ROOT files via:

- `Pass0p1_ROOTfiles/coin_replay_production_<RUN>_-1.root`

### Bigtable metadata
Per run metadata is read from `settings/.../run_metadata.csv` (generated from the bigtable).
Columns used include:
- `BCM2_I` (current)
- `BCM2_Q` (charge, mC)
- `comp_livetime`
- `h_esing_Eff` (tracking efficiency)
- `boil_corr` (multiplicative correction; if `-999`, flagged)
- `ps5`, `ps6` (prescale; the generator selects which is valid)

### Cuts (signal yield)
The signal selection in `YieldVsCurrent.C` uses:
- HMS acceptance: `-8 < H_gtr_dp < 8`
- HMS electron ID: `H_cal_etottracknorm > 0.7`, `H_cer_npeSum > 2.0`
- SHMS acceptance: `-10 < P_gtr_dp < 22`
- SHMS pion/hadron selection: `P_cal_etottracknorm < 0.8`

Coincidence-time random subtraction uses:
- branch: `CTime_ePiCoinTime_ROC2`
- peak-search window: `5 ns < CTime_ePiCoinTime_ROC2 < 70 ns`
- signal window: peak center ± 1 ns
- random windows: side windows with the same half-width, shifted by ±2 and ±3 RF periods with RF period = 4 ns

---

## Outputs (per setting)

All results are written under:

```
results/<same path as settings>/<setting_id>/
  PNGs/
  tables/
  canvases/
  logs/
```

### 1) YieldVsCurrent.C
- PNG:
  - `PNGs/yield_vs_current.png`
- ROOT canvas:
  - `canvases/yield_vs_current.root`
- CSV:
  - `tables/yield_vs_current_signal.csv`
- Log:
  - `logs/YieldVsCurrent.log`

The CSV includes run-by-run values such as:
- `BCM2_I`, `yield`, `yield_err`, `Nnet`, `Nnet_err`, `norm_factor`, etc.
- staged yield columns:
  - `yield_no_CB_corr`
  - `yield_CB_corr`
  - `yield_no_hodo3of4_corr`
  - `yield_hodo3of4_corr`
- correction inputs:
  - `coincidence_blocking_ratio`
  - `hodo3of4_LT`
  - `DPR_ns`
  - `P1X_rate_kHz`, `P1Y_rate_kHz`, `P2X_rate_kHz`, `P2Y_rate_kHz`
and flags `status` for bad/missing values.

The standard `yield` and `yield_err` columns are the final corrected values, currently equal to `yield_hodo3of4_corr` and `yield_hodo3of4_corr_err`, so downstream macros continue to work without schema changes.

Includes constant-fit summary (C, chi2/ndf, Prob) in plot and/or CSV summary lines.

### 2) YieldVsTrigger.C
- PNG:
  - `PNGs/yield_vs_trigger_shms34.png`
- ROOT canvas:
  - `canvases/yield_vs_trigger_shms34.root`
- CSV:
  - `tables/yield_vs_trigger_shms34.csv`
- Log:
  - `logs/YieldVsTrigger.log`

### 3) TriggerVsCurrent.C
- PNG:
  - `PNGs/hms_elclean_chargeNorm_vs_current.png`
- ROOT canvas:
  - `canvases/hms_elclean_chargeNorm_vs_current.root`
- CSV:
  - `tables/hms_elclean_chargeNorm_vs_current.csv`
- Log:
  - `logs/TriggerVsCurrent.log`

This diagnostic uses `HMS_hEL_CLEAN` total counts divided by BCM2_Q, with Poisson uncertainty `sqrt(N)/Q`, plotted against BCM2 current and fit with a constant.

### 4) YieldRatioVsTrigger.C
- PNG:
  - `PNGs/yield_ratio_vs_trigger_shms34.png`
- ROOT canvas:
  - `canvases/yield_ratio_vs_trigger_shms34.root`
- CSV:
  - `tables/yield_ratio_vs_trigger_shms34.csv`
- Log:
  - `logs/YieldRatioVsTrigger.log`

This macro reads `yield_vs_trigger_shms34.csv` and `hms_elclean_chargeNorm_vs_current.csv`.
It builds a raw ratio, normalizes it by the zero-rate intercept, then fits the normalized ratio vs SHMS 3/4 rate to extract a tau estimate.

### 5) CoincidenceBlockingByRun.C
- CSV:
  - `tables/coincidence_blocking_by_run.csv`
- PNG:
  - `PNGs/coincidence_blocking_ratio_vs_run.png`
- ROOT canvas:
  - `canvases/coincidence_blocking_ratio_vs_run.root`
- Log:
  - `logs/CoincidenceBlockingByRun.log`

The blocking ratio is:

```text
coincidence_blocking_ratio = ctime_good_withstart / ctime_raw_withstart
```

where `ctime_raw_withstart` uses raw coincidence time in `[-400,400] ns` plus HMS/SHMS good-start-time cuts, and `ctime_good_withstart` additionally requires `5 < CTime.CoinTime_RAW_ROC2 < 70`.

### 6) Hodo3of4LivetimeByRun.C
- CSV:
  - `tables/hodo3of4_livetime_by_run.csv`
- Log:
  - `logs/Hodo3of4LivetimeByRun.log`

The macro computes the SHMS hodoscope 3-of-4 trigger electronic livetime using the incoherent combinatoric approximation:

```text
D_i = R_i * DPR
L_i = 1 - D_i

LT_3of4 =
  L1*L2*L3*L4
  + D1*L2*L3*L4
  + L1*D2*L3*L4
  + L1*L2*D3*L4
  + L1*L2*L3*D4
```

where `R_i` is the report rate for `P1X`, `P1Y`, `P2X`, or `P2Y`, converted from kHz to Hz. The DPR argument is in ns and defaults to `48.0`.

Single-setting run example:

```bash
root -l -b -q 'macros/Hodo3of4LivetimeByRun.C("settings/pass4/pi+sidis/LH2/z0p36/x0p25/Q23p3/hmsPneg1p531_shmsP2p615_hmsTh29p045_shmsTh7p865_thpq2/manifest.json","results",48.0)'
```

The output CSV is keyed by run number and is read by `YieldVsCurrent.C`.

### 7) OutlierByCurrent.C
- Combined CSV:
  - `results/outliers/outliers_by_current.csv`

Run from `rate_dependence_v1/`:

```bash
root -l -b -q macros/OutlierByCurrent.C
```

It scans all `results/pass4/**/tables/yield_vs_current_signal.csv` and `results/pass5/**/tables/yield_vs_current_signal.csv`.
It reports:
- `outlier_by_pull`: valid fitted runs with `abs_pull > 3`
- `excluded_from_fit`: runs excluded from the fit or with invalid yield/yield_err inputs

where:

```text
abs_pull = |(yield - fit_const) / yield_err|
```

### 8) Tau summary macros
- `TauHistogram.C`
  - reads all `yield_ratio_vs_trigger_shms34.csv` files
  - writes `results/tau/tau_hist_all.png`
- `TauVsRateAndCurrent.C`
  - reads all yield-ratio outputs
  - writes `results/tau/tau_vs_shms34.png`
  - writes `results/tau/tau_vs_current.png`

---

## Batch Running (cdaq machine, tcsh)

These commands run each macro over **all settings** by iterating over every `settings/**/manifest.json`.

**Important:** tcsh quoting is picky. Use the exact commands below.

### A) Run CoincidenceBlockingByRun for all settings
From `rate_dependence_v1/`:

```tcsh
set RESULTS = "$cwd/results"

foreach mf (`find settings -name manifest.json | sort`)
  set rel = `echo "$mf" | sed 's|^settings/||'`
  set rel_dir = `echo "$rel" | sed 's|/manifest.json$||'`

  mkdir -p "$RESULTS/$rel_dir/logs"

  echo "RUN CoincidenceBlockingByRun: $mf"
  root -l -b -q 'macros/CoincidenceBlockingByRun.C("'"$cwd/$mf"'","'"$RESULTS"'")' \
    >&! "$RESULTS/$rel_dir/logs/CoincidenceBlockingByRun.batch.log"
end
```

### B) Run Hodo3of4LivetimeByRun for all settings
The third argument is DPR in ns. If omitted, the macro uses the default `48.0`.

```tcsh
set RESULTS = "$cwd/results"

foreach mf (`find settings -name manifest.json | sort`)
  set rel = `echo "$mf" | sed 's|^settings/||'`
  set rel_dir = `echo "$rel" | sed 's|/manifest.json$||'`

  mkdir -p "$RESULTS/$rel_dir/logs"

  echo "RUN Hodo3of4LivetimeByRun: $mf"
  root -l -b -q 'macros/Hodo3of4LivetimeByRun.C("'"$cwd/$mf"'","'"$RESULTS"'",50.0)' \
    >&! "$RESULTS/$rel_dir/logs/Hodo3of4LivetimeByRun.batch.log"
end
```

### C) Run YieldVsCurrent for all settings
Requires:
- `tables/coincidence_blocking_by_run.csv` from `CoincidenceBlockingByRun.C`
- `tables/hodo3of4_livetime_by_run.csv` from `Hodo3of4LivetimeByRun.C`

```tcsh
set RESULTS = "$cwd/results"

foreach mf (`find settings -name manifest.json | sort`)
  set rel = `echo "$mf" | sed 's|^settings/||'`
  set rel_dir = `echo "$rel" | sed 's|/manifest.json$||'`

  mkdir -p "$RESULTS/$rel_dir/logs"

  echo "RUN YieldVsCurrent: $mf"
  root -l -b -q 'macros/YieldVsCurrent.C("'"$cwd/$mf"'","'"$RESULTS"'")' \
    >&! "$RESULTS/$rel_dir/logs/YieldVsCurrent.batch.log"
end
```

### D) Run YieldVsTrigger for all settings
Requires that `YieldVsCurrent.C` has already produced `tables/yield_vs_current_signal.csv`.

```tcsh
set RESULTS = "$cwd/results"

foreach mf (`find settings -name manifest.json | sort`)
  set rel = `echo "$mf" | sed 's|^settings/||'`
  set rel_dir = `echo "$rel" | sed 's|/manifest.json$||'`

  mkdir -p "$RESULTS/$rel_dir/logs"

  echo "RUN YieldVsTrigger: $mf"
  root -l -b -q 'macros/YieldVsTrigger.C("'"$cwd/$mf"'","'"$RESULTS"'")' \
    >&! "$RESULTS/$rel_dir/logs/YieldVsTrigger.batch.log"
end
```

### E) Run TriggerVsCurrent for all settings
```tcsh
set RESULTS = "$cwd/results"

foreach mf (`find settings -name manifest.json | sort`)
  set rel = `echo "$mf" | sed 's|^settings/||'`
  set rel_dir = `echo "$rel" | sed 's|/manifest.json$||'`

  mkdir -p "$RESULTS/$rel_dir/logs"

  echo "RUN TriggerVsCurrent: $mf"
  root -l -b -q 'macros/TriggerVsCurrent.C("'"$cwd/$mf"'","'"$RESULTS"'")' \
    >&! "$RESULTS/$rel_dir/logs/TriggerVsCurrent.batch.log"
end
```

### F) Run YieldRatioVsTrigger for all settings
Requires `YieldVsTrigger.C` and `TriggerVsCurrent.C` outputs.

```tcsh
set RESULTS = "$cwd/results"

foreach mf (`find settings -name manifest.json | sort`)
  set rel = `echo "$mf" | sed 's|^settings/||'`
  set rel_dir = `echo "$rel" | sed 's|/manifest.json$||'`

  mkdir -p "$RESULTS/$rel_dir/logs"

  echo "RUN YieldRatioVsTrigger: $mf"
  root -l -b -q 'macros/YieldRatioVsTrigger.C("'"$cwd/$mf"'","'"$RESULTS"'")' \
    >&! "$RESULTS/$rel_dir/logs/YieldRatioVsTrigger.batch.log"
end
```

### G) Run summary macros
From `rate_dependence_v1/`:

```bash
root -l -b -q macros/OutlierByCurrent.C
root -l -b -q macros/TauHistogram.C
root -l -b -q macros/TauVsRateAndCurrent.C
```

### Optional: run in a screen session
```tcsh
screen -S rateDep
# paste one of the batch blocks above
```

Detach: `Ctrl-a d`  
Exit screen after job: `exit`

---

## Common Pitfalls / Troubleshooting

### “Unmatched '\"'” or quote errors in tcsh
Use the exact quoting style shown above:
```tcsh
root -l -b -q 'macros/Macro.C("'"$cwd/$mf"'","'"$RESULTS"'")'
```
This is tcsh-safe and avoids nested quote issues.

### Missing report file or missing trigger lines
- The macros will mark runs as `MISSING_REPORT` or `MISSING_*_LINE` in the CSV and log.
- Verify your symlink `Pass0p1_REPORTfiles` points to the correct location and contains:
  `COIN/PRODUCTION/replay_coin_production_<RUN>_-1.report`

### Missing replay ROOT file for coincidence blocking
- `CoincidenceBlockingByRun.C` marks runs as `NOFILE` if the replay ROOT file is missing.
- Verify `Pass0p1_ROOTfiles/coin_replay_production_<RUN>_-1.root` exists or that the symlink points to the correct cdaq location.

### Missing hodo 3-of-4 livetime CSV
- `YieldVsCurrent.C` now requires `tables/hodo3of4_livetime_by_run.csv`.
- Run `Hodo3of4LivetimeByRun.C` before `YieldVsCurrent.C`.
- If a run has missing report lines, it will be marked in `hodo3of4_livetime_by_run.csv` and then excluded from the corrected yield fit by `YieldVsCurrent.C`.

### `boil_corr = -999`
- v1 flags `-999` values in logs and runtime warnings (diagnostic intent).
- These runs are not silently removed; they are labeled so you can investigate.

### Interactive zooming
Open the saved ROOT canvas:
```bash
root -l results/.../<setting_id>/canvases/yield_vs_current.root
```
Then in the GUI: right-click axis → “UnZoom” or use the toolbar “Reset”.

---

## Notes on Analysis Strategy

- The constant fit in Yield vs Current / Yield vs Trigger is a direct test of stability:
  - If the yield is independent of current, the constant fit should describe data well (reasonable chi2/ndf and Prob).
- TriggerVsCurrent checks whether charge-normalized HMS EL-clean counts are stable vs current.
- YieldRatioVsTrigger is the tau-extraction diagnostic built from the yield and HMS EL-clean summaries.
- OutlierByCurrent is intended as a review aid, not an automatic run-removal list.
  - `outlier_by_pull` means statistically far from the setting's constant yield fit.
  - `excluded_from_fit` means the run did not have valid inputs for the fit and should be investigated separately.

---

## Maintainers / Contact

This v1 workflow is maintained by RSIDIS analyzers in the RP_Scripts repository.
For questions, contact the current analyzer responsible for rate dependence studies.
