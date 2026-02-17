# rate_dependence_v1 — RSIDIS Rate Dependence Diagnostic (Yield/Current/Trigger)

This directory contains the **v1** workflow for a **rate-dependence diagnostic** study in the Hall‑C RSIDIS coincidence analysis.

The goal is to check whether the **charge-normalized coincidence yield** (computed with **my coincidence-time random subtraction**) depends on **beam current** and/or correlates with **trigger rates**.

This is a *diagnostic* workflow: we intentionally do **minimal filtering** so we can spot problematic runs/outliers.

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

3. **SHMS_EL_CLEAN Trigger Rate vs Current** (signal only)
   - Parses report files to get `SHMS_pEL_CLEAN` rate (kHz)
   - Fits a **line** (`pol1`) and reports slope, slope error, **slope significance** `b/σ_b`, plus `chi2/ndf`, `Prob`.

All three outputs are saved as **PNG** and **ROOT canvases** (interactive), plus **CSV tables** and **logs**.

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
  - SHMS_EL_CLEAN vs current (kHz vs current) with **linear** trend test
- Standardizes output layout:
  - `PNGs/`, `tables/`, `canvases/`, `logs/`
- Keeps the workflow modular with **three separate macros**:
  - `YieldVsCurrent.C`
  - `YieldVsTrigger.C`
  - `TriggerVsCurrent.C`
- Does **not** draw a top x-axis (we intentionally decided against multi-axis plots for simplicity)

---

## Directory Layout

```
rate_dependence_v1/
  bigtable/
    rsidis_bigtable_pass0.csv   # master metadata table (current, charge, efficiencies, boil_corr, etc.)

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
    TriggerVsCurrent.C          # parse report files; plot SHMS_EL_CLEAN vs current; linear fit

  results/
    pass4/pi+sidis/LH2/z.../x.../Q.../<setting_id>/
      PNGs/
      tables/
      canvases/
      logs/

  Skimmed_ROOTfiles -> /work/hallc/c-rsidis/skimfiles/pass0     (symlink)
  Pass0_REPORTfiles -> ...                                       (symlink; see below)
```

---

## Inputs and Conventions

### ROOT files (skimmed)
We use **skimmed coin ROOT files** via the symlink:

- `Skimmed_ROOTfiles/skimmed_coin_replay_production_<RUN>_-1.root`
- Tree name: `T`
- Branch naming uses underscore format, e.g. `CTime_ePiCoinTime_ROC2`.

### Report files (trigger rates)
We use report files via symlink `Pass0_REPORTfiles`:

- `Pass0_REPORTfiles/COIN/PRODUCTION/replay_coin_production_<RUN>_-1.report`

We parse two lines (kHz values):
- SHMS 3/4 trigger rate:
  - `SHMS 3/4 Trigger Rate         : 864.802 kHz`
- SHMS EL_CLEAN rate:
  - `SHMS_pEL_CLEAN :    11827351    [ 6.555 kHz ]`

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
The signal selection uses the same baseline cuts as used in development:
- HMS dp cut
- HMS calorimeter cut
- HMS Cer cut
- SHMS dp cut
- SHMS calorimeter cut
and coincidence time random subtraction applied to all events passing cuts.

(See `YieldVsCurrent.C` for the exact cut string and random subtraction windows.)

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
and flags `status` for bad/missing values.

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
  - `PNGs/shms_elclean_vs_current.png`
- ROOT canvas:
  - `canvases/shms_elclean_vs_current.root`
- CSV:
  - `tables/shms_elclean_vs_current.csv`
- Log:
  - `logs/TriggerVsCurrent.log`

This plot uses a linear fit. The slope significance `b/σ_b` is a quick “how many sigmas away from zero slope” metric.

---

## Batch Running (cdaq machine, tcsh)

These commands run each macro over **all settings** by iterating over every `settings/**/manifest.json`.

**Important:** tcsh quoting is picky. Use the exact commands below.

### A) Run YieldVsCurrent for all settings
From `rate_dependence_v1/`:

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

### B) Run YieldVsTrigger for all settings
(Requires that YieldVsCurrent has already produced `tables/yield_vs_current_signal.csv`.)

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

### C) Run TriggerVsCurrent for all settings
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
- Verify your symlink `Pass0_REPORTfiles` points to the correct location and contains:
  `COIN/PRODUCTION/replay_coin_production_<RUN>_-1.report`

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
- TriggerVsCurrent uses a linear fit because trigger rate is expected to scale with current.
  - Slope significance `b/σ_b` quantifies whether the observed slope differs from zero (in sigma units).

---

## Maintainers / Contact

This v1 workflow is maintained by RSIDIS analyzers in the RP_Scripts repository.
For questions, contact the current analyzer responsible for rate dependence studies.

