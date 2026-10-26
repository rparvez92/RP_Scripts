# tools/ (v5)

This folder contains the scripts used to generate the v5 settings tree, build per-setting sigma models, create overlay group files, and batch-run ROOT macros.

## 0) Directory layout assumed

From the **rsidis_xs_v5/** top level:

- `settings/` : generated settings tree (per setting leaf directories)
- `groups/`   : overlay-group definition files (e.g. z-overlay)
- `results/`  : output from ROOT macros (written by batch_driver)

You can create these once:
```bash
mkdir -p settings groups results
```

---

## 1) Generate the settings tree

`generate_settings_tree_v5.py` reads the bigtable CSV and creates a directory tree like:
`settings/<pass>/<run_type>/<target>/<z>/<x>/<q2>/<setting_id>/`

Each leaf directory contains:
- `runs_signal.txt`
- `runs_dummy.txt`
- `runs_positron.txt`
- `runs_positron_dummy.txt`
- `run_metadata.csv`
- `manifest.json`
- `README.txt`

Example:
```bash
python3 tools/generate_settings_tree_v5.py \
  --bigtable /path/to/rsidis_bigtable_pass0.csv \
  --outdir   settings \
  --targets  LH2 LD2
```

---

## 2) Build sigma_model.txt for each setting (optional / if your workflow uses it)

`build_sigma_model.py` scans `settings/` for `manifest.json` and runs an external calculator executable per setting.

Example:
```bash
python3 tools/build_sigma_model.py \
  --settings_root settings \
  --calc_exe /path/to/sigma_calc_executable \
  --nphibin 1
```

Use `--dry_run` to preview commands.

---

## 3) Create a z-overlay group file (NEW in v5)

`make_groups_z_overlay.py` writes a group file that lists 3 setting manifests (z0p36, z0p50, z0p67) to be overlaid on the same plot.

Output file is created at:
`groups/<pass>/<run_type>/<target>/<x>/<q2>/<tpq_label>/grp_<TARGET>_zOv_<x>_Q<q2>_<tpq_label>.list`

Example (your case):
```bash
python3 tools/make_groups_z_overlay.py \
  --settings-root settings \
  --groups-root   groups \
  --pass          pass4 \
  --run-type      "pi+sidis" \
  --target        LH2 \
  --x             x0p25 \
  --q2            q23p3 \
  --tpq           2.0 \
  --z-dirs        z0p36 z0p5 z0p67 \
  --require-substr hmsPneg1p531 hmsTh29p045 shmsTh7p865
```

This produces a file like:
`groups/pass4/pi+sidis/LH2/x0p25/q23p3/tpq2p0/grp_LH2_zOv_x0p25_Q23p3_tpq2p0.list`

Notes:
- If multiple settings exist under a given z-dir, you MUST disambiguate using `--require-substr`.
- For readability, z0p5 is labeled as z0p50 in the group file.

---

## 4) Batch-run a ROOT macro over all settings (per-setting mode)

`batch_driver_v2.py` walks the `settings/` tree and runs a ROOT macro per manifest.

Example:
```bash
python3 tools/batch_driver_v2.py \
  --settings settings \
  --results  results \
  --root-macro macros/PlotCoinXsec.C \
  --dry-run
```

Remove `--dry-run` to execute.

---

## 5) Structure local SIMC files for `TableCoinXsec.C`

`structure_SimROOTfiles.py` is for local running when the SIMC ROOT files and
their `.hist` files live outside the project, for example on a Samsung T7.

`TableCoinXsec.C` expects:
```text
Pass0_SimROOTfiles/<setting-id>/<sim-stem>.root
Pass0_SimROOTfiles/<setting-id>/<sim-stem>.hist
```

On the local T7 setup, the payloads are flat:
```text
/Volumes/T7/RSIDIS/Pass0_SimROOTfiles/*.root
/Volumes/T7/RSIDIS/Pass0_SimOUTfiles/*.hist
```

The script scans `settings/**/manifest.json`, derives the expected SIMC stem,
and creates project-local symlinks. It does not copy the large ROOT files.

Preview:
```bash
python3 tools/structure_SimROOTfiles.py
```

Create symlinks:
```bash
python3 tools/structure_SimROOTfiles.py --execute
```

If stale or wrong symlinks already exist:
```bash
python3 tools/structure_SimROOTfiles.py --execute --replace
```

---

## 6) Validate output CSV schema (optional)

After your macro writes an `xsec_phipq.csv`, you can verify it contains required columns:

```bash
python3 tools/validate_output_schema.py results/.../xsec_phipq.csv
```

---

## 7) Extract Pass1 good-coin yields

`RP_get_good_coin_ev.C` uses the Baseline PID + Delta-only acceptance sample to
determine CTmean, CTsigma, the central RF window, and ten random RF windows
(five on each side). The RF period is read from the input bigtable's
`RF_signal` column.

```text
Baseline PID + Delta-only acceptance:
  -8 < H_gtr_dp < 8
  -10 < P_gtr_dp < 22
```

The output reports one `RP_Goodcoin` count and error. Pass1 PDFs are written
under `results/Pass1/PDFs/RP_get_good_coin_ev/`; per-setting, summary, and
problematic CSVs are written under
`results/Pass1/Tables/RP_get_good_coin_ev/`. Filenames use the Pass1
setting-ID schema.

Data PID cuts are a common data-only baseline because standard SIMC does not
simulate calorimeter or Cherenkov detector response. Delta-only refers to the
spectrometer momentum acceptance.

Normalize one completed Pass1 setting with:

```bash
root -l -b -q 'macros/RP_get_coin_normyield.C("results/Pass1/Tables/RP_get_good_coin_ev/<setting>.csv")'
```

Or process all available settings and rebuild the summary/problematic tables:

```bash
python3 tools/run_coin_normyield_batch.py
```

This stage uses each run's bigtable `comp_livetime`; trigger and PID
efficiencies remain provisional constants of 1.0. Outputs use the same setting
ID filenames under `results/Pass1/{PDFs,Tables}/RP_get_coin_normyield/`.

---

## 8) Generate and run the RP SIMC input manifest

From `rsidis_xs_v5/`, preview the inputs derived from the settings tree:

```bash
python3 tools/generate_simc_inputs.py --dry-run
```

Generate the `.inp` files and manifest in the macOS SIMC checkout:

```bash
python3 tools/generate_simc_inputs.py
```

Preview all manifest jobs with a short validation run size:

```bash
python3 tools/run_simc_batch.py --dry-run --ngen 100
```

Run pending jobs serially, skipping outputs already complete for the requested
event count:

```bash
python3 tools/run_simc_batch.py --run --ngen 100
```

The batch runner writes outputs under
`/Volumes/T7/RSIDIS/Phase<1|2>/Simulation/{outfiles,runout,ROOTfiles}/<run-type>/`
and atomically updates `infiles/RP_Simc/simc_batch_status.csv`. Use filters such
as `--phase`, `--reaction`, `--run-type`, `--target`, or `--limit` for a smaller
batch. Use `--overwrite` only when intentionally replacing incomplete or
mismatched output sets.

---

## 9) Extract SIMC metadata, normalization, and QA

From `rsidis_xs_v5/`, inspect the full manifest inventory and the simulation
products currently available on T7:

```bash
python3 tools/RP_extract_simc_info.py
```

The Pass1 job manifest is read from
`JLab_upstream_forks/simc_gfortran/util/batch_job/simc_job_list_pass1.csv`.
Products are resolved under
`/Volumes/T7/RSIDIS/Pass1/Simulation/{ROOTfiles,outfiles,runout}/coin/`.

The extractor writes the complete catalog, a warning/error-only table, and a
multipage diagnostic PDF under:

```text
results/Pass1/Tables/RP_extract_simc_info/
results/Pass1/PDFs/RP_extract_simc_info/
```

Limit a rerun when needed:

```bash
python3 tools/RP_extract_simc_info.py --phase phase1 --reaction sidis
```

Use `--no-pdf` for a CSV-only validation. Missing future products remain
`PENDING`, and intentionally unsupported manifest rows remain `SKIPPED`.
Warnings do not fail the command; a structural `ERROR` produces exit status 2.

For selected reconstructed events, the normalized simulation yield and its
Monte Carlo uncertainty are:

```text
Y_SIMC = sum(fWeight)
sigma_MC = sqrt(sum(fWeight^2))
```

`fWeight` already includes the SIMC normalization factor divided by the
generated-event count and must not be normalized a second time.

The catalog reports no-cut bookkeeping and one untiered baseline-acceptance
result using `-8 < hsdelta < 8` and `-10 < ssdelta < 22`. File sizes are
reported in decimal MB (`bytes / 1e6`). Reconstructed `Ngen` and `normfac` are
reported as single values, while the extractor still rejects either branch
when it is not constant across the reconstructed tree. Kinematic weighted
means, RMS values, and recon-minus-SIMC residuals use untiered column names.

Kinematic availability distinguishes usable branches (`available`), branches
containing only sentinel values such as `-999` (`sentinel_only`), and absent
branches (`missing`). Sentinel-only and missing metrics remain blank, with the
reason printed in the PDF panel. Entry-count diagnostics use at most six
settings per page so Raw/Hist/Recon markers remain distinguishable.

The Pass1 manifest contains exclusive `pip` production only; no exclusive
`pim` points are synthesized by the extractor.

For delta and exclusive samples, the QA layer derives common original and
reconstructed coordinates when their reaction-native `xbj`, `z`, or `pt2`
branches are sentinel-only:

```text
xbj = Q2 / (2 Mp nu)
z   = phad / nu
pt2 = (phad sin(thetapq))^2
```

These coordinates support uniform reconstruction QA and do not imply that the
variables have the same reaction interpretation as SIDIS. Missing mass is
reported from `missmass` and `missmass_recon`.

---

## 10) Build setting-wise data-to-MC comparisons

Process one setting from its target normalized-yield CSV:

```bash
python3 tools/RP_data_mc_compare.py \
  results/Pass1/Tables/RP_get_coin_normyield/<setting>.csv
```

Process every discovered LH2/LD2 setting:

```bash
python3 tools/RP_data_mc_compare.py --all
python3 tools/RP_data_mc_compare_Summary.py
```

Use `--bigtable-dir <path>` when the setting-wise bigtable leaf hierarchy is
not under `rsidis_xs_v5/bigtable`.  The comparison joins `BCM2_I` from the
matching leaf by run number for the Signal normalized-yield-versus-current
diagnostic.

Outputs are written under:

```text
results/Pass1/Tables/RP_data_mc_compare/
results/Pass1/PDFs/RP_data_mc_compare/
```

The data calculation rebuilds differential random subtraction from the saved
run-specific CT windows, normalizes each run with its saved `Norm_factor`, and
combines runs by beam-charge fraction. It then subtracts positrons and the
phase/target-specific thickness-scaled dummy yield. `S_dummy` is retained in
every binned output row.

Each SIMC reaction uses `fWeight` directly. SIDIS, rho, and delta are always
shown and summed without fitted scaling; exclusive is included for `pip` and
is explicitly zero for `pim`, where that production channel is physically
inapplicable. Data use the current
baseline PID and momentum-acceptance cuts; SIMC uses the matching momentum
acceptance. Acceptance QA pages remove only the plotted HMS or SHMS momentum
cut.

Results remain visibly provisional while trigger or PID efficiency is set to
one. Computer livetime comes from the Pass1 bigtable. A sentinel setting
identity is `SKIPPED`. An unavailable positron or dummy correction class is
represented by a zero histogram and identified with `ZERO_CORRECTION_CLASS`
instead of leaving the setting `PENDING`; a missing required SIMC component is
an error.

In `RP_get_good_coin_ev`, a finite negative random-subtracted `RP_Goodcoin` is
floored to zero for positron runs only. `RP_Goodcoin_err` is also set to zero
for those floored runs, so both zeroes propagate to `RP_Normyield` and
`RP_Normyield_err`; electron runs are not floored.

All per-setting 1D comparison variables use 15 fixed bins with variable-specific
ranges. `P_t` spans 0--0.6 GeV. The PDF places a 5-`P_t`-by-15-`phipq` polar
`P_t`-versus-phi page immediately after the normalized-yield correction page;
corrected data and total MC share one color scale, while Data/MC is shown only
in bins with positive MC support. The page also reports the setting's nominal
`z` and `thpq`. Missing mass is labeled `M_X` and spans 0--4 GeV, with reconstructed negative values counted
as underflow. Data or MC flow above 2% receives a `*_FLOW_GT_2PCT_*` warning.
PDF pages use the whole wide canvas for the baseline comparison and
report `S_dummy` to four decimal places. CSV
outputs retain the full configured `S_dummy` precision.  SIMC effective-statistics
warnings remain part of SIMC QA but are not propagated into the downstream
Data-to-MC status.

---

## 11) Extract setting-integrated SIDIS cross sections

Build the Pass1 central-model catalog from the current SIMC job manifest:

```bash
python3 tools/RP_build_sidis_model.py
```

Then extract every setting and rebuild the summary tables:

```bash
python3 tools/RP_extract_sidis_xsec.py --all
python3 tools/RP_extract_sidis_xsec_Summary.py
```

Outputs use the Pass1 setting-ID schema under
`results/Pass1/<Tables,PDFs>/RP_extract_sidis_xsec/`. The extraction integrates
the untiered `hdelta` rows from `RP_data_mc_compare`, corresponding to Baseline
PID + Delta-only acceptance. No full-cut or tier suffix is used. The yield
correction is `C_Y = Y_Data/Y_MC`, where `Y_MC` is the sum of SIDIS, rho,
delta, and the physically applicable exclusive contribution. The same
correction scales the central SIDIS cross-section and multiplicity models.

---

## 12) Build the setting-integrated SIDIS overview

```bash
python3 tools/RP_plot_sidis_setting_overview.py
```

The overview joins the Pass1 extraction summary to the setting-ID model catalog
and writes its catalog and two multipage PDFs under
`results/Pass1/<Tables,PDFs>/RP_sidis_setting_overview/`. It produces controlled
central-`z` scans and central-`P_t` scans of the extracted and model SIDIS
cross section and multiplicity. The cross-section and multiplicity PDFs place
their respective `z` and `P_t` scans side by side with plot grids. `P_t` is the magnitude
`|p_pi sin(theta_pq)|` in GeV; neither `p_T` nor `P_t^2` is used as the scan
coordinate. The overview axes use fixed ranges `0 <= z <= 1` and
`0 <= P_t <= 0.6 GeV`, with labeled major ticks every 0.1. Overview columns
use the concise `Q2` and `thpq` labels, and PDF
provenance is relative to `RP_Scripts/`. Inputs and outputs contain no full-cut
or tier labels.

---

## 13) Combine nominal-theta-pq settings and study azimuthal coverage

```bash
python3 tools/RP_build_thpq_combination.py
python3 tools/RP_combine_thpq_yields.py
python3 tools/RP_azimuthal_coverage_study.py
```

The first step assigns every available setting to a fixed-kinematics group
defined by phase, pass, pion charge, target, `x`, `Q2`, and shared nominal `z`.
Nominal `thpq` is intentionally excluded from the key. Group IDs use fixed
numeric formatting, for example `1501x0p25Q23p3z0p36`. The audit map, one-row
group summary, and problematic rows are written under
`results/Pass1/Tables/RP_thpq_combination/`.

The second step sums corrected data and reaction-separated SIMC across every
contributing nominal-`thpq` setting, bin by bin, on the shared grid of five
`P_t` bins by 15 `phipq` bins. It propagates independent statistical errors in
quadrature and writes all 75 bins per group under
`results/Pass1/Tables/RP_combine_thpq_yields/`. Provenance columns retain all
contributing setting IDs and nominal `thpq` values.

Coverage is evaluated only after those sums. Complete bin-level and curve-level
CSVs under `results/Pass1/Tables/RP_azimuthal_coverage_study/` report the
independent-bin masks, longest contiguous circular passing range,
harmonic-design rank and condition, and explicit decision reasons. PDFs are
written as `z_*/<combination_group_id>.pdf` directly under the matching
coverage-study PDF directory, replacing the former setting-wise results. Page one places polar Data/MC maps before and
after all coverage requirements side by side; rejected bins are absent from
the after map. Page two compares the combined corrected-data and total-MC
`phipq` projections before and after selection, with Data/MC ratio subpads.
The after projection sums only surviving `P_t`-by-`phipq` cells. Page three
retains the five individual `P_t` pads.

The coverage policy is versioned as `RP_AZIMUTHAL_COVERAGE_V6_THPQ_COMBINED`.
Its defaults require positive corrected data and MC, MC relative uncertainty
at most 30%, data/MC relative uncertainty at most 30%, at least `pi` radians of
contiguous coverage, design rank three, and weighted-design condition at most
100. A data/MC ratio outside `[0.5,1.5]` is a warning rather than a rejection.
All thresholds have command-line options.

Curve rows also report the unweighted and weighted harmonic-design conditions, coefficient
uncertainties, and coefficient correlations. These are diagnostics for studying
fit stability versus accepted angular span and do not yet alter approval.

These steps use the shared analysis binning in `tools/RP_sidis_binning.py`.

---

## 14) Extract and fit SIDIS azimuthal distributions

```bash
python3 tools/RP_sidis_azimuthal_analysis.py
```

This post-combination, coverage-controlled workflow writes cross-section and multiplicity point,
problematic, fit-parameter, and fit-parameter-problematic CSVs under
`results/Pass1/Tables/RP_sidis_azimuthal_analysis/`, plus one multipage PDF per
observable and nominal `z` under matching `PDFs/z_*` directories. It consumes
the combined-`thpq` yields and combined coverage decisions from section 13;
each output row retains the contributing setting IDs and nominal `thpq` values.

Each page contains the five `P_t` bins `[0,0.1)`, `[0.1,0.2)`, `[0.2,0.3)`,
`[0.3,0.4)`, and `[0.4,0.5]` GeV. Fifteen half-open phi bins cover `[0,2 pi)`;
an angle at exactly `2 pi` is canonicalized to zero. The current
`calc_semi_xsec` is evaluated with 15 phi bins at each `P_t`-bin center. Each
extracted point uses the model cross section and multiplicity returned at its
matching phi-bin center.

The fitted forms are `C_0 + C_1 cos(phi) + C_2 cos(2 phi)` for cross section and
`M_0 + M_1 cos(phi) + M_2 cos(2 phi)` for multiplicity. The indices 0, 1, and 2
identify the constant, `cos(phi)`, and `cos(2 phi)` basis coefficients. Only bins in an approved curve's
selected contiguous coverage mask enter the fit. Excluded finite points remain
visible as open gray markers. Every point and fit row records the coverage
binning version, policy version, mask/status, and decision reasons. Missing or
rejected coverage decisions fail closed and cannot silently enter a fit.
