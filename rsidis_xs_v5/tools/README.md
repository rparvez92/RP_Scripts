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

## 7) Extract Delta-only and Full-cut data yields

`RP_get_good_coin_ev.C` uses the PID-selected Delta-only sample to determine
CTmean, CTsigma, and the common central/random RF windows. It then integrates
those same fixed windows for two matched acceptance tiers:

```text
Delta-only:
  -8 < H_gtr_dp < 8
  -10 < P_gtr_dp < 22

Full-cut:
  Delta-only plus
  -0.15 < H_gtr_th < 0.15
  -0.10 < H_gtr_ph < 0.10
  -0.15 < P_gtr_th < 0.15
  -0.10 < P_gtr_ph < 0.10
```

The corresponding SIMC branches use the same numerical limits:
`hsxptar`/`ssxptar` are restricted to `[-0.15, 0.15]`, and
`hsyptar`/`ssyptar` are restricted to `[-0.10, 0.10]`.

The output reports explicit `RP_Goodcoin_delta` and `RP_Goodcoin_full`
counts and errors. `RP_get_coin_normyield.C` applies the same run
normalization factor to both, producing explicit `_delta` and `_full`
normalized yields. Replay `ransubcoin` and `normyield` comparisons remain
against the Delta-only tier.

Data PID cuts are a common data-only baseline because standard SIMC does not
simulate calorimeter or Cherenkov detector response. The Delta-only and
Full-cut labels refer to geometrical acceptance matched between data and SIMC.

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

The extractor writes the complete catalog, a warning/error-only table, and a
multipage diagnostic PDF under:

```text
results/Tables/RP_extract_simc_info/
results/PDFs/RP_extract_simc_info/
```

Limit a rerun when needed:

```bash
python3 tools/RP_extract_simc_info.py --phase 1 --reaction sidis
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

The catalog uses `nocut`, `delta`, and `full` acceptance tiers. File sizes are
reported in decimal MB (`bytes / 1e6`). Reconstructed `Ngen` and `normfac` are
reported as single values, while the extractor still rejects either branch
when it is not constant across the reconstructed tree.

Kinematic weighted means, RMS values, and recon-minus-SIMC residuals are saved
with explicit `delta` and `full` tier names. The QA PDF groups target/charge
points beneath one label per setting and compares Delta-only with Full-cut
yield and MC-precision diagnostics.

Kinematic availability distinguishes usable branches (`available`), branches
containing only sentinel values such as `-999` (`sentinel_only`), and absent
branches (`missing`). Sentinel-only and missing metrics remain blank, with the
reason printed in the PDF panel. Entry-count diagnostics use at most six
settings per page so Raw/Hist/Recon markers remain distinguishable.

Exclusive PIMINUS production from LH2 is omitted entirely from the SIMC input
manifest and therefore from batch-running and QA catalogs; exclusive PIMINUS
from LD2 remains supported.

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
  results/Tables/RP_get_coin_normyield/<setting>.csv
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
results/Tables/RP_data_mc_compare/
results/PDFs/RP_data_mc_compare/
```

The data calculation rebuilds differential random subtraction from the saved
run-specific CT windows, normalizes each run with its saved `Norm_factor`, and
combines runs by beam-charge fraction. It then subtracts positrons and the
phase/target-specific thickness-scaled dummy yield. `S_dummy` is retained in
every binned output row.

Each SIMC reaction uses `fWeight` directly. SIDIS, rho, delta, and exclusive
are shown individually and summed without fitted scaling. Delta-only and
Full-cut use identical data/SIMC acceptance definitions. Acceptance QA pages
label cut-removed (N-1) distributions explicitly: every Full-cut requirement
except the cut on the plotted variable is active.

Results remain visibly provisional while computer live time, trigger
efficiency, or PID efficiency are provisionally set to one. A sentinel setting
identity is `SKIPPED`; a valid setting with missing data classes or MC
components is `PENDING`.

All per-setting comparison variables use 20 fixed bins with variable-specific
ranges.  PDF pages use the same wide layout, place Delta-only and Full-cut
comparisons side by side, and report `S_dummy` to four decimal places.  CSV
outputs retain the full configured `S_dummy` precision.  SIMC effective-statistics
warnings remain part of SIMC QA but are not propagated into the downstream
Data-to-MC status.
