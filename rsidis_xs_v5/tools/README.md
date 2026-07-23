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

## 7) Generate and run the RP SIMC input manifest

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
