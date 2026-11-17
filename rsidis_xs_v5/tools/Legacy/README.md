# Archived tools

These scripts were moved from `tools/` on 2026-11-17 after checking script
imports, command invocations, output consumers, and documentation in `RP_Scripts/`.
They are not required by the current Pass1 SIDIS analysis.

| Script | Role and reason for archival |
| --- | --- |
| `RP_plot_sidis_setting_overview.py` | Redundant setting-integrated plots; no downstream analysis consumes its output. |
| `filter_bigtable.py` | Older phase-bigtable splitter; Pass1 uses `tools/filter_pass1_sidis_bigtable.py`. |
| `generate_simc_inputs.py` | Older SIMC production inputs and `infiles/RP_Simc/simc_input_manifest.csv`; current Pass1 extraction/model stages use `JLab_upstream_forks/simc_gfortran/util/batch_job/simc_job_list_pass1.csv`. |
| `generate_settings_tree_v5.py` | Historical settings manifests and run lists for the older ROOT workflow. |
| `build_sigma_model.py` | Historical per-manifest `sigma_model.txt`; current Pass1 builds its catalog with `tools/RP_build_sidis_model.py`. |
| `make_groups_z_overlay.py` | Historical manifest group lists for ROOT z-overlay plots. |
| `batch_driver_v2.py` | Generic historical manifest-based ROOT macro runner. |
| `structure_SimROOTfiles.py` | Pass0 SIMC symlinks for `TableCoinXsec.C`; the existing archived copy was byte-identical. |
| `validate_output_schema.py` | Checks the older `xsec_phipq.csv` schema. |

The older ROOT macros still contain readers for settings run lists, model text
files, and group lists. Archiving the generators does not remove those existing
inputs or results. They remain available here for reproducing the historical
workflow.

`tools/run_simc_batch.py` remains in place as requested. It reads the older
generator's manifest, but does not import or invoke that generator. Regenerating
that older manifest requires the archived `generate_simc_inputs.py`.

Run archived commands from the `rsidis_xs_v5/` project directory, for example:

```bash
python3 tools/Legacy/generate_simc_inputs.py --dry-run
python3 tools/Legacy/structure_SimROOTfiles.py
```

The archived splitter, SIMC generator, and overview plotter retain their original
project-relative defaults after relocation. Their existing tests remain in
`tools/tests/`, with imports updated to `Legacy`.
