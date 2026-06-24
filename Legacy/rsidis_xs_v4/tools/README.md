# tools/ (v4)

Kept **file naming consistent with earlier versions** where practical.

- `generate_settings_tree_v2_1.py`:
  Generates the settings tree + per-setting manifests from the bigtable.
  (Ported for rsidis_xs_v4 paths/wording; file name kept for compatibility.)

- `generate_settings_tree_v4.py`:
  Alias copy of the above for clarity.

- `build_sigma_model.py`:
  Generates `sigma_model.txt` inside each setting directory by calling `calc_semi_xsec`
  and parsing `sigsemi`.

- `batch_driver_v2.py`:
  Walks all `manifest.json` under settings root and runs a ROOT macro per setting,
  mirroring output into a parallel results directory.
  (File name kept for compatibility.)

- `make_symlinks.sh`:
  Creates/updates symlinks to large data directories.

- `validate_output_schema.py`:
  Checks that your extracted CSV table contains the required schema columns.
