# v4 workflow

## 1) Generate settings
Use `tools/generate_settings_tree_v4.py` to create:
- run lists
- per-setting manifests

## 2) Configure analysis
Create a config JSON per setting in `config/`, based on `analysis_template.json`.

## 3) Extract table
Run extraction to produce:
- `results/<setting_id>/tables/xsec_phipq.csv`
- `results/<setting_id>/tables/metadata.json`

## 4) Plot
Plotting should only read the CSV table (no dependence on ROOT trees).

## 5) Systematics
Introduce PID variants via `cut_tag`, producing multiple tables for comparison.
