# Output schema (v4)

This document defines the canonical table schema used by v4.

## File
Per setting:
- `results/<setting_id>/tables/xsec_phipq.csv`
- `results/<setting_id>/tables/metadata.json`

One row corresponds to one point:
**(pt bin, z bin, phipq bin)**.

## Required columns (MVP)
| Column | Type | Meaning |
|---|---:|---|
| setting_id | string | Unique identifier of the setting |
| cut_tag | string | PID/acceptance cut tag used to produce the row |
| pt_bin | int | pt bin index (0..Npt-1) |
| z_bin | int | z bin index (0..Nz-1) |
| phipq_bin | int | phipq bin index (0..Nphi-1) |
| pt_low, pt_high, pt_center | float | pt bin edges + center |
| z_low, z_high, z_center | float | z bin edges + center |
| phipq_low, phipq_high, phipq_center | float | phipq bin edges + center |
| xsec | float | extracted cross section value |
| xsec_stat | float | statistical uncertainty on xsec |

## Strongly recommended columns (debug + audit)
| Column | Type | Meaning |
|---|---:|---|
| yield_data | float | random-subtracted data yield in this bin |
| yield_dummy | float | dummy contribution |
| yield_pos | float | positron contribution |
| yield_pos_dummy | float | positron dummy contribution |
| yield_sim | float | SIM yield (acceptance-weighted) |
| sigma_model | float | model cross section used for normalization (if applicable) |
| norm_factor | float | overall normalization applied |
| charge | float | charge used |
| eff_total | float | total efficiency correction used |
| run_list_tag | string | identifies run grouping used |

## metadata.json (per setting)
Should include:
- input file lists (data/sim/dummy/positron)
- tree name(s)
- bin edges arrays for pt/z/phipq
- cut thresholds + cut_tag
- analysis git hash/branch, timestamp
