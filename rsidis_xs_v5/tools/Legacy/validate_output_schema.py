#!/usr/bin/env python3
import csv
import sys
from pathlib import Path

REQUIRED = [
  "setting_id","cut_tag","pt_bin","z_bin","phipq_bin",
  "pt_low","pt_high","pt_center",
  "z_low","z_high","z_center",
  "phipq_low","phipq_high","phipq_center",
  "xsec","xsec_stat"
]

def main() -> int:
  if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <xsec_phipq.csv>")
    return 2
  p = Path(sys.argv[1])
  with p.open() as f:
    r = csv.reader(f)
    header = next(r)
  missing = [c for c in REQUIRED if c not in header]
  if missing:
    print("Missing required columns:")
    for c in missing:
      print("  -", c)
    return 1
  print("OK: schema contains all required columns.")
  return 0

if __name__ == "__main__":
  raise SystemExit(main())
