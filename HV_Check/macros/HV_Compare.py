#!/usr/bin/env python3
"""
Compare RSIDIS 2026 and 2025 HV snapshots.

Inputs:
  HV_Check/rsidis_2026_HV_27087.txt
  HV_Check/rsidis_2025_HV_24425.txt

Outputs:
  HV_Check/results/HV_Compare.html
  HV_Check/results/HV_Compare.csv

The detector key is normalized only for the known EPICS naming difference
V0Setr -> V0Set. Voltage values are matched when |2026 - 2025| <= 5.
Mismatched or missing value cells are highlighted with a red border in HTML.
"""

from __future__ import annotations

import csv
import html
import re
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Dict, Iterable, Optional


ROOT = Path(__file__).resolve().parents[1]
FILE_2026 = ROOT / "rsidis_2026_HV_27087.txt"
FILE_2025 = ROOT / "rsidis_2025_HV_24425.txt"
OUT_DIR = ROOT / "results"
OUT_HTML = OUT_DIR / "HV_Compare.html"
OUT_CSV = OUT_DIR / "HV_Compare.csv"
MATCH_TOLERANCE = Decimal("5")


@dataclass(frozen=True)
class HVValue:
    key_raw: str
    value_text: str
    value: Decimal
    comment: str


@dataclass
class Row:
    detector: str
    value_2026: str
    value_2025: str
    diff: str
    status: str
    raw_2026: str
    raw_2025: str
    comment_2026: str
    comment_2025: str


LINE_RE = re.compile(r"^\s*(\S+)\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(?:#\s*(.*))?$")


def normalize_key(key: str) -> str:
    # 2025 snapshots use V0Setr where 2026 uses V0Set for the same setting.
    if key.endswith(":V0Setr"):
        return key[: -len(":V0Setr")] + ":V0Set"
    return key


def parse_file(path: Path) -> Dict[str, HVValue]:
    values: Dict[str, HVValue] = {}
    with path.open() as f:
        for line_number, line in enumerate(f, start=1):
            match = LINE_RE.match(line)
            if not match:
                continue
            key_raw, value_text, comment = match.groups()
            try:
                value = Decimal(value_text)
            except InvalidOperation:
                continue
            key = normalize_key(key_raw)
            if key in values:
                raise ValueError(f"Duplicate normalized key {key!r} in {path} at line {line_number}")
            values[key] = HVValue(
                key_raw=key_raw,
                value_text=value_text,
                value=value,
                comment=comment or "",
            )
    return values


def decimal_to_text(value: Optional[Decimal]) -> str:
    if value is None:
        return "MISSING"
    return format(value, "f")


def make_rows(v2026: Dict[str, HVValue], v2025: Dict[str, HVValue]) -> list[Row]:
    rows: list[Row] = []
    for key in sorted(set(v2026) | set(v2025)):
        a = v2026.get(key)
        b = v2025.get(key)
        if a is None:
            rows.append(
                Row(
                    detector=key,
                    value_2026="MISSING",
                    value_2025=decimal_to_text(b.value if b else None),
                    diff="",
                    status="MISSING_IN_2026",
                    raw_2026="",
                    raw_2025=b.key_raw if b else "",
                    comment_2026="",
                    comment_2025=b.comment if b else "",
                )
            )
            continue
        if b is None:
            rows.append(
                Row(
                    detector=key,
                    value_2026=decimal_to_text(a.value),
                    value_2025="MISSING",
                    diff="",
                    status="MISSING_IN_2025",
                    raw_2026=a.key_raw,
                    raw_2025="",
                    comment_2026=a.comment,
                    comment_2025="",
                )
            )
            continue

        diff = a.value - b.value
        status = "MATCH" if abs(diff) <= MATCH_TOLERANCE else "MISMATCH"
        rows.append(
            Row(
                detector=key,
                value_2026=decimal_to_text(a.value),
                value_2025=decimal_to_text(b.value),
                diff=decimal_to_text(diff),
                status=status,
                raw_2026=a.key_raw,
                raw_2025=b.key_raw,
                comment_2026=a.comment,
                comment_2025=b.comment,
            )
        )
    return rows


def count_status(rows: Iterable[Row]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row.status] = counts.get(row.status, 0) + 1
    return counts


def write_csv(rows: list[Row], path: Path) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Detector", "2026 rsidis", "2025 rsidis", "diff", "status"])
        for row in rows:
            writer.writerow([row.detector, row.value_2026, row.value_2025, row.diff, row.status])


def td(text: str, *, bad: bool = False) -> str:
    cls = ' class="bad"' if bad else ""
    return f"<td{cls}>{html.escape(text)}</td>"


def write_html(rows: list[Row], path: Path) -> None:
    counts = count_status(rows)
    total = len(rows)
    mismatched = total - counts.get("MATCH", 0)

    body_rows = []
    for row in rows:
        bad = row.status != "MATCH"
        title_2026 = f' title="raw: {html.escape(row.raw_2026)}; {html.escape(row.comment_2026)}"' if row.raw_2026 or row.comment_2026 else ""
        title_2025 = f' title="raw: {html.escape(row.raw_2025)}; {html.escape(row.comment_2025)}"' if row.raw_2025 or row.comment_2025 else ""
        cls = "bad" if bad else ""
        body_rows.append(
            "<tr>"
            + td(row.detector)
            + f'<td class="{cls}"{title_2026}>{html.escape(row.value_2026)}</td>'
            + f'<td class="{cls}"{title_2025}>{html.escape(row.value_2025)}</td>'
            + td(row.diff)
            + td(row.status)
            + "</tr>"
        )

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>RSIDIS HV Compare</title>
  <style>
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      margin: 24px;
      color: #111;
    }}
    h1 {{
      margin-bottom: 4px;
    }}
    .meta {{
      margin: 2px 0;
      color: #444;
      font-size: 14px;
    }}
    .summary {{
      margin: 18px 0;
      border-collapse: collapse;
    }}
    .summary td, .summary th {{
      border: 1px solid #ccc;
      padding: 6px 10px;
    }}
    table.data {{
      border-collapse: collapse;
      font-size: 13px;
      width: 100%;
    }}
    table.data th {{
      position: sticky;
      top: 0;
      background: #f3f3f3;
      border: 1px solid #ccc;
      padding: 6px 8px;
      text-align: left;
    }}
    table.data td {{
      border: 1px solid #ddd;
      padding: 5px 8px;
      white-space: nowrap;
    }}
    table.data td.bad {{
      border: 2px solid #c62828;
      background: #fff7f7;
    }}
  </style>
</head>
<body>
  <h1>RSIDIS HV Comparison</h1>
  <div class="meta">2026 sample: {html.escape(FILE_2026.name)}</div>
  <div class="meta">2025 sample: {html.escape(FILE_2025.name)}</div>
  <div class="meta">Matched when |2026 - 2025| &le; {html.escape(format(MATCH_TOLERANCE, "f"))}. Known key normalization: V0Setr &rarr; V0Set.</div>

  <table class="summary">
    <tr><th>Total detectors</th><td>{total}</td></tr>
    <tr><th>Matched</th><td>{counts.get("MATCH", 0)}</td></tr>
    <tr><th>Mismatched or missing</th><td>{mismatched}</td></tr>
    <tr><th>Mismatched values</th><td>{counts.get("MISMATCH", 0)}</td></tr>
    <tr><th>Missing in 2026</th><td>{counts.get("MISSING_IN_2026", 0)}</td></tr>
    <tr><th>Missing in 2025</th><td>{counts.get("MISSING_IN_2025", 0)}</td></tr>
  </table>

  <table class="data">
    <thead>
      <tr>
        <th>Detector</th>
        <th>2026 rsidis</th>
        <th>2025 rsidis</th>
        <th>diff</th>
        <th>status</th>
      </tr>
    </thead>
    <tbody>
      {"".join(body_rows)}
    </tbody>
  </table>
</body>
</html>
"""
    path.write_text(html_text)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    v2026 = parse_file(FILE_2026)
    v2025 = parse_file(FILE_2025)
    rows = make_rows(v2026, v2025)
    write_csv(rows, OUT_CSV)
    write_html(rows, OUT_HTML)

    counts = count_status(rows)
    print(f"Wrote {OUT_HTML}")
    print(f"Wrote {OUT_CSV}")
    print(f"Total detectors: {len(rows)}")
    print(f"Matched: {counts.get('MATCH', 0)}")
    print(f"Mismatched values: {counts.get('MISMATCH', 0)}")
    print(f"Missing in 2026: {counts.get('MISSING_IN_2026', 0)}")
    print(f"Missing in 2025: {counts.get('MISSING_IN_2025', 0)}")


if __name__ == "__main__":
    main()
