#!/usr/bin/env python3
"""
build_sigma_model.py  (robust v5)

Fixes:
  1) Runs calc_semi_xsec with cwd=directory containing the executable, so relative
     files like f1f2tables/... are found regardless of where you launch this script.
  2) Parses ONLY the numeric line immediately after the header:
       "#phi, sighad, sigsemi (ub/GeV2/sr2)"
     This prevents accidentally grabbing numbers from other lines (like input echo
     lines or a crash backtrace).
  3) Checks return code; if calc exits non-zero, we treat it as a failure and show
     a short debug snippet.

Usage:
  python3 tools/build_sigma_model.py --settings_root settings \
    --calc_exe ../../simc_gfortran/util/sidisxsec/calc_semi_xsec --nphibin 1
"""

from __future__ import annotations
import argparse
import json
import re
import subprocess
from pathlib import Path
from typing import Optional, Dict, Any, Tuple, List

PASS_TO_EBEAM = {4: 8.5831, 5: 10.6716}
TARGET_AZ = {"lh2": (1, 1), "ld2": (2, 1)}
RUN_TYPES_SIDIS = {"pi+sidis": +1, "pi-sidis": -1}

# float with optional exponent, supports Fortran D exponent
FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")

def load_manifest(path: Path) -> Dict[str, Any]:
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}

def norm(s: str) -> str:
    return str(s).strip().lower()

def parse_dir_number(token: str, prefix: str) -> Optional[float]:
    t = norm(token)
    pfx = norm(prefix)
    if not t.startswith(pfx):
        return None
    body = t[len(pfx):]
    sign = 1.0
    if body.startswith("neg"):
        sign = -1.0
        body = body[3:]
    body = body.replace("p", ".")
    if body == "":
        return None
    try:
        return sign * float(body)
    except Exception:
        return None

def find_first_token(tokens: List[str], startswith: str) -> Optional[str]:
    sw = norm(startswith)
    for t in tokens:
        if norm(t).startswith(sw):
            return str(t)
    return None

def find_embedded(tokens: List[str], regex: str) -> Optional[str]:
    pat = re.compile(regex, re.IGNORECASE)
    for t in tokens:
        m = pat.search(str(t))
        if m:
            return m.group(1)
    return None

def resolve_setting_inputs(setting_dir: Path) -> Optional[Tuple[float,float,float,float,float,int,int,int]]:
    manifest = load_manifest(setting_dir / "manifest.json")
    tokens: List[str] = []

    if isinstance(manifest.get("path_keys"), list):
        tokens.extend([str(x) for x in manifest["path_keys"]])

    tokens.extend(list(setting_dir.parts))
    tokens_n = [norm(t) for t in tokens]

    run_type = None
    for t in tokens_n:
        if t in RUN_TYPES_SIDIS:
            run_type = t
            break
    if run_type is None:
        return None
    pion_charge = RUN_TYPES_SIDIS[run_type]

    beam_pass = None
    for t in tokens_n:
        m = re.search(r"pass(\d+)", t)
        if m:
            beam_pass = int(m.group(1))
            break
    if beam_pass not in PASS_TO_EBEAM:
        return None
    ebeam = PASS_TO_EBEAM[beam_pass]

    target = None
    for t in tokens_n:
        if t in TARGET_AZ:
            target = t
            break
    if target is None:
        return None
    A, Z = TARGET_AZ[target]

    tok_x  = find_first_token(tokens, "x")
    tok_z  = find_first_token(tokens, "z")
    tok_q2 = find_first_token(tokens, "q2")

    x  = parse_dir_number(tok_x, "x") if tok_x else None
    z  = parse_dir_number(tok_z, "z") if tok_z else None
    Q2 = parse_dir_number(tok_q2, "q2") if tok_q2 else None

    tok_thpq = find_first_token(tokens, "thpq")
    if tok_thpq:
        thpq = parse_dir_number(tok_thpq, "thpq")
    else:
        embedded = find_embedded(tokens, r"(thpq(?:neg)?\d+(?:p\d+)?)")
        thpq = parse_dir_number(embedded, "thpq") if embedded else None

    # skip placeholders like neg999
    if any(v is None for v in (x, z, Q2, thpq)):
        return None
    if any(abs(v) > 100 for v in (x, z, Q2, thpq)):
        return None

    return float(ebeam), float(Q2), float(x), float(z), float(thpq), int(A), int(Z), int(pion_charge)

def parse_sigsemi_from_output_strict(output: str) -> Optional[float]:
    """
    Strict: only accept the numeric line right after the header.
    """
    lines = output.splitlines()
    for i, line in enumerate(lines):
        if "#phi" in line.lower() and "sigsemi" in line.lower():
            if i + 1 >= len(lines):
                return None
            cand = lines[i + 1]
            nums = FLOAT_RE.findall(cand)
            if len(nums) < 3:
                return None
            vals = [float(n.replace("D", "E").replace("d", "E")) for n in nums]
            return vals[2]
    return None

def run_calc(calc_exe: Path, ebeam: float, Q2: float, x: float, z: float, thpq: float,
             A: int, Z: int, pion_charge: int, nphibin: int) -> float:
    inp = f"{ebeam} {Q2} {x} {z} {thpq} {A} {Z} {pion_charge} {nphibin}\n"
    p = subprocess.run([str(calc_exe)],
                       input=inp.encode("utf-8"),
                       stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT,
                       cwd=str(calc_exe.parent),
                       check=False)
    out = p.stdout.decode("utf-8", errors="replace")

    if p.returncode != 0:
        snippet = "\n".join(out.splitlines()[:120])
        raise RuntimeError(f"calc_semi_xsec exit code={p.returncode}\n--- calc output (first 120 lines) ---\n{snippet}")

    sigma = parse_sigsemi_from_output_strict(out)
    if sigma is None:
        snippet = "\n".join(out.splitlines()[:120])
        raise RuntimeError("could not parse sigsemi from calc output (strict)\n--- calc output (first 120 lines) ---\n" + snippet)

    return sigma

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--settings_root", required=True)
    ap.add_argument("--calc_exe", required=True)
    ap.add_argument("--nphibin", type=int, default=1)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--dry_run", action="store_true")
    args = ap.parse_args()

    settings_root = Path(args.settings_root).expanduser().resolve()
    calc_exe = Path(args.calc_exe).expanduser().resolve()

    if not settings_root.is_dir():
        raise SystemExit(f"settings_root not found: {settings_root}")
    if not calc_exe.exists():
        raise SystemExit(f"calc_exe not found: {calc_exe}")

    manifests = sorted(settings_root.rglob("manifest.json"))

    ok = skip = fail = 0
    for mp in manifests:
        setting_dir = mp.parent
        out_txt = setting_dir / "sigma_model.txt"

        inputs = resolve_setting_inputs(setting_dir)
        if inputs is None:
            skip += 1
            continue

        if out_txt.exists() and not args.overwrite:
            skip += 1
            continue

        ebeam, Q2, x, z, thpq, A, Z, pion_charge = inputs

        if args.dry_run:
            print(f"DRY setting={setting_dir} ebeam={ebeam} Q2={Q2} x={x} z={z} thpq={thpq} A={A} Z={Z} pion_charge={pion_charge} nphibin={args.nphibin}")
            ok += 1
            continue

        try:
            sigma = run_calc(calc_exe, ebeam, Q2, x, z, thpq, A, Z, pion_charge, args.nphibin)
            out_txt.write_text(
                f"{sigma:.10g}\n"
                f"# sigma_model (sigsemi) in ub/GeV^2/sr^2\n"
                f"# ebeam={ebeam} Q2={Q2} x={x} z={z} thpq={thpq} A={A} Z={Z} pion_charge={pion_charge} nphibin={args.nphibin}\n"
            )
            ok += 1
        except Exception as e:
            fail += 1
            print(f"FAIL setting={setting_dir} err={e}")

    print(f"DONE ok={ok} skip={skip} fail={fail}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
