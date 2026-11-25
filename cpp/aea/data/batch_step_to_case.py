#!/usr/bin/env python3
"""Convert every STEP/STP file in a directory to bench_cert_step text format."""
import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--pattern", default="*.stp")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(args.input_dir.glob(args.pattern))
    if not files:
        raise SystemExit(f"no files matching {args.pattern} in {args.input_dir}")
    for f in files:
        out = args.out_dir / (f.stem + ".txt")
        if out.exists() and out.stat().st_size > 0:
            continue
        print("convert", f.name)
        subprocess.run([sys.executable, str(HERE / "step_to_case.py"), str(f), str(out)], check=True)


if __name__ == "__main__":
    main()
