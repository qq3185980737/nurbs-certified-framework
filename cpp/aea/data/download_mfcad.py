#!/usr/bin/env python3
"""Fetch the first N MFCAD++ STEP files from the public GitHub mirror.

The mirror is sparse-cloned (blob filter + sparse checkout) so the full history is
not downloaded. Pass --repo to reuse an existing clone.
"""
import argparse
import shutil
import subprocess
from pathlib import Path

URL = "https://github.com/hducg/MFCAD.git"


def ensure_repo(repo: Path):
    if (repo / ".git").exists():
        return
    repo.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["git", "clone", "--depth", "1", "--filter=blob:none",
                    "--no-checkout", URL, str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "sparse-checkout", "init", "--cone"], check=True)
    subprocess.run(["git", "-C", str(repo), "sparse-checkout", "set", "dataset/step"], check=True)
    subprocess.run(["git", "-C", str(repo), "checkout"], check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--repo", type=Path, default=Path("data/raw/mfcad_repo"))
    ap.add_argument("--out", type=Path, default=Path("data/raw/mfcad"))
    args = ap.parse_args()
    ensure_repo(args.repo)
    src = args.repo / "dataset" / "step"
    files = sorted(src.glob("*.step"))[: args.n]
    if not files:
        raise SystemExit(f"no .step files under {src}")
    args.out.mkdir(parents=True, exist_ok=True)
    for f in files:
        dest = args.out / f.name
        if not dest.exists() or dest.stat().st_size == 0:
            shutil.copy2(f, dest)
            print("copied", dest.name, dest.stat().st_size)
    print(f"selected {len(files)} STEP files in {args.out}")


if __name__ == "__main__":
    main()
