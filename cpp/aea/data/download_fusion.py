#!/usr/bin/env python3
"""Download the first N STEP files from the Fusion 360 Gallery segmentation archive.

The archive is read remotely; the 483 MB zip is not downloaded in full.
"""
import argparse
from pathlib import Path

from remotezip import RemoteZip

DEFAULT_URL = (
    "https://fusion-360-gallery-dataset.s3.us-west-2.amazonaws.com/"
    "segmentation/s2.0.1/s2.0.1_extended_step.zip"
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--out", type=Path, default=Path("data/raw/fusion"))
    ap.add_argument("--url", default=DEFAULT_URL)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    z = RemoteZip(args.url)
    names = sorted(n for n in z.namelist() if n.lower().endswith(".stp"))
    if not names:
        raise SystemExit("no .stp files in archive")
    chosen = names[: args.n]
    for name in chosen:
        dest = args.out / Path(name).name
        if dest.exists() and dest.stat().st_size > 0:
            continue
        dest.write_bytes(z.read(name))
        print("downloaded", dest.name, dest.stat().st_size)
    print(f"selected {len(chosen)} / {len(names)} STEP files in {args.out}")


if __name__ == "__main__":
    main()
