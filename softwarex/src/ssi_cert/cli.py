"""Command-line interface for SSI-Cert."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import sys
from pathlib import Path

from . import __version__, verifier


def _load(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _run_one(path: Path, chunked: bool, jobs: int, quiet: bool) -> int:
    try:
        cert = _load(path)
    except Exception as exc:
        print(f"ERROR: {path}: {exc}", file=sys.stderr)
        return 2
    use_chunked = chunked or "chunks" in cert
    sink = io.StringIO() if quiet else None
    try:
        if sink is not None:
            with contextlib.redirect_stdout(sink):
                result = verifier.verify_chunked(cert, jobs=jobs, verbose=True) if use_chunked else verifier.verify(cert, verbose=True)
        else:
            result = verifier.verify_chunked(cert, jobs=jobs, verbose=True) if use_chunked else verifier.verify(cert, verbose=True)
        ok = result[0] if isinstance(result, tuple) else result
    except Exception as exc:  # malformed input is a verification error, not a rejection
        print(f"ERROR: {path}: {exc}", file=sys.stderr)
        return 2
    print(f"{'ACCEPT' if ok else 'REJECT'} {path}")
    return 0 if ok else 1


def cmd_verify(args: argparse.Namespace) -> int:
    return _run_one(Path(args.certificate), args.chunked, args.jobs, args.quiet)


def cmd_verify_dir(args: argparse.Namespace) -> int:
    files = sorted(Path(args.directory).glob("*.json"))
    if not files:
        print(f"ERROR: no .json certificates in {args.directory}", file=sys.stderr)
        return 2
    counts = {0: 0, 1: 0, 2: 0}
    for path in files:
        rc = _run_one(path, args.chunked, args.jobs, args.quiet)
        counts[rc] = counts.get(rc, 0) + 1
    print(f"SUMMARY files={len(files)} accept={counts.get(0, 0)} reject={counts.get(1, 0)} error={counts.get(2, 0)}")
    if counts.get(2, 0):
        return 2
    return 0 if counts.get(1, 0) == 0 else 1


def cmd_info(args: argparse.Namespace) -> int:
    path = Path(args.certificate)
    try:
        raw = path.read_bytes()
        cert = json.loads(raw.decode("utf-8"))
    except Exception as exc:
        print(f"ERROR: {path}: {exc}", file=sys.stderr)
        return 2
    chunks = cert.get("chunks")
    print(f"file: {path}")
    print(f"sha256: {hashlib.sha256(raw).hexdigest()}")
    print(f"size_bytes: {len(raw)}")
    print(f"keys: {', '.join(sorted(cert))}")
    if isinstance(chunks, list):
        print(f"chunks: {len(chunks)}")
    else:
        tree = cert.get("tree") or {}
        print(f"tree_root: {'yes' if tree else 'no'}")
    patches = cert.get("patches")
    if isinstance(patches, list):
        print(f"patches: {len(patches)}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ssi-cert",
        description="Independent exact-rational verifier for SSI output certificates.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_verify = sub.add_parser("verify", help="verify one certificate")
    p_verify.add_argument("certificate", help="path to certificate JSON")
    p_verify.add_argument("--chunked", action="store_true", help="force chunked verification")
    p_verify.add_argument("--jobs", type=int, default=1, help="workers for chunked verification")
    p_verify.add_argument("--quiet", action="store_true", help="suppress verifier diagnostics")
    p_verify.set_defaults(func=cmd_verify)

    p_dir = sub.add_parser("verify-dir", help="verify all *.json files in a directory")
    p_dir.add_argument("directory", help="directory containing certificates")
    p_dir.add_argument("--chunked", action="store_true", help="force chunked verification")
    p_dir.add_argument("--jobs", type=int, default=1, help="workers for chunked verification")
    p_dir.add_argument("--quiet", action="store_true", help="suppress verifier diagnostics")
    p_dir.set_defaults(func=cmd_verify_dir)

    p_info = sub.add_parser("info", help="show certificate metadata without verifying")
    p_info.add_argument("certificate", help="path to certificate JSON")
    p_info.set_defaults(func=cmd_info)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
