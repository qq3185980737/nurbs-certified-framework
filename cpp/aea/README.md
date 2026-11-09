# AEA reproduction package

Reproduction code for:

> Verifiable Certificates for Surface–Surface Intersection Output in CAD/CAE Kernels

This directory contains the AEA-specific verifier, data pipeline, instrumented
benchmark and result summaries used in E10/E11. The canonical C++ sources stay in
their existing locations:

- `cpp/include/nurbs/cert_ssi.hpp`
- `cpp/tests/bench_cert_step.cpp`
- `cpp/tests/bench_cert_batch.cpp`
- `cpp/tests/cert_verify.cpp`

The package does **not** redistribute the ABC, Fusion 360 Gallery or MFCAD++ STEP
files.

## Layout

```text
cpp/aea/
  README.md
  build.sh                 build benchmark + C++/GMP verifier
  requirements.txt
  data/
    model_ids.tsv          model IDs / file names used for E10/E11
    download_fusion.py     first N Fusion 360 Gallery STEP files (remote zip)
    download_mfcad.py      first N MFCAD++ STEP files (sparse GitHub clone)
    step_to_case.py        STEP -> rational Bézier patch data
    batch_step_to_case.py  batch wrapper
    run_bench.sh           run bench_cert_step over a case directory
    occt_compare.py        face-level OCCT section comparison
    summarize.py           aggregate RESULT lines
  results/                 E10/E11 summaries and OCCT pair tables
  verify/
    ssi_verifier.py        Python reference verifier (stdlib only)
    certs/                 small example certificates
  data/raw/                local cache: raw STEP files (gitignored)
  data/case/               local cache: converted .txt cases (gitignored)
  certificates/            local cache: full E1 set + generated certs (gitignored)
```

## Local data cache (not committed)

The large inputs and certificate archives are stored under `cpp/aea` for local
reproduction, but are excluded from git:

```text
data/raw/fusion/           100 Fusion STEP files
data/raw/mfcad/            100 MFCAD STEP files
data/raw/mfcad_repo/       sparse MFCAD clone (~190 MB)
data/raw/abc/              ABC sample + metadata
data/raw/samples/          small Fusion/MFCAD samples
data/case/fusion/          100 converted Fusion cases
data/case/mfcad/           100 converted MFCAD cases
certificates/e1/           645 E1 certificates
certificates/samples/      individual cert_*.json files
certificates/generated/    acis + random certificate sets from the generator
```

These paths are in `.gitignore`; they are local caches, not part of the git
history. The full certificate archive is also on Zenodo (`10.5281/zenodo.23607833`).

## Environment

- Python 3.12
- `cadquery-ocp==8.0.1.1.0` (OCP bindings used by `step_to_case.py` and
  `occt_compare.py`)
- `remotezip`, `requests`, `numpy`
- C++17 compiler (g++), GMP for `cert_verify.cpp`
- The benchmark binary is built for Linux/WSL. On Windows, run `build.sh` under WSL.

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/WSL: source .venv/bin/activate
pip install -r requirements.txt
```

## Quick start

Run from `cpp/aea`:

```bash
./build.sh
```

This builds:

- `build/bench_cert_step` from `../tests/bench_cert_step.cpp`
- `build/cert_verify` from `../tests/cert_verify.cpp`

### 1. Download a 100-model sample

Fusion 360 Gallery segmentation:

```bash
python data/download_fusion.py --n 100 --out data/raw/fusion
```

MFCAD++ (sparse clone of the public mirror):

```bash
python data/download_mfcad.py --n 100 --out data/raw/mfcad
```

The selection is the first 100 filenames in lexicographic order, matching the paper.

### 2. Convert STEP to benchmark data

```bash
python data/batch_step_to_case.py \
  --input-dir data/raw/fusion --pattern '*.stp' \
  --out-dir data/case/fusion
python data/batch_step_to_case.py \
  --input-dir data/raw/mfcad --pattern '*.step' \
  --out-dir data/case/mfcad
```

### 3. Run the benchmark

```bash
bash data/run_bench.sh data/case/fusion
bash data/run_bench.sh data/case/mfcad
python data/summarize.py --input results/bench_summary.tsv
```

### 4. OCCT face-level comparison (optional)

```bash
python data/occt_compare.py \
  --input-dir data/raw/fusion --pattern '*.stp' \
  --out results/occt_fusion.tsv
```

### 5. Verify a certificate

```bash
python verify/ssi_verifier.py verify/certs/cert_demo.json
./build/cert_verify verify/certs/cert_demo.json
```

## Expected E11 results

For the first 100 Fusion 360 and first 100 MFCAD++ filenames (τ = 0.02):

| dataset | patch pairs | disjoint | 1D certified | coplanar | tangent | cover-limit | patch-boundary | undecided |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Fusion 100 | 1,855,412 | 1,850,650 | 419 | 4,056 | 76 | 1 | 128 | 82 |
| MFCAD 100 | 7,954 | 7,208 | 4 | 742 | 0 | 0 | 0 | 0 |

The Fusion/MFCAD samples are single-solid parts; non-adjacent pairs test false
positives and adjacency handling. ABC is the multi-body intersecting case in E10.

## E10 ABC model

The ABC model used in the paper is:

```text
00000050_80d90bfdd2e74e709956122a
```

Its STEP file is not included; use the ABC dataset page to obtain it, then:

```bash
python data/step_to_case.py /path/to/model.step data/case/abc/model.txt
bash data/run_bench.sh data/case/abc
```

Expected (after excluding 98 adjacent face pairs, τ = 0.02): 2,750 patch pairs,
2,610 disjoint, 79 1D certified, 43 coplanar, 3 cover-limit, 14 patch-boundary,
1 undecided, about 39.4 s.

## Data and licence notes

- ABC models: copyright belongs to their creators; research use under the ABC
  dataset terms; do not redistribute the STEP files.
- Fusion 360 Gallery: research dataset; see the Autodesk license.
- MFCAD++: CC BY; see the Queen's University Belfast dataset page.
- This directory ships **scripts, model IDs and results only**, not the raw STEP data.
- Full certificate archive and both verifiers: Zenodo `10.5281/zenodo.23607833`.

## Scope

The certificate contract covers untrimmed tensor-product patches. The B-rep pipeline
converts each face to rational Bézier patches and excludes face pairs that share a
B-rep edge; trimming, p-curves and B-rep topology are outside the current certificate
model. E10/E11 are classification/coverage runs, not full certificate serialization
for every patch pair.
