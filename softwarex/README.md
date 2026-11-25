# SSI-Cert

Independent, exact-rational verification of surface–surface intersection (SSI) output
certificates.

SSI-Cert is the reference software accompanying:

> **Verifiable Certificates for Surface–Surface Intersection Output in CAD Kernels**

It replays a certificate produced by a CAD kernel or SSI solver and decides, without
trusting the generator, whether the certificate proves the claimed coverage, exclusion,
existence, miss, or branch-coverage statement. All geometric predicates are evaluated
with Python `Fraction`; the verifier does not use floating-point tolerances.

## Why

Kernel-internal checks and sampled comparisons can show that an SSI result looks
plausible, but they do not produce evidence that a third party can replay. SSI-Cert
turns post-hoc evidence into a deterministic predicate:

```text
certificate + verifier  ->  ACCEPT / REJECT
```

The certificate is a JSON document. The verifier is a small Python package with no
runtime dependencies outside the standard library.

## Install

```bash
git clone https://github.com/qq3185980737/nurbs-certified-framework.git
cd nurbs-certified-framework/softwarex
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
python -m pip install -e .
```

## Quick start

```bash
ssi-cert verify examples/certs/cert_demo.json
```

Expected output ends with:

```text
FINAL: ACCEPT
ACCEPT examples/certs/cert_demo.json
```

Exit codes:

```text
0  certificate accepted
1  certificate rejected
2  malformed input or verification error
```

Other commands:

```bash
ssi-cert info examples/certs/cert_demo.json
ssi-cert verify-dir examples/certs
ssi-cert verify examples/certs/cert_chunked.json --chunked --jobs 2
```

## What is verified

The verifier implements the contract described in the accompanying manuscript:

- tree structure and dyadic subdivision;
- exact AABB exclusion of leaf boxes;
- exact Krawczyk existence checks for certified boxes;
- miss and coverage claims;
- branch-coverage claims;
- chunked certificates.

The certificate schema is documented in `docs/certificate-schema.md`.

## Scope and limitations

- The verifier checks a certificate, not the kernel that produced it.
- The current contract covers untrimmed tensor-product Bézier/NURBS patches.
- Trimmed B-rep faces, negative weights, and mixed-weight rational patches are outside
  the current contract.
- The verifier is intentionally conservative: `REJECT` means "not proved by this
  certificate", not "the geometric statement is false".
- Large certificates are verified in Python; chunking is available, but this is a
  reference implementation, not a high-performance production checker.

## Tests

```bash
python -m pip install pytest
pytest
```

The suite covers acceptance, rejection, malformed input, CLI exit codes, and the
chunked path.

## Repository layout

```text
softwarex/
  src/ssi_cert/        Python package and CLI
  examples/certs/      small certificates for a first run
  tests/               pytest suite
  docs/                schema and usage notes
  cpp/                 optional C++/GMP verifier source
```

The canonical C++ sources remain in the parent repository:

```text
cpp/include/nurbs/cert_ssi.hpp
cpp/tests/cert_verify.cpp
cpp/tests/bench_cert_batch.cpp
```

## Data and artifact

The full certificate corpus and benchmark artifact are archived on Zenodo:

```text
https://doi.org/10.5281/zenodo.24126465
```

Raw ABC, Fusion 360 Gallery, and MFCAD++ STEP files are not redistributed because of
their licences. The artifact contains model IDs, conversion scripts, checksums, and the
derived certificate data.

## License

MIT. See `LICENSE`. Third-party dataset terms are summarized in `LICENSE_NOTES.md`.

## Citation

If you use this software, cite the accompanying manuscript and the archived artifact.
See `CITATION.cff`.
