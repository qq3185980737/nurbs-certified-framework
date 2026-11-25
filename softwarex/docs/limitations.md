# Scope and limitations

SSI-Cert is a certificate checker, not a geometry kernel and not a complete SSI solver.

## Covered

- Rational Bézier patches with positive weights.
- Untrimmed tensor-product parameter domains.
- Exact rational tree/leaf exclusion.
- Krawczyk existence boxes.
- Coverage, miss, and branch-coverage claims.
- Chunked certificates.

## Not covered

- Trimmed B-rep faces and p-curves.
- Negative or mixed-sign weights.
- Self-intersection certification.
- Approximate or floating-point-only certificates.
- Kernel correctness or algorithm completeness.
- Production-scale performance.

## Meaning of REJECT

`REJECT` means the supplied certificate does not prove the claimed property under this
contract. It does not prove the geometric statement is false. A different certificate
may still verify.

## Data

The full benchmark corpus is archived on Zenodo. Raw ABC, Fusion 360 Gallery, and
MFCAD++ STEP files are not redistributed because of their licences; the artifact
contains model IDs, conversion scripts, and checksums.
