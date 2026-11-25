# Reproducing the benchmark

The small certificates in `examples/certs/` are enough to check installation and the
accept/reject contract. They are not the full benchmark.

The full corpus, including the 645-certificate E1 set, the E12 challenge cases, and the
E10/E11 summaries, is archived on Zenodo:

```text
https://doi.org/10.5281/zenodo.24126465
```

Download `aea-artifact-v1.6.0.zip`, verify its `CHECKSUMS.sha256`, and use the scripts in
its `data/` directory to reproduce E10/E11. The raw STEP files are not redistributed;
`dataset_versions.tsv` pins the dataset versions and selection rule.

## Quick batch check

```bash
ssi-cert verify-dir examples/certs
```

Expected summary:

```text
SUMMARY files=... accept=... reject=0 error=0
```

For a chunked certificate:

```bash
ssi-cert verify examples/certs/cert_chunked.json --chunked --jobs 1
```

`--jobs` greater than 1 uses process-based parallelism and is supported on POSIX systems.
On Windows, use `--jobs 1`.
