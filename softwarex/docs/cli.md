# Command-line interface

```text
ssi-cert verify <certificate.json> [--chunked] [--jobs N] [--quiet]
ssi-cert verify-dir <directory> [--chunked] [--jobs N] [--quiet]
ssi-cert info <certificate.json>
ssi-cert --version
```

Exit codes:

```text
0  certificate accepted
1  certificate rejected
2  malformed input, missing file, or verification error
```

Examples:

```bash
ssi-cert verify examples/certs/cert_demo.json
ssi-cert verify examples/certs/cert_chunked.json --chunked --jobs 2
ssi-cert verify-dir examples/certs
ssi-cert info examples/certs/cert_demo.json
```

`info` reports the certificate SHA-256, size, top-level keys, number of patches, and
whether the certificate is chunked. It does not verify the certificate.

`verify-dir` checks every `*.json` file in the directory and prints one line per file
plus a summary. It returns `0` only when every certificate is accepted.
