# Optional C++ verifier

`cert_verify.cpp` is a standalone C++17 port of the reference verifier. It is provided
for cross-checking and batch runs; the Python package remains the reference
implementation and the one exercised by `pytest`.

## Build

```bash
sudo apt install libgmp-dev
g++ -std=c++17 -O2 cert_verify.cpp -lgmpxx -lgmp -o cert_verify
```

## Run

```bash
./cert_verify ../examples/certs/cert_demo.json
./cert_verify ../examples/certs/cert_chunked.json --chunked --jobs 1
./cert_verify --batch ../examples/certs
```

Expected for the demo certificate:

```text
FINAL: ACCEPT
```

The C++ verifier uses the same certificate contract as the Python verifier. See
`../docs/certificate-schema.md`.
