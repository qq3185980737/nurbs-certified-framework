# Install

## Requirements

- Python 3.9 or newer.
- No third-party Python runtime dependencies.
- Optional: `pytest` for the test suite.
- Optional: a C++17 compiler and GMP for the standalone C++ verifier.

## From the repository

```bash
git clone https://github.com/qq3185980737/nurbs-certified-framework.git
cd nurbs-certified-framework/softwarex
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
python -m pip install -e .
```

The `ssi-cert` command should now be on `PATH`:

```bash
ssi-cert --version
ssi-cert verify examples/certs/cert_demo.json
```

## Without installation

```bash
PYTHONPATH=src python -m ssi_cert verify examples/certs/cert_demo.json
```

On PowerShell:

```powershell
$env:PYTHONPATH = "src"
python -m ssi_cert verify examples/certs/cert_demo.json
```

## C++ verifier (optional)

```bash
sudo apt install libgmp-dev
g++ -std=c++17 -O2 cpp/cert_verify.cpp -lgmp -o cert_verify
./cert_verify examples/certs/cert_demo.json
```

The Python verifier is the reference implementation. The C++ verifier is provided for
cross-checking and batch use; it is not required by the Python package.
