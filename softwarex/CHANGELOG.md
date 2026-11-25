# Changelog

## 0.1.1 — 2026-11-25

- C++ verifier now validates every tree box against `[0,1]^4`, checks axis/children structure, and rejects malformed nodes instead of skipping them.
- Fixed duplicate `cert_ok` accumulation in chunked parallel statistics.
- Unified numeric-string parsing: certificates accept `n/d` or hexadecimal floats only; decimal strings are rejected consistently.
- Corrected C++ build instructions to link `-lgmpxx -lgmp` and to use `-O3 -march=native`.

## 0.1.0 — 2026-11-25
