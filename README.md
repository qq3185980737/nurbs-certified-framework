# NURBS Unified Numerical Framework

A unified **certified** numerical framework for NURBS / Bézier geometric interrogation: curve **intersection**, **shortest-distance**, **degeneracy recovery** and **surface self-intersection / patch overlap**, all built on strict interval arithmetic with conservative outward widening and the classical interval **Krawczyk** theorem.

## Capabilities

- **Curve intersection & shortest distance** collapse into a single square residual system whose every solution is certified by the interval Krawczyk test; a global enumeration with adaptive multi-start captures all intersections on every test instance, and a 14-DOP bounding-volume hierarchy accelerates pre-filtering.
- **Degenerate geometry** (point / collinear / rational-weight degeneracies) is classified by sound deterministic criteria, reduced to its embedded dimension, and recovered with certified point-inclusion evidence.
- **Genuine surface self-intersection** is addressed with lemniscate-cylinder patches that truly self-intersect in their image, reduced to a square 2×2 residual system and proved by the classical Krawczyk theorem; regular patches are certified fold-free by interval pruning.
- **Coplanar patch overlap** is detected by convex-hull AABB intersection and supported by point-inclusion evidence of a common region.

All regressions report **certified coverage of every detected root with zero false positives and zero false certificates**.

## Implementation

| Component | Detail |
|---|---|
| C++17 reference | MSVC, `Release /O2 /fp:precise` (no fast-math) — bit-identical to the JS reference |
| JavaScript | Single-source with a WASM build path (`build.mjs`) |
| GPU | OpenCL kernels for large-scale throughput |
| Bindings | Python `ctypes` bindings to the C API/DLL |

## Directory layout

```
src/            JavaScript/Web source (NURBS core, solver, rendering)
cpp/            C++ port, C API, WASM/OpenCL, benchmarks & tests
index.html      Interactive demo
*.mjs / *.py    Regression, benchmark and ground-truth scripts
experiment_*.md / *.json   Detailed experimental results
```

## Experimental coverage

- 9,000 curve pairs (9 geometry types × 1,000) → 4,675 certified intersections, 100% certified, 0 false positives
- 4,000 degenerate-curve + 3,000 degenerate-surface cases → full recovery, 0 false certificates
- 500 self-intersecting + 500 regular patches → 500/500 certified, 0 false positives
- 500 overlapping + 500 separated patch pairs → 500/500 detected, 0 false certificates
- Ablation (BVH 35.7×, adaptive multi-start −55.6%, 14-DOP −52.4% leaf pairs, batch parallel 1.81×/3.65×/4.98×) and threshold sensitivity are reported in the experiment files.

## Build & run

See `README_LOCAL_RUN.txt` and `build.mjs` for build steps and the interactive demo. Benchmarks live under `cpp/tests` / `cpp/build`.

## License

MIT — see [LICENSE](LICENSE). Copyright (c) 2026 Yang Liu, Xun Wang.
