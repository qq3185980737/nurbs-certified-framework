# NURBS Unified Numerical Framework

A unified **certified** numerical framework for NURBS / Bézier geometric interrogation. Curve **intersection**, **shortest distance**, **degeneracy recovery**, **genuine surface self-intersection**, and **coplanar patch overlap** are all solved from a single mathematical core built on strict interval arithmetic with conservative outward widening and the classical interval **Krawczyk** theorem.

Every reported root is certified — the framework emits **zero false positives and zero false certificates** on its full regression suite.

## Highlights

- **One residual system, two outcomes.** Curve intersection and shortest distance collapse into a single square residual system whose solutions are certified by the interval Krawczyk test. A global enumeration with adaptive multi-start provably captures all intersections; a 14-DOP bounding-volume hierarchy accelerates pre-filtering.
- **Degenerate geometry handled explicitly.** Point / collinear / rational-weight degeneracies are classified by sound deterministic criteria, reduced to their embedded dimension, and recovered with certified point-inclusion evidence.
- **Genuine surface self-intersection.** Lemniscate-cylinder patches that truly self-intersect in their image reduce to a square 2×2 residual system proved by the classical Krawczyk theorem; regular patches are certified fold-free by interval pruning.
- **Coplanar patch overlap.** Convex-hull AABB intersection plus certified point-inclusion evidence of a common region.
- **Identical numerics everywhere.** The C++17 reference is bit-identical to the JavaScript reference (`/fp:precise`, no fast-math).

## Experimental coverage

| Regression | Result |
|---|---|
| 9,000 curve pairs (9 geometry types × 1,000) | 4,675 certified intersections, **100% certified, 0 false positives** |
| 4,000 degenerate curves + 3,000 degenerate surfaces | full recovery, **0 false certificates** |
| 500 self-intersecting + 500 regular patches | 500/500 certified, **0 false positives** |
| 500 overlapping + 500 separated patch pairs | 500/500 detected, **0 false certificates** |

Ablation (BVH 35.7×, adaptive multi-start −55.6% iterations, 14-DOP −52.4% leaf pairs, batch parallelism 1.81× / 3.65× / 4.98×) and threshold-sensitivity results are detailed in `experiment_*.md`.

## Quick start

Interactive Web demo — open `index.html` in a browser (WASM build optional, see `build.mjs`).

C++ benchmark / test binaries live under `cpp/build` (MSVC `Release /O2 /fp:precise`):

```sh
cd cpp
build/test_equivalence.exe      # JS/C++ bit-identical equivalence (eval 302160, solve 4208, find 68 — all pass)
build/bench_cs.exe              # certified curve–surface intersection
build/bench_true_selfintersect_r4.exe  # genuine self-intersecting patches
```

Python `ctypes` bindings to the C API are under `cpp/python`. WASM frontend test: `node cpp/wasm/wasm_test.mjs`.

## Directory layout

```
src/            JavaScript/Web source (NURBS core, solver, rendering)
cpp/            C++17 port, C API, WASM/OpenCL, benchmarks & tests
index.html      Interactive demo
*.mjs / *.py    Regression, benchmark and ground-truth scripts
experiment_*.md / *.json   Detailed experimental results
```

## Reproduce

- Build & run instructions: `README_LOCAL_RUN.txt`, `cpp/CMakeLists.txt`, `build.mjs`.
- Regression / ablation scripts: `*.mjs`, `*.py` at repository root.
- Ground-truth verification: `_ground_truth.py` (e.g. lemniscate self-intersection coordinates to ~1e-16, `min‖∂S/∂s×∂S/∂v‖ = √2·π` independently verified).

## License

MIT — see [LICENSE](LICENSE). Copyright (c) 2026 Yang Liu, Xun Wang.

## Related work

This repository is the open-source companion of a research manuscript on certified unified interrogation of NURBS/Bézier geometry (Computer-Aided Design, under review).
