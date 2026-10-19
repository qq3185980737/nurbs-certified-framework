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

Interactive Web demo — **must be served over HTTP; do not open `index.html` directly as `file://`.** The math core is WebAssembly (`cpp/wasm/nurbs.wasm`) loaded via `fetch`, which browsers block under the `file://` protocol — the canvas and control-point inputs stay blank and you cannot type into them. Start a local server and open `index.html`:

```sh
# from the repository root
python -m http.server 8000
# then open in your browser
http://localhost:8000/index.html
```

On Windows, double-click `启动服务器.bat` to start the server and open the browser automatically (a browser window will open; keep the console window open while using the app). WASM rebuild / front-end build: see `cpp/wasm/README.md` and `build.mjs`.

### Web demo features

The demo exposes the certified machinery directly, matching the manuscript's claims:

- **All intersections, not just the closest** — the result panel lists every root inside the tolerance band (`u`, `v`, residual `|P1(u)−P2(v)|²`, point coordinates, iteration count).
- **Condition number κ** — the root Jacobian condition number is shown; tangency / singularity reports `+∞` (the paper's κ > 10⁴ transverse-vs-tangent classifier).
- **Certified proof summary** — a Krawczyk interval enumeration reports how many transverse roots are *proved* (`proven_count`), how many are suspected tangent / degenerate (`tangent_count`), how many boxes were evaluated, and whether the full parameter domain was covered (`complete`). Exposed in the C API as `nurbs_enumerate_intersections` and in the WASM frontend as `enumerateIntersections`.
- Eight preset scenarios (transverse, multiple, tangent, higher-order tangent, near/far separated, endpoint, circle-vs-curve) exercise all of these outputs.
- **Built-in stress test** — the "运行压力测试" button runs 100 random curve pairs with a **fixed seed** (fully reproducible) and reports success rate, average iterations, condition-number finite/singular counts, and a Krawczyk proof summary (proved / suspected-tangent / full-domain coverage) over a sampled subset.

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
- Frontend pressure regression: serve `stress_test.html` over HTTP (60 curated combos + 300 random pairs; captures any alert/exception and reports κ/proof stats).
- Ground-truth verification: `_ground_truth.py` (e.g. lemniscate self-intersection coordinates to ~1e-16, `min‖∂S/∂s×∂S/∂v‖ = √2·π` independently verified).

## License

MIT — see [LICENSE](LICENSE). Copyright (c) 2026 Yang Liu, Xun Wang.

## Related work

This repository is the open-source companion of a research manuscript on certified unified interrogation of NURBS/Bézier geometry (Computer-Aided Design, under review).
