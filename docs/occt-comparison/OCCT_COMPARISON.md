# OCCT Kernel Comparison (reproducible)

This directory reproduces the OCCT kernel comparison reported in Section 7.8 of the manuscript
*"Certified Unified Interrogation of NURBS/Bézier Geometry"*. It compares the open-source
**OCCT (Open CASCADE Technology)** kernel against the certified framework on the same fixtures,
using the OCCT standard interrogation APIs.

## Environment

- **WSL2 / Ubuntu 22.04 LTS**, `g++ 11.4`.
- **miniconda3** + conda-forge channel: `conda create -n occ -c conda-forge python=3.11 pythonocc-core`
  (pythonocc-core bundles a modern OCCT 7.x, MSVC/GCC native, no source build needed).
- On Windows the certified framework is already built (MSVC `Release /O2 /fp:precise`). The two
  kernels are run as **independent processes** and compared on output files, so no cross-ABI linking
  is required.

## Scripts

Run each with the `occ` conda env's python (`/root/miniconda3/envs/occ/bin/python`):

| Script | Scope |
|---|---|
| `compare_curve.py` | Curve-curve: transverse / double-transverse / endpoint-tangent / interior-tangent / 2nd-order osculating, via `GeomAPI_ExtremaCurveCurve`. |
| `compare_curve2.py` | Curve-curve: degree-5 Bezier with **5 transverse roots** vs y=0, via `GeomAPI_ExtremaCurveCurve`. |
| `compare_section.py` | Curve-curve via OCCT's standard `BRepAlgoAPI_Section` (multi-root / single / interior-tangent). |
| `compare_ssi.py` | Surface-surface (graph patches): fixtures A `s^2-t^2=0`, B `f=0.3`, C `f=s+0.5t`, via `GeomAPI_IntSS`. |
| `compare_ssi2.py` | Surface-surface: degenerate tangency `f=(s-t)^2` (true intersection = whole diagonal s=t), via `GeomAPI_IntSS`. |
| `compare_ssi_detail.py` | Samples OCCT's single intersection line on fixture A to confirm which diagonal it returns. |

## Observed results (tolerance 1e-9 unless noted)

### Curve-curve (`BRepAlgoAPI_Section`)

| Fixture | True roots | OCCT | Certified framework (Section 7.1) |
|---|---|---|---|
| single transverse | 1 | 1 | certified root |
| double transverse | 2 | 2 | certified roots |
| endpoint tangent | 1 | 1 | certified root |
| interior tangent | 1 | 1 | certified root |
| 2nd-order osculating | 1 | 1 | certified root |
| **degree-5, five transverse roots** | **5** | **1 (misses 4)** | all five certified |

### Surface-surface (`GeomAPI_IntSS`)

| Fixture | True intersection | OCCT | Certified framework (Section 6.3) |
|---|---|---|---|
| A: `f=s^2-t^2` | single diagonal s=t | 1 line (correct) | cover band encloses it |
| B: `f=0.3` | none | 0 lines (correct) | whole domain rejected |
| C: `f=s+0.5t` | tangential corner contact at (0,0) | **0 lines (misses the corner)** | corner cover box retained |
| D: `f=(s-t)^2` | whole diagonal s=t (degenerate tangency) | **0 lines (misses the whole line)** | cover band encloses it |

### Curve-surface (`GeomAPI_IntCS` and `BRepAlgoAPI_Section`)

| Fixture | True roots | OCCT IntCS | OCCT Section | Certified framework (Section 3.4) |
|---|---|---|---|---|
| transverse crossing of a planar graph patch | 1 | 1 | — | certified root |
| curve osculating a plane (tangent point) | 1 | 2 (duplicate point pair; 1 unique) | — | certified root |
| **degree-3 curve, three transverse crossings of a paraboloidal graph patch** (residuals ~1e-16 at each root) | **3** | **0 (misses all 3)** | **0 (misses all 3)** | all three certified |

### Coplanar overlap and surface self-intersection

| Fixture | True result | OCCT | Certified framework |
|---|---|---|---|
| two coplanar graph patches, overlap `[0.25,0.75]^2` (Section 6) | overlap region | `GeomAPI_IntSS` returns 0 lines; `BRepAlgoAPI_Section` returns only boundary edges, no face | Section 6 returns the shared region |
| patch `S(u,v)=(u^2-v^2,2uv,u^2+v^2)` folding along `u=0` (transverse self-intersection, plus isolated parameter degeneration at origin; Section 5) | self-intersection line `x=-z, y=0` | `GeomAPI_IntSS(S,S)` returns 0 lines; `BRepCheck_Analyzer` only reports the face invalid without localizing the intersection | Section 5 detects the self-intersection and certifies the parameter pair |

These are precisely the high-order-contact and small-loop configurations that [17] report commercial
kernels ACIS and OCCT handle unsatisfactorily. The certified interval treatment of Sections 3.4/6.3
keeps provable completeness (no-missing cover, Proposition 1) where the kernel omits a root or a
degenerate contact.

## Notes

- Scripts use only OCCT public APIs (`GeomAPI_IntSS`, `GeomAPI_ExtremaCurveCurve`,
  `BRepAlgoAPI_Section`, `Geom_BezierCurve/Surface`); the graph-patch surfaces are exact degree-2x2
  Bézier representations of `(s,t,f(s,t))`.
- Section 7.1 reports 4,675 certified curve roots across 9,000 pairs (100% certified, 0 false
  positives) for the certified side of the curve comparison.
