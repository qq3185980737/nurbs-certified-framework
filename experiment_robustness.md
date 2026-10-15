# NURBS 统一数值框架 —— 对抗/退化/尺度鲁棒性实验结果

**seed**: 20261012  **M/族**: 400

| family | detail | method | n | familyAcc% | convRate% | avgIter | p50Time(ms) | p99Time(ms) | distP50 | distMax |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SHALLOW | angle=6° | system | 400 | 100.00 | 100.0 | 3.4 | 0.75 | 3.00 | 0.00e+0 | 8.39e-9 |
| SHALLOW | angle=6° | GN-baseline | 400 | 100.00 | 100.0 | 1.4 | 0.19 | 0.54 | 0.00e+0 | 0.00e+0 |
| SHALLOW | angle=10° | system | 400 | 100.00 | 100.0 | 3.4 | 0.74 | 2.01 | 0.00e+0 | 9.07e-9 |
| SHALLOW | angle=10° | GN-baseline | 400 | 100.00 | 100.0 | 1.4 | 0.18 | 0.45 | 0.00e+0 | 0.00e+0 |
| SHALLOW | angle=12° | system | 400 | 100.00 | 100.0 | 3.4 | 0.68 | 1.31 | 0.00e+0 | 9.84e-9 |
| SHALLOW | angle=12° | GN-baseline | 400 | 100.00 | 100.0 | 1.4 | 0.18 | 0.35 | 0.00e+0 | 0.00e+0 |
| SHALLOW | angle=14° | system | 400 | 100.00 | 100.0 | 3.8 | 0.70 | 1.34 | 0.00e+0 | 8.94e-9 |
| SHALLOW | angle=14° | GN-baseline | 400 | 100.00 | 100.0 | 1.4 | 0.18 | 0.37 | 0.00e+0 | 0.00e+0 |
| SHALLOW | angle=18° | system | 400 | 100.00 | 100.0 | 3.4 | 0.70 | 1.35 | 0.00e+0 | 9.27e-9 |
| SHALLOW | angle=18° | GN-baseline | 400 | 100.00 | 100.0 | 1.5 | 0.18 | 0.44 | 0.00e+0 | 0.00e+0 |
| SCALE_1e3 |  | system | 400 | 100.00 | 100.0 | 4.6 | 0.78 | 2.21 | 0.00e+0 | 9.29e-9 |
| SCALE_1e3 |  | GN-baseline | 400 | 100.00 | 100.0 | 1.8 | 0.18 | 0.41 | 0.00e+0 | 0.00e+0 |
| EXTREME_W |  | system | 400 | 100.00 | 100.0 | 3.0 | 0.45 | 1.69 | 1.11e-16 | 9.11e-9 |
| EXTREME_W |  | GN-baseline | 400 | 100.00 | 100.0 | 1.1 | 0.15 | 0.45 | 0.00e+0 | 4.44e-16 |
| DEGENERATE |  | system | 400 | 100.00 | 100.0 | 3.3 | 0.26 | 1.08 | 5.99e-12 | 9.94e-9 |
| DEGENERATE |  | GN-baseline | 400 | 100.00 | 100.0 | 1.0 | 0.10 | 0.33 | 0.00e+0 | 0.00e+0 |
| NEAR_PARALLEL |  | system | 400 | 100.00 | 100.0 | 32.2 | 1.28 | 10.22 | 4.76e-3 | 9.73e-3 |
| NEAR_PARALLEL |  | GN-baseline | 400 | 100.00 | 100.0 | 2.4 | 0.18 | 0.43 | 6.17e-3 | 9.94e-3 |
| COINCIDENT |  | system | 400 | 100.00 | 100.0 | 1.5 | 0.25 | 0.66 | 0.00e+0 | 9.68e-9 |
| COINCIDENT |  | GN-baseline | 400 | 100.00 | 100.0 | 1.0 | 0.14 | 0.36 | 0.00e+0 | 0.00e+0 |

---
*由 experiment_robustness.mjs 生成*
