# NURBS 统一数值框架 —— 论文实验结果（CAD 方向）

**seed**: 20261012   **运行**: 2026-10-12T17:55:16.133Z

## 表1：基线对照与消融

| family | method | n | familyAcc% | convRate% | avgIter | p50Time(ms) | p99Time(ms) | distP50 | distMax |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| INTERSECT | system(prefilter+LM) | 2000 | 100.00 | 100.0 | 2.6 | 0.63 | 1.55 | 0.00e+0 | 9.55e-9 |
| INTERSECT | LM(no prefilter) | 2000 | 100.00 | 100.0 | 2.5 | 0.59 | 1.13 | 0.00e+0 | 9.87e-9 |
| INTERSECT | GN(no prefilter, baseline) | 2000 | 100.00 | 100.0 | 3.0 | 0.16 | 0.42 | 0.00e+0 | 4.44e-16 |
| TANGENT | system(prefilter+LM) | 400 | 100.00 | 100.0 | 4.2 | 0.60 | 1.91 | 0.00e+0 | 9.95e-9 |
| TANGENT | LM(no prefilter) | 400 | 100.00 | 100.0 | 3.8 | 0.58 | 2.03 | 0.00e+0 | 9.83e-9 |
| TANGENT | GN(no prefilter, baseline) | 400 | 100.00 | 100.0 | 3.7 | 0.18 | 0.41 | 0.00e+0 | 0.00e+0 |
| HIGH_ORDER | system(prefilter+LM) | 400 | 100.00 | 100.0 | 2.3 | 0.88 | 1.63 | 0.00e+0 | 6.95e-9 |
| HIGH_ORDER | LM(no prefilter) | 400 | 100.00 | 100.0 | 2.3 | 0.82 | 1.37 | 0.00e+0 | 4.57e-9 |
| HIGH_ORDER | GN(no prefilter, baseline) | 400 | 100.00 | 100.0 | 3.0 | 0.19 | 0.42 | 0.00e+0 | 2.24e-16 |
| SEPARATED | system(prefilter+LM) | 1000 | 100.00 | 100.0 | 32.0 | 3.73 | 9.37 | 4.39e+0 | 8.98e+0 |
| SEPARATED | LM(no prefilter) | 1000 | 100.00 | 100.0 | 15.5 | 0.73 | 10.18 | 4.56e+0 | 1.08e+1 |
| SEPARATED | GN(no prefilter, baseline) | 1000 | 100.00 | 100.0 | 2.4 | 0.14 | 0.32 | 4.81e+0 | 1.03e+1 |
| NEAR | system(prefilter+LM) | 400 | 100.00 | 100.0 | 36.8 | 1.38 | 10.24 | 2.95e-2 | 5.75e-2 |
| NEAR | LM(no prefilter) | 400 | 100.00 | 100.0 | 36.3 | 1.19 | 10.75 | 2.95e-2 | 5.75e-2 |
| NEAR | GN(no prefilter, baseline) | 400 | 100.00 | 100.0 | 3.1 | 0.16 | 0.40 | 3.94e-2 | 5.99e-2 |

## 表2：规模扩展（system）

| degree | nControlPoints | n | avgIter | p50Time(ms) | maxTime(ms) |
| --- | --- | --- | --- | --- | --- |
| 2 | 6 | 40 | 2.3 | 0.70 | 1.70 |
| 2 | 10 | 40 | 2.0 | 0.74 | 1.27 |
| 2 | 15 | 40 | 2.3 | 0.84 | 1.64 |
| 3 | 6 | 40 | 2.9 | 0.69 | 1.26 |
| 3 | 10 | 40 | 2.3 | 0.96 | 2.33 |
| 3 | 15 | 40 | 2.8 | 1.05 | 2.71 |
| 5 | 6 | 40 | 2.4 | 1.00 | 1.99 |
| 5 | 10 | 40 | 2.2 | 1.11 | 3.03 |
| 5 | 15 | 40 | 3.0 | 1.29 | 3.62 |
| 7 | 10 | 40 | 2.1 | 1.45 | 3.76 |
| 7 | 15 | 40 | 1.8 | 1.61 | 2.26 |

---
*由 experiment_benchmark.mjs 生成，可复现：node experiment_benchmark.mjs*
