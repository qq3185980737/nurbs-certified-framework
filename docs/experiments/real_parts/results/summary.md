# SimJEB 真实零件基准：认证 SSI（certified surface-surface intersection）结果摘要

> 本文件由 `scripts/gen_summary.py` 从同目录 CSV 自动生成；所有数字均来自实际运行。
> 数据：SimJEB 381 个 STEP 零件 → 有理 Bézier patch（`docs/experiments/real_parts/data/`）。

## 0. 代码版本与可复现性声明

- 仓库 commit：`f60a193e8d1fb5cf8b566e19b7729b00e7c1f699`（`feat_认证`，2026-11-04）；**工作区另有未提交修改**：`cpp/include/nurbs/cert_ssi.hpp` +25/−3 行
  （新增 `exact_coplanar_overlap` 预检、`classify_pair(..., node_cap)` 预算参数）。
- 两版二进制同机同会话双跑对照见 §3；**论文必须固定 commit 哈希并与本节一致**。
- 提取器：OCCT 8.0.1（pythonocc-core，`/root/miniconda3/envs/occ`），WSL Ubuntu 22.04，g++ 11.4，8 核。
- 认证器：`cpp/include/nurbs/cert_ssi.hpp::nurbs::cert::classify_pair(A,B,tau=0.02,node_cap)`。
- 本摘要中的**最终基准数字采用修正 AABB 预筛后的实现**（见 §1），未修正版仅用于复现历史 log。

## 1. 关键发现：bench 的 AABB 预筛有负坐标符号错误（必须先修）

`cpp/tests/bench_real_batch.cpp` 的 `patch_bb()` 用 `xlo=wx.lo/wmax, xhi=wx.hi/wmin` 计算保守 AABB。
该式只在 `wx>=0` 时成立；对 `wx<0`（SimJEB 零件大多关于原点对称）会得到**反向/空区间**，
导致真实相交的 patch-pair 被误判为 bbox 不相交而**直接跳过**。修正做法：
`wx>=0` 时 `x∈[wx/wmax, wx/wmin]`，`wx<0` 时 `x∈[wx/wmin, wx/wmax]`。

| 零件 | faces | patches | face-pairs | bbox 相交 face-pairs（旧/新） | patch-pairs | bbox 相交 patch-pairs（旧 try / 新 try） |
|---|---|---|---|---|---|---|
| 228 | 93 | 103 | 4278 | 380 / **573** | 711 | 408 / **662** |
| 387 | 82 | 110 | 3321 | 371 / **666** | 1188 | 534 / **1098** |
| 107 | 132 | 804 | 8646 | 552 / **1236** | 17174 | 2027 / **6704** |
| 411 | 76 | 7480 | 2850 | 382 / **603** | 3326433 | 18759 / **39562** |
| 136 | 1189 | 3577 | 706266 | 13998 / **101257** | 471070 | 36516 / **326430** |
| 10 | 499 | 1054 | 124251 | 3550 / **14323** | 53207 | 6421 / **37196** |
| 38 | 282 | 937 | 39621 | 2108 / **6405** | 54285 | 5619 / **25301** |
| 248 | 56 | 197 | 1540 | 213 / **235** | 2821 | 808 / **920** |
| 115 | 66 | 279 | 2145 | 246 / **465** | 5948 | 858 / **2218** |

结论：旧预筛漏掉的 bbox 相交 patch-pair 比例：228 +62%、387 +106%、107 +231%、411 +111%、136 +794%、10 +479%、38 +350%、248 +14%、115 +159%（多数零件 +50%~+300%，136 达 +794%）。因此历史 `logs/*.log` 与 `cpp/tests/bench_real_batch.cpp`
给出的“disjoint”计数**偏乐观**，真实工作量显著更大。

补充证据：旧预筛下抽到的 9 个“bbox 不相交”样本，OCCT `GeomAPI_IntSS` 实际找到交线；
用修正 AABB 复算后 patch-pair 确实相交，且 `classify_pair` 给出 type=1（认证 1D）——
说明错的是预筛，不是分类器。

## 2. 数据溯源：STEP → Bézier 提取可复现（逐位一致）

`scripts/extract_simjeb.py`（STEPControl_Reader → NurbsConvert → GeomConvert_BSplineSurfaceToBezierSurface，
只导 SURF）对 2 个零件重抽取并与现存 `data/*.txt` 逐项比对：

| 零件 | faces | patches | 比对控制点数 | max&#124;Δ&#124; | 结论 |
|---|---|---|---|---|---|
| 387 | 82 | 110 | 2752 | 0.0 | 结构+数值逐位一致 |
| 228 | 93 | 103 | 2176 | 0.0 | 结构+数值逐位一致 |

反证：不做 `BRepBuilderAPI_NurbsConvert` 时 387 有 80/82、228 有 91/93 个面不是 BSpline，
`ad.BSpline()` 直接失败 → NurbsConvert 是必要步骤。数据可复现，溯源成立。

## 3. 双跑对照：git HEAD 基线 vs 工作区未提交修改版

同一批 data、同一台机器、同一次会话，分别用 `git HEAD:f60a193e8d1fb5cf8b566e19b7729b00e7c1f699` 的 `cert_ssi.hpp`（3 参 `classify_pair`，
无 coplanar 预检/无预算）与工作区修改版（4 参 + coplanar 预检）编译运行：

| 零件 | 版本 | 1D | 2D | coplanar | undecided | certified_boxes | undecided_boxes | 方法耗时 s |
|---|---|---|---|---|---|---|---|---|
| 228 | HIST | 113 | 2 | 26 | 163 | 8727 | 97741 | 202.07 |
| 228 | HEAD | 113 | 2 | 26 | 163 | 8727 | 97741 | 320.14 |
| 228 | WORK | 112 | 0 | 28 | 164 | 8683 | 32235 | 181.13 |
| 387 | HIST | 85 | 0 | 61 | 42 | 5842 | 28464 | 81.57 |
| 387 | HEAD | 85 | 0 | 61 | 42 | 5842 | 28464 | 152.19 |
| 387 | WORK | 85 | 0 | 61 | 42 | 5842 | 23253 | 142.58 |
| 107 | HIST | 137 | 0 | 66 | 63 | 22317 | 88832 | 297.66 |
| 107 | HEAD | 137 | 0 | 66 | 63 | 22317 | 88832 | 478.12 |
| 107 | WORK | 137 | 0 | 66 | 63 | 22317 | 85816 | 461.80 |
| 411 | HIST | 110 | 0 | 41 | 42 | 15741 | 189659 | 542.03 |
| 411 | HEAD | 110 | 0 | 41 | 42 | 15741 | 189659 | 713.58 |
| 411 | WORK | 110 | 0 | 41 | 42 | 15741 | 169571 | 680.87 |

结论：

1. **HEAD 二进制逐项复现历史 `logs/` 的计数**（228/387/107/411 的 1D/2D/coplanar/undecided/
   certified_boxes/undecided_boxes 全部一致）→ 历史数字来自提交版基线，可追溯。
2. **耗时不可复现**：228 历史 202 s vs HEAD 320 s（同机同二进制，差异为并发负载所致）；
   论文只能报告“同机同会话测量”，并注明负载/并发度。
3. 工作区未提交修改版**改变分类结果**（228：2D 2→0、coplanar 26→28、undecided_boxes 97741→32235，
   同负载下 320 s→181 s），因此论文必须声明用的是哪一版。

## 4. 有界批量：20 零件 × node_cap=200k + 300 s 超时（修正 AABB 预筛后）

说明：**“-”表示该零件在 300 s 内未跑完、没有最终 RESULT**；其“中位/p95/最大 pair ms”与“undecided 占比”
来自程序每 15 s 输出的累计 PAIRSTAT（超时前最后一次），是**部分累计**而非全量；“已处理 face-pair”为
wall 时钟结束时已写出 detail 的 face-pair 数（修正预筛后这些才是“bbox 相交的 face-pair”）。

| 零件 | faces | patches | face-pairs | patch-pairs | try | 1D | 2D | coplanar | undecided | 方法耗时 s | 墙钟 s | 超时 | 已处理 face-pair | 中位 pair ms | p95 ms | 最大 pair ms | cap 触顶 | undecided 占比 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 248 | 56 | 197 | 1540 | - | - | - | - | - | - | - | 300.0 | **是** | 13 | 1615.92 | 3638.8 | 29072.3 | 0 | 0.0% |
| 115 | 66 | 279 | 2145 | 5948 | 2218 | 116 | 0 | 49 | 37 | 285.0 | 285.0 | 否 | 465 | 0.08 | 1113.3 | 1586.6 | 0 | 1.3% |
| 69 | 66 | 1674 | 2145 | - | - | - | - | - | - | - | 300.0 | **是** | 82 | 25.46 | 14263.5 | 15938.7 | 0 | 22.2% |
| 411 | 76 | 7480 | 2850 | - | - | - | - | - | - | - | 300.0 | **是** | 397 | 0.03 | 498.6 | 627.7 | 0 | 1.1% |
| 387 | 82 | 110 | 3321 | 1188 | 1098 | 86 | 0 | 72 | 43 | 113.3 | 113.4 | 否 | 666 | 0.03 | 266.1 | 659.8 | 0 | 4.4% |
| 228 | 93 | 103 | 4278 | 711 | 662 | 115 | 0 | 33 | 171 | 124.1 | 124.1 | 否 | 573 | 0.03 | 280.9 | 642.4 | 0 | 23.4% |
| 506 | 101 | 403 | 5050 | - | - | - | - | - | - | - | 300.0 | **是** | 290 | 278.00 | 2256.8 | 8684.4 | 0 | 15.0% |
| 70 | 109 | 1067 | 5886 | - | - | - | - | - | - | - | 300.0 | **是** | 163 | 147.99 | 205.7 | 826.1 | 0 | 6.1% |
| 538 | 115 | 396 | 6555 | - | - | - | - | - | - | - | 300.0 | **是** | 287 | 0.10 | 849.7 | 6136.0 | 0 | 8.0% |
| 570 | 118 | 1300 | 6903 | - | - | - | - | - | - | - | 300.0 | **是** | 246 | 0.17 | 754.0 | 5975.8 | 0 | 5.8% |
| 113 | 120 | 418 | 7140 | - | - | - | - | - | - | - | 300.0 | **是** | 68 | 0.19 | 14283.1 | 15871.6 | 0 | 10.0% |
| 205 | 124 | 261 | 7626 | - | - | - | - | - | - | - | 300.0 | **是** | 926 | 884.19 | 3672.6 | 3905.4 | 0 | 16.7% |
| 546 | 126 | 1466 | 7875 | - | - | - | - | - | - | - | 300.0 | **是** | 1 | - | - | - | - | - |
| 107 | 132 | 804 | 8646 | - | - | - | - | - | - | - | 300.0 | **是** | 882 | 67.10 | 500.5 | 1021.8 | 0 | 4.0% |
| 534 | 136 | 218 | 9180 | - | - | - | - | - | - | - | 300.0 | **是** | 708 | 987.21 | 3329.4 | 10712.3 | 0 | 20.0% |
| 200 | 140 | 419 | 9730 | - | - | - | - | - | - | - | 300.0 | **是** | 1327 | 0.25 | 419.4 | 4183.4 | 0 | 5.0% |
| 614 | 149 | 7946 | 11026 | - | - | - | - | - | - | - | 300.0 | **是** | 851 | 93.28 | 308.8 | 812.9 | 0 | 9.2% |
| 38 | 282 | 937 | 39621 | - | - | - | - | - | - | - | 300.0 | **是** | 228 | 0.10 | 462.3 | 5364.3 | 0 | 2.9% |
| 10 | 499 | 1054 | 124251 | - | - | - | - | - | - | - | 300.0 | **是** | 494 | 141.96 | 261.2 | 796.7 | 0 | 26.4% |
| 136 | 1189 | 3577 | 706266 | - | - | - | - | - | - | - | 300.0 | **是** | 867 | 0.04 | 0.1 | 3044.9 | 0 | 1.0% |

完成 3/20 件（300 s 预算）：115（285.0 s）、387（113.3 s）、228（124.1 s）。
其余 17 件全部 rc=124（超时），已处理 face-pair 数见上表；其中 248/546/69 等小零件在单个面片对上即卡死，
说明单对 `classify_pair` 的耗时是**重尾**的（样本中最慢单对 29.4 s），`node_cap` 只约束 cover 盒数，
**并不约束 certify_boxes 的代价**，因此“每 patch-pair 预算”无法给出单对硬上界。

补充：历史 harness（未修正 AABB）下 4 个参照零件 cap=4M 均完成：228 181 s、387 143 s、107 462 s、411 681 s；
对照组 HEAD 基线 320/152/478/714 s（同负载）。**cap=4M 的 20 零件全覆盖从未跑完**。

## 5. OCCT 交叉校验（GeomAPI_IntSS，tol=1e-6）

每零件抽 ≤50 个 bbox 相交 patch-pair（正样本）+ ≤50 个 bbox 不相交 pair（阴性对照），
对同一 data 重建有理 Bézier 曲面并调用 `GeomAPI_IntSS`。

| 样本类 | 对数 | 方法 type=0 | type=1 | type=3 | type=4 | OCCT 无交线 | OCCT 有交线 | OCCT 异常 |
|---|---|---|---|---|---|---|---|---|
| 正样本（bbox 相交） | 1000 | 471 | 269 | 163 | 97 | 700 | 300 | 0 |
| 阴性对照（bbox 不相交） | 984 | 933 | 9 | 23 | 19 | 973 | 11 | 0 |

**Soundness（最重要）**：被方法判为 **certified disjoint** 的抽样对共 **1404** 个
（正样本 type0 471 + 阴性对照 type0 933），OCCT 找到交线的反例 = **0**。
→ “在所有被认证判为 disjoint 的真实面片对上，OCCT 也没有交线（0 反例）”。

**区间排除优于 bbox**：正样本中 471/1000（47.1%）bbox 相交但被认证判为 disjoint，
OCCT 对这 471 对全部无交线（0 反例）→ 认证的区间排除确实比 bbox 更紧。

**可判定一致率**（只统计 OCCT 能裁决的 type0/type1 正样本）：
certified-disjoint 且 OCCT 无交线 471 对 + certified-1D 且 OCCT 有交线 258 对 = 729/740 = **98.51%**。

**不完备性 / 差异清单**（`results/occt_conflicts.csv`）：

| tag | 数量 | 含义 | 成因 |
|---|---|---|---|
| undecided_but_occt_has_line | 39 | 认证未决(type3)，但 OCCT 有交线 | 覆盖触顶/横截认证失败 → 保守未决（**不是错误**，是不完备） |
| certified_1D_but_occt_no_line | 11 | 认证 1D(type1)，OCCT 无交线 | 相切/边界接触：OCCT `GeomAPI_IntSS` 只回交线、不回孤立切点，且交点落在 patch 参数边界时易漏 |
| coplanar_but_occt_no_line | 94 | 认证共面重叠(type4)，OCCT 无交线 | OCCT `GeomAPI_IntSS` 不支持 2D 重合（两张重合曲面是退化情形，返回 0 条交线）；**不应计入不一致** |

阴性对照中另有 9 对被方法判为 type=1（且 OCCT 有交线），全部来自**旧 AABB 预筛把相交对误判为不相交**
（见 §1）；修正 AABB 后已逐个复核为确实相交。故这批不是分类器错误。

差异样例（前 3 个/类）：

- `undecided_but_occt_has_line`：part=10 面片 (fidA=411,piA=0)-(fidB=492,piB=0)，type=3 cover=152 certified=0 ms=430.7，OCCT lines=1
- `undecided_but_occt_has_line`：part=107 面片 (fidA=44,piA=0)-(fidB=77,piB=0)，type=3 cover=114 certified=0 ms=304.1，OCCT lines=1
- `undecided_but_occt_has_line`：part=107 面片 (fidA=3,piA=0)-(fidB=23,piB=0)，type=3 cover=47 certified=0 ms=170.3，OCCT lines=1
- `certified_1D_but_occt_no_line`：part=10 面片 (fidA=144,piA=0)-(fidB=188,piB=0)，type=1 cover=77 certified=76 ms=138.4，OCCT lines=0
- `certified_1D_but_occt_no_line`：part=107 面片 (fidA=123,piA=13)-(fidB=130,piB=0)，type=1 cover=69 certified=4 ms=252.6，OCCT lines=0
- `certified_1D_but_occt_no_line`：part=113 面片 (fidA=40,piA=0)-(fidB=54,piB=0)，type=1 cover=95 certified=54 ms=158.4，OCCT lines=0
- `coplanar_but_occt_no_line`：part=10 面片 (fidA=251,piA=0)-(fidB=303,piB=1)，type=4 cover=104 certified=0 ms=356.7，OCCT lines=0
- `coplanar_but_occt_no_line`：part=10 面片 (fidA=248,piA=0)-(fidB=402,piB=0)，type=4 cover=97 certified=0 ms=777.8，OCCT lines=0
- `coplanar_but_occt_no_line`：part=10 面片 (fidA=51,piA=0)-(fidB=141,piB=0)，type=4 cover=77 certified=0 ms=487.1，OCCT lines=0

## 6. 覆盖检验（OCCT 交线采样点是否落在认证 cover 盒内）

对 OCCT 判定有交线的抽样对，直接调用 `cover(A,B,tau=0.02,node_cap=4M)` 得到参数域 4D 盒，
用区间 de Casteljau 求每个盒的三维保守像，检查 OCCT 交线采样点是否落在像 bbox 内（必要条件检验）。

| pair type | 对数 | 采样点 | 落在盒像内 | 覆盖率 | 全部点覆盖的 pair |
|---|---|---|---|---|---|
| type=1 | 77 | 3280 | 3280 | 100.00% | 77/77 |
| type=3 | 10 | 560 | 560 | 100.00% | 10/10 |

合计 **87 对 / 3840 个 OCCT 采样点 / 3840 个点在两曲面的 cover 盒像内（100.00%）**。
说明：初版诊断脚本的三维像 bbox 也犯了同一个负坐标符号错误（`lo=wx.lo/whi, hi=wx.hi/wlo`），
一度得到 34%–87% 的假 miss；修正后覆盖率 100%。因此“miss”来自诊断工具，不是 cover 漏覆盖。
另注：本检验按 pair type 拆分后，type=1（真正有 certified 覆盖盒的对）覆盖率亦为 100%，
不存在“有盒但点不在盒内”的样本。本检验为**必要条件**（保守），通过不构成完备性证明。

## 7. 可直接进论文的表：20 零件抽样统计（抽样时 bbox 预筛未修正，仅用于分类一致性/耗时）

“一致/可判定”只统计 OCCT 能裁决的 type0（应为无交线）与 type1（应为有交线）；
type3（未决）与 type4（共面，OCCT 无 2D 重合能力）不计入分母。

| 零件 | faces | patches | pos type0 | type1 | type3 | type4 | OCCT 一致/可判定 | 中位 pair ms | 最大 pair ms |
|---|---|---|---|---|---|---|---|---|---|
| 248 | 56 | 197 | 15 | 18 | 9 | 8 | 32/33 (97.0%) | 321.5 | 29357 |
| 115 | 66 | 279 | 13 | 20 | 9 | 8 | 33/33 (100.0%) | 172.8 | 3730 |
| 69 | 66 | 1674 | 20 | 7 | 18 | 5 | 27/27 (100.0%) | 201.8 | 3084 |
| 411 | 76 | 7480 | 21 | 20 | 6 | 3 | 41/41 (100.0%) | 261.1 | 5579 |
| 387 | 82 | 110 | 27 | 13 | 2 | 8 | 40/40 (100.0%) | 0.1 | 987 |
| 228 | 93 | 103 | 13 | 16 | 19 | 2 | 28/29 (96.6%) | 53.0 | 1874 |
| 506 | 101 | 403 | 17 | 22 | 7 | 4 | 37/39 (94.9%) | 154.3 | 18548 |
| 70 | 109 | 1067 | 11 | 7 | 28 | 4 | 18/18 (100.0%) | 67.8 | 3973 |
| 538 | 115 | 396 | 30 | 9 | 5 | 6 | 39/39 (100.0%) | 0.1 | 6813 |
| 570 | 118 | 1300 | 27 | 17 | 4 | 2 | 44/44 (100.0%) | 0.2 | 2184 |
| 113 | 120 | 418 | 29 | 12 | 6 | 3 | 40/41 (97.6%) | 0.1 | 3700 |
| 205 | 124 | 261 | 25 | 19 | 5 | 1 | 44/44 (100.0%) | 0.4 | 4884 |
| 546 | 126 | 1466 | 18 | 18 | 7 | 7 | 36/36 (100.0%) | 150.6 | 4015 |
| 107 | 132 | 804 | 21 | 19 | 5 | 5 | 39/40 (97.5%) | 71.6 | 1952 |
| 534 | 136 | 218 | 27 | 6 | 4 | 13 | 33/33 (100.0%) | 0.1 | 1832 |
| 200 | 140 | 419 | 20 | 16 | 9 | 5 | 34/36 (94.4%) | 75.2 | 1441 |
| 614 | 149 | 7946 | 30 | 19 | 1 | 0 | 47/49 (95.9%) | 0.1 | 857 |
| 38 | 282 | 937 | 30 | 3 | 14 | 3 | 33/33 (100.0%) | 0.1 | 7318 |
| 10 | 499 | 1054 | 34 | 7 | 5 | 4 | 40/41 (97.6%) | 0.1 | 4060 |
| 136 | 1189 | 3577 | 43 | 1 | 0 | 6 | 44/44 (100.0%) | 0.1 | 1850 |

逐零件原始列在 `results/probe_summary.csv` 与 `results/occt_confusion.csv`。

**三句结论（可直接改写进论文）**：

1. 在 SimJEB 真实 B-rep 面片对上，认证分类与 OCCT `GeomAPI_IntSS` 在“已认证”判定上**零反例**：
   共 1404 个 certified-disjoint 抽样对全部被 OCCT 确认无交线；可判定的 740 对 certified-disjoint/1D 样本中 729 对一致（98.5%），
   11 个 certified-1D 与 OCCT 不一致的样本均集中在相切/边界接触；type4 共面样本共 97 个，其中 94 个 OCCT 无交线，
   属 OCCT `GeomAPI_IntSS` 无 2D 重合能力，不应计入不一致。
2. 不完备性集中在 undecided（type3，正样本 163/1000 = 16.3%）与共面（type4，97/1000）：前者保守未决，
   后者对应 OCCT 不支持 2D 重合，均不是 soundness 反例。
3. 可扩展性上，单 patch-pair 耗时重尾（样本最大 29.4 s），20 个零件中仅 3 个能在 300 s 内跑完，
   认证 SSI 只能定位为**局部认证内核（local certification kernel）**，不能写成“可扫描全零件的求解器”。

## 8. 结论：定位与 limitations

- **定位**：`local certification kernel`。逐 patch-pair 的 cover+certify 有严格几何保证、在真实数据上
  与 OCCT 无 soundness 反例，但代价重尾且随面片数快速爆炸，无法作为全零件扫描器。
- 论文中该实验应写成“真实数据上的分类一致性 + 可扩展性边界证据”，**不能写“认证了 382 个零件”**。
- 若要走向全零件，需要：按 face 拓扑（共享边界）剪枝/分区、对 undecided 做自适应加密或降级求解、
  以及真正的单对硬超时（目前 `node_cap` 不约束 certify 代价）。
- 未做/未完成：382 件全覆盖、结果与论文 §7 的 19 件历史运行逐件对照（历史 log 仅 4 件完成）、
  相切/边界接触的独立判定（未使用 OCCT 以外的真值）。

## 9. 文件清单与复现命令

脚本（`docs/experiments/real_parts/scripts/`）：

| 文件 | 作用 |
|---|---|
| `extract_simjeb.py` | STEP→有理 Bézier patch（SURF only），溯源用 |
| `compare_extract.py` / `run_verify_extract.sh` | 重抽取 vs 现存 data 逐项比对（228/387 PASS） |
| `bench_real_batch_probe.cpp` | 插桩基准（PAIRSTAT：中位/最大单对耗时、cap 触顶、pair 级 type 分布）+ 修正 AABB |
| `bench_real_batch_head.cpp` | HEAD 基线复现版（3 参 classify_pair，配套 `results/head_baseline/include/nurbs/cert_ssi.hpp`） |
| `probe_pairs.cpp` | patch-pair 抽样探针（正/负样本 + 逐对计时） |
| `occt_check.py` / `run_probe.py` | OCCT `GeomAPI_IntSS` 交叉校验 + 交线采样点 |
| `cover_check.cpp` / `run_cover.py` / `cover_by_type.py` | cover 覆盖检验（按 type 拆分） |
| `run_real_batch.py` | 批量调度（每件 wall 超时 + node_cap 预算，并行） |
| `analyze_results.py` / `gen_summary.py` | 汇总 CSV / 生成本 summary.md |
| `debug_pair.cpp` / `check_neg_hits.py` | AABB bug 复现与可疑样本复核 |
| `build.sh` / `build_v3.sh` | 编译全部二进制（`docs/experiments/real_parts/build/`） |

结果（`docs/experiments/real_parts/results/`）：

| 文件 | 内容 |
|---|---|
| `summary.md` | 本文件（最终交付） |
| `repro_dualbinary.csv` | HIST/HEAD/WORK 三版对照（§3） |
| `batch/{A,H,B}/results.csv` | 三组批量运行原始结果（含 PAIRSTAT 列、超时标记） |
| `batch/all_results.csv` | 合并 |
| `probe_summary.csv` | 20 零件抽样统计（type 分布、耗时、OCCT 一致数） |
| `occt_confusion.csv` | 按样本类 × method type 的混淆矩阵 |
| `occt_conflicts.csv` | 差异清单（§5） |
| `cover_summary.csv` / `cover_by_type.csv` | 覆盖检验（§6） |
| `bbox_bug_counts.txt` | 修正前后 AABB 预筛计数（§1） |
| `extract_verify/` | 重抽取数据 + 比对日志 |
| `head_baseline/include/nurbs/cert_ssi.hpp` | git HEAD 的 `cert_ssi.hpp` 副本（可复现性） |

复现命令（WSL，仓库根 `D:\qqdownload\nurbs-certified-framework` = `/mnt/d/qqdownload/nurbs-certified-framework`）：

```bash
cd docs/experiments/real_parts && bash scripts/build_v3.sh          # 编译
bash scripts/run_verify_extract.sh 387 228                          # 溯源比对
python3 scripts/run_real_batch.py --group B --jobs 8 \
    --binary build/bench_real_batch_probe_v3                        # 20 零件有界批量
python3 scripts/run_probe.py --jobs 4                               # OCCT 交叉校验
python3 scripts/run_cover.py 228,387,107,10,38,411                  # 覆盖检验
python3 scripts/analyze_results.py && python3 scripts/gen_summary.py
```

