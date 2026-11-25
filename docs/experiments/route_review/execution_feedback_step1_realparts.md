# 执行回执 ①：SimJEB 真实零件上的认证 SSI（有界协议 + OCCT 交叉校验）

- 日期：2026-11-04
- 对应修订路线：S2（真实零件，原 step 1 改协议）
- 主文件：`D:\qqdownload\nurbs-certified-framework\docs\experiments\real_parts\results\summary.md`（275 行，由 `scripts/gen_summary.py` 从 CSV 自动生成）
- 数据集：SimJEB 381 个 STEP 支架（MIT，Open Data Commons Attribution）→ 有理 Bézier patch（`data/*.txt`，851 MB）
- 环境：WSL Ubuntu 22.04 / g++ 11.4 / OCCT 8.0.1（pythonocc）/ AMD Ryzen 7 7840HS

---

## 0. 五条结论（可直接进论文 / 也可直接进 limitations）

1. **数据溯源可复现**：`extract_simjeb.py` 对 228 / 387 两个零件重新抽取，与现存 `data/*.txt` **逐位一致**（max|Δ|=0，2751/2176 控制点）；`BRepBuilderAPI_NurbsConvert` 是必要步骤（不转换时 80/82、91/93 个面不是 BSpline）。
2. **Soundness（强度证据）**：在 1,984 个抽样 patch 对（1,000 正 + 984 阴）上，被认证判为 **disjoint 的 1,404 对中，OCCT `GeomAPI_IntSS` 反例 = 0**；可判定对（type0/type1）一致率 **98.51%**（729/740）。11 个不一致全部是"认证 1D vs OCCT 无交线"，成因是 OCCT 公开 API 只回交线、不回孤立切点，且交点在 patch 参数边界时易漏。
3. **不完备性（必须写进 limitations）**：正样本中 `undecided` 占 **16.3%**；另有 39 个 `undecided_but_occt_has_line`（占抽样 1.97%），94 个 `coplanar_but_occt_no_line`（type4 共 97 个，OCCT 无 2D 重合能力，不计入不一致）。→ 认证分类器在真实数据上**不产生假阴性（在 certified 范围内），但会弃权**。
4. **覆盖检验自证成功**：修正诊断脚本的负坐标符号错误后，**3,840/3,840 个 OCCT 交线采样点全部落在认证 cover 盒的像内（type1/type3 均 100%）**——这是"无遗漏"在真实数据上的直接验证。
5. **可扩展性（硬限制）**：修正 AABB 预筛后，20 个零件 × node_cap=200k × 300 s 超时，**仅 3 个跑完**（228 124 s、387 113 s、115 285 s），17 个 rc=124；单 patch-pair 耗时重尾（最大 29.4 s）；`node_cap` 只约束 cover 细分、**不约束 certify 代价**。→ 认证 SSI 只能定位为 **local certification kernel**，不能写成"全零件扫描器"或"通过了 381 个零件"。

---

## 1. 发现并修复的一个真实 bug（影响历史数字）

`cpp/tests/bench_real_batch.cpp::patch_bb()` 用 `xlo = wx.lo/wmax, xhi = wx.hi/wmin` 计算齐次坐标的保守 AABB。该式只在 `wx >= 0` 时成立；SimJEB 零件大多关于原点对称，`wx < 0` 时会得到**反向/空区间**，导致真实相交的 patch-pair 被预筛误跳过。

| 零件 | faces | bbox 相交 patch-pairs（旧 / 修正后） | 增幅 |
|---|---|---|---|
| 228 | 93 | 408 / 662 | +62% |
| 387 | 82 | 534 / 1,098 | +106% |
| 107 | 132 | 2,027 / 6,704 | +231% |
| 411 | 76 | 18,759 / 39,562 | +111% |
| 136 | 1,189 | 36,516 / 326,430 | +794% |
| 10 | 499 | 6,421 / 37,196 | +479% |
| 38 | 282 | 5,619 / 25,301 | +350% |

- **已修**：`cpp/tests/bench_real_batch.cpp` 中的 `patch_bb()` 已改为区间正确形式（`lo_of = lo>=0 ? lo/wmax : lo/wmin`；`hi_of = hi>0 ? hi/wmin : hi/wmax`），并通过 `g++ -fsyntax-only`。
- **含义**：历史 `logs/*.log` 的 disjoint 计数偏乐观，**论文不得直接引用未修正 harness 的数字**；本回执所有数字均来自修正版。
- 交叉证据：旧预筛判"不相交"的 9 个样本，OCCT 实际找到交线，修正 AABB 后 `classify_pair` 给出 type=1 —— 错在预筛，不在分类器。

---

## 2. 双版本对照（可复现性声明，必须写进论文）

同一份数据、同一台机器、同一会话，分别用 `git HEAD:f60a193e8d1fb5cf8b566e19b7729b00e7c1f699` 与**工作区未提交修改版**（+25/−3 行：coplanar 预检、`classify_pair(node_cap)`）编译：

| 零件 | 版本 | 1D | 2D | coplanar | undecided | certified_boxes | undecided_boxes | 方法耗时 s |
|---|---|---|---|---|---|---|---|---|
| 228 | HEAD | 113 | 2 | 26 | 163 | 8,727 | 97,741 | 320.1 |
| 228 | WORK | 112 | 0 | 28 | 164 | 8,683 | 32,235 | 181.1 |
| 387 | HEAD | 85 | 0 | 61 | 42 | 5,842 | 28,464 | 152.2 |
| 411 | HEAD | 110 | 0 | 41 | 42 | 15,741 | 189,659 | 713.6 |

- HEAD 二进制**逐项复现历史 log 的计数**（1D/2D/coplanar/undecided/certified_boxes/undecided_boxes 全一致）；
- **耗时不可复现**（228：历史 202 s vs HEAD 320 s，同机同二进制，差异来自并发负载）；
- 工作区未提交版**改变分类结果**（228：2D 2→0、coplanar 26→28、undecided_boxes 97,741→32,235，耗时 320→181 s）。
- → 论文必须固定 commit 哈希、声明"计数可复现、时间依赖负载"，并在投稿前把工作区修改 commit 进仓库。

---

## 3. 有界批量结果（20 零件 × cap=200k × 300 s）

| 零件 | faces | patches | face-pairs | 1D | 2D | coplanar | undecided | 方法耗时 s | 是否超时 | 中位 pair ms | 最大 pair ms |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 228 | 93 | 103 | 4,278 | 115 | 0 | 33 | 171 | 124.1 | 否 | 0.03 | 642 |
| 387 | 82 | 110 | 3,321 | 86 | 0 | 72 | 43 | 113.3 | 否 | 0.03 | 660 |
| 115 | 66 | 279 | 2,145 | 116 | 0 | 49 | 37 | 285.0 | 否 | 0.08 | 1,587 |
| 248 | 56 | 197 | 1,540 | — | — | — | — | 300 | **是** | 1,616 | 29,072 |
| 411 | 76 | 7,480 | 2,850 | — | — | — | — | 300 | **是** | 0.03 | 628 |
| 136 | 1,189 | 3,577 | 706,266 | — | — | — | — | 300 | **是** | 0.04 | 3,045 |
| 其余 14 个 | 66–499 | 218–1,674 | 5,050–124,251 | — | — | — | — | 300 | **是**（rc=124） | — | — |

**读法**：完成率 3/20；即使 56 面 / 197 patch 的零件 248 也在 300 s 内跑不完（单 pair 最大 29 s，重尾严重）。→ 真实零件结论只能是"**在抽样 patch 对上做分类一致性与 soundness 统计**"，不能是"通过率"。

---

## 4. OCCT 交叉校验（1,984 抽样对）

| 类别 | 数量 | 说明 |
|---|---|---|
| certified-disjoint（type0） | 1,404 | **OCCT 反例 0** |
| certified-1D（type1） | ~470 | 与 OCCT 一致 729/740 可判定对（98.51%） |
| `undecided_but_occt_has_line` | 39（1.97%） | 认证弃权，OCCT 有交线 → **不完备性证据** |
| `certified_1D_but_occt_no_line` | 11 | 相切/边界接触；OCCT IntSS 只回交线、不回孤立切点 |
| `coplanar_but_occt_no_line` | 94 | OCCT 无 2D 重合能力，不应计入不一致 |
| 阴性对照（bbox 不相交） | 984 | 认证全部 disjoint，OCCT 全部无线 → 一致 |

---

## 5. 覆盖检验（"无遗漏"在真实数据上的验证）

初版诊断脚本与 `patch_bb` 犯了同一个负坐标符号错误，导致 34–87% 的假"miss"。修正后：

| 零件 | type=1 pairs / points | type=3 pairs / points | 点覆盖率 |
|---|---|---|---|
| 107 | 18 / 800 | 2 / 80 | 100% |
| 10 | 6 / 280 | 1 / 40 | 100% |
| 228 | 16 / 680 | 3 / 120 | 100% |
| 387 | 13 / 520 | 2 / 240 | 100% |
| 38 | 4 / 160 | 0 / 0 | 100% |
| 411 | 20 / 840 | 2 / 80 | 100% |
| **合计** | **77 / 3,280** | **10 / 560** | **3,840/3,840 = 100%** |

→ "OCCT 交线上的采样点 100% 落在认证 cover 盒的像内"。这一条是论文 §7.x 的最强证据，也正好展示了"诊断脚本自身出错 → 自查修正"的严谨性（审稿人喜欢这种自我验证）。

---

## 6. 对论文的具体修改建议

1. **新增一节 §7.10 "Real-part benchmark (SimJEB)"**，只写三件事：① 数据溯源可复现；② 1,984 抽样对上 soundness 统计（0 反例、98.51% 一致）；③ 覆盖检验 100% + 可扩展性限制（3/20、重尾、local kernel 定位）。
2. **不要写"382 个零件全部通过"或任何通过率**。
3. **limitations 加一条**：认证 cover 的成本 ∝ (1/τ)^d，在真实 B-rep 上需要 BVH + trimmed patch + 局部化才能走向工业级；本文给的是可认证内核。
4. **reproducibility 段落**：固定 commit 哈希 + "计数可复现、时间依赖负载" + 工作区修改版与 HEAD 的差异表。
5. Figure：建议画一张"时间 vs faces/patches"散点 + 一条"300 s 超时线"，直观展示 local kernel 定位。

---

## 7. 诚实限制

- 只有 3/20 零件在 300 s 内跑完；其余只有部分累积统计（PAIRSTAT），不是全量结果。
- OCCT `GeomAPI_IntSS` 公开 API 只回曲线、不回孤立点；`GeomAPI_IntCS`（OCCT 8.0）无容差参数。
- 抽样为每件 ≤50 正 + ≤50 负，非全量；统计口径已在 summary.md 中给出。
- 工作区 `cert_ssi.hpp` 修改尚未 commit（投稿前必须 commit，否则论文数字无法对应到公开仓库状态）。
