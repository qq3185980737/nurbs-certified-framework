# rev19 → rev20 修改清单（2026-11-04）

- 基准稿：`01_文章正文/P1_CAD_rev19_清理版.docx`（189 段 / 11 表 / 31 refs / 20 页）
- 新稿：`01_文章正文/P1_CAGD_rev20_草案.docx`（209 段 / 12 表 / 41 refs）
- 所有改动均在原文基础上**插入或追加**，未删除任何原有段落与表格；公式、图、表 1–11 原样保留。

## 一、结构与正文改动

| # | 位置 | 改动 | 目的 |
|---|---|---|---|
| 1 | §1.1 相关工作 | 新增 1 段（约 1,800 字符）：矩阵表示线 [29]→[32]、自交线 [33]/[18]、重叠线 [34]、认证线 [36][37][38]、相邻片决策 [35]、小环/相切 [39]、最近点 [40] | 回应主编"未引用 SOTA"；把真正的 SOTA（Jia/Cheng/Yang 组）接进来 |
| 2 | §2.5（新增，Heading 2） | 证书保证等级：定理级（2×2 Krawczyk）/ 包含级（graph-patch cover）/ 证据级（仿射共面重叠、κ 标签、退化分类）；并说明 `Y` 是任意实矩阵，Krawczyk 不要求区间逆 | 解决 "certified" 口径不一致这一最容易被 CAGD 审稿人攻击的点 |
| 3 | 贡献 (4) | 追加：共面重叠为 evidence-level 且仅限仿射片，一般 NURBS 重叠由 [34] 解决，本文不声称 | 主动降级会与 Yang & Jia TOG 2025 冲突的卖点 |
| 4 | 贡献 (5) | 追加：§2.5 三级约定贯穿 §3–§7，定理级与证据级不混用 | 同上 |
| 5 | §6.2 | 追加共面重叠的适用范围边界，明确引 [34] | 防止 over-claim |
| 6 | §7.8 定量对比 | 新增 1 段：SISL 4.7 / OCCT 8.0.1 / 矩阵法 三工具在 A–F 同 fixture 的三段式结论（横截一致；C 角点 SISL 优于 OCCT；D 退化切触 OCCT 漏、SISL 拓扑碎片化、认证 85,344 盒 145 ms） | 直接回应 "benchmarked against"；把"别人算不出来"改成"别人没有证书" |
| 7 | §7.8 caveats | "Two caveats" → "Three caveats"：新增"曲面级 SISL/OCCT 结果是 geometry-only，不含存在唯一性证书，因此比较的是输出的完备性与稳健性，不是保证强度" | 可比性边界 |
| 8 | §7.10（新增，Heading 2） | 真实零件基准：SimJEB 381 件溯源、1,984 抽样 patch 对、1,404 个 certified-disjoint 零反例、729/740（98.51%）一致、39（1.97%）弃权、cover 覆盖 3,840/3,840（100%）、3/20 在 300 s 内完成、单对重尾 | 回应 "examples on challenging geometries"，同时如实交代可扩展性边界 |
| 9 | Table 12（新增） | 真实零件基准汇总表（11 行 × 2 列，Table Grid） | 同上 |
| 10 | §8 结论 | 追加真实零件结论 + "local certification kernel" 定位 + §2.5 保证等级声明 | 与新增实验一致 |
| 11 | Abstract | 末尾追加 1 句（约 40 词）真实零件 soundness 头条 | 让主编在第一段就看到新证据 |
| 12 | Data availability | 追加 SimJEB/基线产物路径、commit `1ff0e1f`、per-pair 计数可复现而 wall-clock 依赖负载 | 可复现性声明 |

## 二、参考文献改动

- **修正 [19]**：原 "Li G, Wu X." → "Li G, Wu X, Luo Z, Na L, Cheng JS."，补 DOI 与勘误 DOI 10.1145/3839107。
- **新增 [32]–[40]**（9 条，均在正文被引用，无孤立引用）：
  - [32] Jia, Li, Cheng, CAD 2022（矩阵表示法现状）
  - [33] Park, Hong, Kim, Elber, CAGD 2021（自交最新后续）
  - [34] Yang & Jia, TOG 2025（重叠域提取）
  - [35] Liu et al., CAGD 2026（相邻片求交检测）
  - [36] Breiding, Rose, Timme, TOMS 2023（区间认证零点）
  - [37] Elber, CAGD 2022（保守核）
  - [38] Wang et al., CGF 2025（watertightness）
  - [39] Yang & Jia, TOG 2026（小环/相切）
  - [40] Botbol, Busé, Chardin, Yildirim, SIAA 2020（正交投影/最近点）
- **新增 [41]**：Whalen, Beyene, Mueller, SimJEB, arXiv:2105.03534（数据集引用）。
- 参考文献连续 [1]–[41]，无重复。

## 三、配套产物

- `03_投稿配套材料/Cover_letter_CAGD_rev20.docx`：含"对 CAD 拒稿意见的逐条回应"与"与最近先验工作的关系"。
- 仓库 `docs/experiments/`：route_review / literature_2026 / real_parts / baselines 四个目录 + README。
- 仓库 commit：`1ff0e1f`（代码与全部实验产物）、`97de9a7`（README）。

## 四、尚未处理（需作者决策）

1. 是否按评估报告 §4(c) **收窄主线**（把曲线/最近点/重叠降为实例）——会改标题、摘要、贡献列表。
2. 目标期刊最终确认（CAGD vs CAD 重投）。
3. Table 9 仍为 rev19 的 16 行定性表，未按新增文献扩行（建议 S5 阶段扩到 18–20 行）。
4. 摘要 292 词，若 CAGD 要求 ≤250 需再压缩（可删 §2.5 相关从句）。
5. 图件（Fig. 1–5）未更新；新增真实零件建议补 1 张"时间 vs faces"散点图。
