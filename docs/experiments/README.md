# docs/experiments — 实验与路线评估产物索引（2026-11-04）

本目录是论文「Certified Unified Interrogation of NURBS/Bézier Geometry」（rev20）在 CAD 拒稿、改投 CAGD 过程中新增的实验与评估产物。所有数字均来自本机实际运行，脚本与结果 CSV 一并保留，可按下列命令复现。

## 目录

| 目录 | 内容 |
|---|---|
| `route_review/` | 四步路线评估（`route_assessment_2026-11-04.md`）、证书等级审计（`certificate_levels.md`）、第 ①/②④ 步执行回执 |
| `literature_2026/` | 2020–2026 结构化文献检索（Crossref/OpenAlex，逐条 DOI 可核验），必补引用 12 条 |
| `real_parts/` | SimJEB 381 件 STEP 支架 → 有理 Bézier patch，有界批量认证 + OCCT 交叉校验 + cover 检验；主文件 `results/summary.md` |
| `baselines/` | SISL 4.7 / OCCT 8.0.1 / 矩阵表示法最小复现 的同 fixture 对比；主文件 `comparison.md` |

## 快速复现

```sh
# 真实零件（需要 WSL + OCCT 8.0.1 pythonocc）
bash docs/experiments/real_parts/scripts/run_verify_extract.sh   # 提取溯源（逐位一致性）
python3 docs/experiments/real_parts/scripts/gen_summary.py       # 由 CSV 重新生成 summary.md

# 基线对比
bash docs/experiments/baselines/scripts/build_run_sisl.sh        # SISL A–F 主表
bash docs/experiments/baselines/scripts/run_sisl_failure.sh      # 失效构型
python3 docs/experiments/baselines/scripts/make_compare_csv.py   # 汇总 ssi_baselines_compare.csv
```

## 重要说明

- **未入库的大文件**：`real_parts/stp/`（931 MB）、`real_parts/data/`（852 MB）、`*.zip`、`baselines/sisl_build/` 均在 `.gitignore` 中；从 SimJEB 原始 STEP 用 `extract_simjeb.py` 可逐位重建 `data/`。
- **代码版本**：论文 §7 数字对应 commit `1ff0e1f`。历史 `real_parts/logs/` 由修复 AABB bug 之前的 harness 产生，其 disjoint 计数偏乐观，不应引用；见 `real_parts/results/bbox_bug_counts.txt`。
- **可复现性边界**：per-pair 计数可逐项复现；wall-clock 时间依赖并发负载，论文只报告同机同会话的测量值。
