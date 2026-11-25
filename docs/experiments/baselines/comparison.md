# 可执行基线汇总：SISL / OCCT IntPatch / 矩阵表示法 / 认证覆盖

本目录是论文 "NURBS 认证式统一求交（interval Krawczyk）" 的 **同输入定量对比基线**。
所有数字都来自本会话实际运行的程序，原始输出在 `results/`，脚本在 `scripts/`，源码在 `src/`。
**没有编造任何数字**；跑不通/未实现的项都在文中明确标注。

---

## 0. 基线一览（方法 / 保证类型 / 是否开源 / 同 fixture 表现 / 备注）

| 基线 | 具体实现 | 保证类型 | 同 fixture（A–F）表现 | 开源 | 运行时间（本机） | 备注 |
|---|---|---|---|---|---|---|
| **SISL 4.7** | `s1859`（曲面×曲面拓扑）+ `s1310`（追踪交线）；`s1858`（曲面×曲线）；`s1857`（曲线×曲线） | **无**：几何求解器，近似曲线/点；不证明存在/唯一/完备 | A/B/C/E/F 正确；**D 退化切触炸成 1628 条微曲线 + 346 点、1536 次追踪失败** | 是（AGPL/商业双许可） | 0.002–0.8 ms（A/B/C/E/F）；**2.9 s（D）** | 2D 重合只返回边界环（§4）；无曲面自交算子；输出拓扑随 epsge 变化 |
| **OCCT 8.0.1（IntPatch）** | `GeomAPI_IntSS`（S/S，tol=1e-6/1e-4）、`GeomAPI_IntCS`（S/C，无 tol 参数） | **无**：几何求解器，同上 | A/B/E/F 正确；**C 返回 0（公开 API 无孤立点）；D 返回 0 条曲线（整条交线漏掉）** | 是（LGPL） | 0.2–11 ms | `GeomAPI_IntSS` 的公开结果只有曲线，没有孤立交点 |
| **认证覆盖（本文）** | `cpp/tests/bench_ssi.cpp` 的严格区间排除 + 细分覆盖（本目录 `src/certified_cover_probe.cpp` 复刻同一 routine） | **有**：区间排除严格 ⇒ 覆盖盒带 ⊇ 解集（无遗漏）；参数盒 tol=1e-3 | A: 3070 覆盖盒 / 4052 排除；B: 0/1；C: 1/20；**D: 85344/40188（覆盖对角线）**；E: 951/1417 | 是（本仓库） | 6 ms（A）、**145 ms（D）**、2.8 ms（E） | 输出是盒带（不是显式曲线）；D 的代价比 A 高 ~24× |
| **矩阵表示法（minimal）** | TB–M 2009 curve/surface 特化 + Jia 2022 的 2×2/3×3 graph-patch 版本（`scripts/matrix_baseline.py`） | **无**（无存在性/重数证书） | curve×surface：1/3/5 根正确、二重根需去重、有理清分母产生多余特征值；2×2 S/S：A/B/C/E 正确，**D 每 t 返回 2 个重特征值 → 点数翻倍** | 是（本目录脚本） | 0.03–0.23 ms（CS）；0.7–26 ms（S/S 扫描） | **minimal reproduction，非 Jia 2022 完整实现**；只做"矩阵→特征值"这一步 |

> 认证覆盖的 `tol=1e-3` 是**参数域盒尺寸**，与 SISL 的 `epsge`、OCCT 的 `tol`
> （三维几何容差）**不是同一个量纲**，因此"覆盖盒数"不能与"交线数"直接比较，
> 只能与"是否漏解/是否有证书"比较。

---

## 1. 主表：同一组 fixture、同一输入

fixture（graph patch `S1=(s,t,f1)`、`S2=(s,t,f2)`，交集 = `f1−f2=0`）：

- **A** `f=s²−t²`：1 条对角线 `s=t`（原点相切、其余横截）
- **B** `f=0.3`：无交
- **C** `f=s+0.5t`：仅角点 `(0,0)` 切触（1 个孤立点）
- **D** `f=(s−t)²`：整条对角线 `s=t`，**退化切触**
- **E** 双三次×双三次 `u³+v³=0.1`：1 条内部交线
- **F** 有理四分之一圆柱 × 水平直线 `y=z=0.5`：1 个交点 `(0.866025,0.5,0.5)`

完整 CSV：`results/ssi_baselines_compare.csv`（SISL/OCCT/认证三方并列）。
SISL 原始：`results/sisl_ssi.csv`；OCCT 原始：`results/occt_ssi_raw.txt`；
认证原始：`results/certified_cover_probe.txt`。

| fixture | 真值 | SISL(1e-6) | OCCT(1e-6) | 认证覆盖(tol=1e-3) | 一致性判断 |
|---|---|---|---|---|---|
| A | 1 曲线 | 1 曲线，残差 1.1e-16，0.77 ms | 1 曲线，残差 5.4e-15，2.9 ms | 3070 覆盖盒，6 ms | **三方一致** |
| B | 无交 | 0，<0.01 ms | 0，0.55 ms | 全域排除（排除 1 盒） | **三方一致** |
| C | 1 点 | **1 点**（(0,0)），0.16 ms | **0（漏）**，0.3 ms | 1 个覆盖盒（角点），0.01 ms | SISL/认证对；OCCT 公开结果无法表达孤立点 |
| D | 1 条切触对角线 | **1628 微曲线 + 346 点，1536 次追踪失败，2.9 s** | **0 曲线（整条漏）**，0.5 ms | 85344 覆盖盒覆盖对角线，145 ms | **三者各不相同**：OCCT 漏，SISL 碎片化，认证覆盖无遗漏 |
| E | 1 曲线 | 1 曲线，残差 4.1e-7，0.76 ms | 1 曲线，残差 2.7e-8，7.2 ms | 951 覆盖盒，2.8 ms | **三方一致** |
| F | 1 点 | 1 点，残差 1.1e-16，0.06 ms | 1 点，残差 1.6e-16，4.6 ms | n/a（曲线-曲面，不在 bench_ssi 场景） | **两方一致** |

> OCCT 的 C/D 结果不是脚本 bug：OCCT 8.0.1 的 `GeomAPI_IntSS` 公开接口只有
> `NbLines()/Line()`，**没有孤立点输出**；D 的 0 条是 IntPatch 在退化切触上的
> 已知失效模式（仓库 `docs/occt-comparison/OCCT_COMPARISON.md` 亦有记录，但本表
> 数字是按同一 fixture、同一容差重新跑出来的，未抄旧值）。

---

## 2. 结论（建议论文这样写，三段）

**① 横截构型上三方一致。** A/B/E/F 上 SISL、OCCT、认证覆盖给出相同的交线/交点
数量与位置（残差 1e-7 以内），SISL 甚至更快（亚毫秒级）。因此论文不能写
"开源内核普遍算不出来"，也不能声称速度优势。

**② 退化/相切构型上基线失效，但失效方式不同。**

- OCCT `GeomAPI_IntSS`：D 返回 **0 条曲线（整条交线漏掉）**，C 返回空
  （公开 API 不能表达孤立切点）——**漏解（false negative）**。
- SISL `s1859`：D **不返回 0，而是把一条解析上光滑的对角线炸成 1628 条微曲线 +
  346 个孤立点、1536 次追踪失败，耗时 2.9 s**（epsge=1e-4 时降为 70 条+54 点、
  116 ms）。几何点都在对角线上（`max|u−v|≈2e-7`），但**拓扑不可用、输出随容差
  剧烈变化、没有任何"这条曲线存在且唯一"的证书**。
- 认证覆盖：D 用 **85344 个覆盖盒**覆盖整条对角线，排除 40188 盒，6 ms→145 ms
  （A 的 24×），**保证"覆盖盒带 ⊇ 解集"（无遗漏）**，代价是盒带不是显式曲线。

**③ 本文的差异是"证书与无遗漏"，不是"别人都算不出来"。** 论文应把卖点从
"OCCT 全面失败"改成："现有开源几何内核在退化切触构型上或漏解（OCCT）或给出
无证书、拓扑不稳定的碎片（SISL）；本文的贡献是把这一分支变成**可判定 + 可证明
无遗漏**的分类/覆盖结果"。这一定位比原来的表述更抗审，且与上面的实测数据一致。

---

## 3. 可比性边界（必须写进论文的 caveat）

1. **SISL / OCCT 是几何求解器**：给近似几何对象，**不对根的存在性/唯一性/完备性
   作证明**；两者都接受一个容差参数（SISL `epsge`，OCCT `tol`），但容差只影响
   数值精度，不构成证书。`GeomAPI_IntCS`（曲线-曲面）在 OCCT 8.0 中**没有容差
   参数**，只能用内部默认值。
2. **"交线数"不是可靠指标**：SISL 在 D 上交线数 70（epsge=1e-4）~1628（epsge=1e-6），
   OCCT 是 0，认证是"85344 覆盖盒"。三者的输出**语义不同**（曲线对象 / 无 /
   盒带），只能比较"是否漏解 / 是否有证书"。
3. **认证覆盖的 tol 是参数域盒尺寸**，与三维几何容差量纲不同，不能把
   "覆盖盒数"换算成"交线数"。
4. **SISL 没有 2D 重合输出类别**：两张曲面 2D 重合时 `s1859` 返回**重合区域的
   边界环**（全重合 4 条边界、部分重合 4 条边界、仅共边 1 条），`stat=0` 不报警。
   见 `results/sisl_overlap_probe.txt`。→ 论文的 coplanar-overlap 分类信息
   在 SISL 输出契约中不存在，这项对比只能做"能力有/无"，不能做定量。
5. **SISL 不提供曲面自交算子**：对它调 `s1859(S,S)` 返回参数方形 4 条边界；
   本文 self-intersection 与 SISL 的差异是**能力差异**。
6. **矩阵基线是 minimal reproduction**（见 `matrix_baseline.md`）：只复现了
   "矩阵→特征值"一步；完整 Jia 2022 的 moving-surface/Dixon 消元未实现。
   其耗时**不能**与 SISL/OCCT/认证比较，只能做机制（根数/重数/证书）对照。

---

## 4. 与最新文献的可比性（Yang & Jia, ACM TOG 2025）

文献：J. Yang, X. Jia, "Overlap Region Extraction of Two NURBS Surfaces",
*ACM Transactions on Graphics* 44 (2025), DOI
[10.1145/3763308](https://dl.acm.org/doi/10.1145/3763308)。

- **它做什么**：直接**提取两张 NURBS 曲面的重叠区域**（overlap region）——
  即把"二维重合"作为一等公民输出（本会话未能取得正文，仅据题目/文档定位；
  实施细节未读，不作定量陈述）。
- **本文的 coplanar-overlap 只对仿射片成立**：本文此前的"共面重叠"结论是在
  `S1=(s,t,f1)`、`S2=(s,t,f2)` 的仿射 graph-patch 上、通过 `f1−f2` 的恒等/区间
  判定得到的**分类结果**，不是一般 NURBS 曲面的重叠区域提取算法。
- **不可比的原因**：① 问题层次不同（分类 vs 区域提取；输出是标签 vs 参数化区域）；
  ② 输入范围不同（仿射片/同参数域 vs 一般 NURBS）；③ 本文没有重叠区域提取的
  输出契约（只有 `coplanar` 标签）。因此**不能**把本文与 Yang & Jia 2025 放在
  同一张"精度/时间"表里比较。
- **对论文的含义（必须降级/重新定位）**：coplanar-overlap 这一卖点应改为
  "本文能**判定**仿射片共面重合（并给出覆盖/排除）"，并在 related work 中明确
  引用 Yang & Jia 2025 作为"一般 NURBS 重叠区域提取"的 SOTA；否则会被审稿人
  指出"重叠区域提取你们没做、也没比"。
- 本目录的 SISL 重合探测（§3.4）恰好提供了旁证：**主流开源内核连"这是 2D 重合"
  都不输出**，只给边界环；这可以作为"为什么需要专门的重叠分类"的论据，但**不能**
  作为"我们比 Yang & Jia 更好"的证据。

---

## 5. 复现命令

```bash
# 0) 环境：WSL Ubuntu 22.04, g++ 11.4, cmake 3.22, pythonocc env:
#    /root/miniconda3/envs/occ/bin/python  (OCCT 8.0.1.2, numpy 2.4.6)
cd /mnt/d/qqdownload/nurbs-certified-framework/docs/experiments/baselines

# 1) 构建 SISL 静态库（源码复用 /mnt/d/sisl_dl/src/SISL-master, SISL 4.7）
bash scripts/build_sisl.sh

# 2) SISL A-F 主表 + 明细
bash scripts/build_run_sisl.sh

# 3) SISL 失效构型证据（多根 / 自交 / D 的容差行为）
bash scripts/run_sisl_failure.sh

# 4) SISL 2D 重合探测
g++ -O2 -std=c++17 -I/mnt/d/sisl_dl/src/SISL-master/include \
    src/sisl_overlap_probe.cpp sisl_build/libsisl.a -lm \
    -o sisl_build/sisl_overlap_probe
./sisl_build/sisl_overlap_probe | grep -v 'Zone equality' \
    | tee results/sisl_overlap_probe.txt

# 5) OCCT IntPatch 同 fixture（同容差 1e-6 / 1e-4）
/root/miniconda3/envs/occ/bin/python scripts/occt_ssi_fixtures.py \
    > results/occt_ssi_raw.txt

# 6) 认证覆盖（复刻 cpp/tests/bench_ssi.cpp 的区间覆盖）
g++ -O2 -std=c++17 src/certified_cover_probe.cpp -o sisl_build/certified_cover_probe
./sisl_build/certified_cover_probe | tee results/certified_cover_probe.txt

# 7) 矩阵表示法 minimal reproduction
/root/miniconda3/envs/occ/bin/python scripts/matrix_baseline.py \
    > results/matrix_baseline_raw.txt

# 8) 汇总 CSV
/root/miniconda3/envs/occ/bin/python scripts/make_compare_csv.py
```

## 6. 文件索引

| 路径 | 内容 |
|---|---|
| `sisl_baseline.md` | SISL 构建、API、A–F 表、失效构型、重合探测、可比性边界 |
| `matrix_baseline.md` | 矩阵表示法 minimal reproduction 与结果 |
| `comparison.md` | 本文件：汇总表、三段结论、可比性边界、Yang & Jia 2025 / Jia 2022 定位 |
| `src/sisl_ssi_harness.cpp` | SISL A–F harness（s1859+s1310+s1858） |
| `src/sisl_failure_probe.cpp` | 多根（s1857）/ 自交 / D 容差行为 |
| `src/sisl_overlap_probe.cpp` | 2D 重合探测（G/G2/H/H2/I/J） |
| `src/certified_cover_probe.cpp` | 认证区间覆盖（复刻 bench_ssi） |
| `scripts/matrix_baseline.py` | 矩阵表示法最小复现 |
| `scripts/occt_ssi_fixtures.py` | OCCT IntPatch 同 fixture |
| `scripts/build_sisl.sh`, `build_run_sisl.sh`, `run_sisl_failure.sh` | 构建/运行脚本 |
| `scripts/make_compare_csv.py` | 汇总 CSV（只解析已有结果，不产生新数字） |
| `results/ssi_baselines_compare.csv` | **主表（SISL / OCCT / 认证并列）** |
| `results/sisl_ssi.csv`, `sisl_ssi_detail.txt` | SISL A–F（epsge 1e-6 / 1e-4） |
| `results/sisl_failure_probe.txt` | 失效构型证据 |
| `results/sisl_overlap_probe.txt` | 2D 重合探测 |
| `results/occt_ssi_raw.txt` | OCCT 原始输出 |
| `results/certified_cover_probe.txt` | 认证覆盖原始输出 |
| `results/matrix_baseline_raw.txt` | 矩阵方法原始输出 |
| `results/matrix_curve_surface.csv`, `matrix_surface_surface.csv`, `matrix_certified_1d.csv` | 矩阵方法分表 |
| `results/build_sisl.log` | SISL 构建日志 |
