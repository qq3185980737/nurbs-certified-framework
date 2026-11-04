# 近年文献检索报告（2020–2026）
## 面向改投 CAGD 的 NURBS/Bézier 认证式统一求交论文

**检索日期**：2026-11-04
**检索执行**：Crossref REST API + OpenAlex API（结构化元数据，逐条可核验）
**目标**：回应 CAD 拒稿意见中"缺少挑战性几何算例、缺少与 state of the art 的定性/定量比较、未引用 Elber / Kim / Buse / Mourrain 等相关工作"的指控。

> **可核验性声明**：下文每一条的作者、标题、期刊、卷期页、DOI 均取自 Crossref/OpenAlex 返回的结构化记录，非凭记忆生成。凡元数据存在不确定性（例如 ACM TOG 的"提前分页"），已在条目内明确标注 **【待核】**。

---

## (A) 必补引用清单（12 篇，按优先级）

### A1. Jia, Li, Cheng (2022) — 矩阵表示求交的**直接后续**，加权最高

- **作者**：Xiaohong Jia, Kai Li, Jinsan Cheng
- **标题**：Computing the Intersection of Two Rational Surfaces Using Matrix Representations
- **期刊**：Computer-Aided Design, **150**, 103303 (2022-09)
- **DOI**：[10.1016/j.cad.2022.103303](https://doi.org/10.1016/j.cad.2022.103303)
- **定位句（英文草稿）**：*Jia et al. revive the matrix-representation approach of Thang, Busé and Mourrain for rational surface–surface intersection, which is precisely the algebraic baseline our interval certificate must be compared against rather than merely cited.*
- **为什么必补**：这是 Thang–Busé–Mourrain 2009 之后该学派路线**唯一**的实质性求交后续（Busé/Mourrain 本人 2020–2026 未再发表求交论文，见 (C)）。审稿人点名 Busé/Mourrain，若不引本文，等于承认不了解该路线的现状。

### A2. Park, Hong, Kim, Elber (2021) — Kim/Elber 自交线的**正牌后续**

- **作者**：Youngjin Park, Q Youn Hong, Myung-Soo Kim, Gershon Elber
- **标题**：Self-intersection computation for freeform surfaces based on a regional representation scheme for miter points
- **期刊**：Computer Aided Geometric Design, **86**, 101979 (2021-03)
- **DOI**：[10.1016/j.cagd.2021.101979](https://doi.org/10.1016/j.cagd.2021.101979)
- **定位句（英文草稿）**：*Park et al. continue the Elber–Kim self-intersection program by representing miter points regionally, and their geometric decomposition stands as the direct predecessor that our certified self-intersection test must either match in robustness or explicitly trade off against.*
- **为什么必补**：已引的 Pekerman–Elber–Kim 2008 是 18 年前的工作；本文是该组自交检测的最新代表（被引 11）。CAGD 编辑部会直接核对这一条。

### A3. Yang, Jia (2025) — 共面重叠的 state of the art

- **作者**：Jieyin Yang, Xiaohong Jia
- **标题**：Overlap Region Extraction of Two NURBS Surfaces
- **期刊**：ACM Transactions on Graphics, **44**, 1–20 (2025-12)
- **DOI**：[10.1145/3763308](https://doi.org/10.1145/3763308)
- **定位句（英文草稿）**：*Yang and Jia address exactly the coplanar-overlap configuration that degenerate cases of surface–surface intersection produce, so our overlap module must be benchmarked against their region extraction rather than presented as an untested add-on.*
- **为什么必补**：论文四大卖点之一是"共面重叠"。这是目前唯一专门处理 NURBS 曲面重叠域提取的 TOG 论文，不引则该卖点无法成立。

### A4. Burr, Hauenstein, Lee (2026) — **最接近的方法学先验工作（interval Krawczyk）**

- **作者**：Michael Burr, Jonathan Hauenstein, Kisun Lee
- **标题**：Certified Surface Approximations Using the Interval Krawczyk Test
- **出处**：ISSAC '26 — Proc. 2026 Int. Symp. on Symbolic and Algebraic Computation, pp. **77–86** (2026-07-12)
- **DOI**：[10.1145/3815436.3815444](https://doi.org/10.1145/3815436.3815444)
- **定位句（英文草稿）**：*Burr et al. already deploy the interval Krawczyk test to certify surface approximations, so our contribution cannot be the operator itself but must be its geometric unification across intersection, closest point, self-intersection and overlap with a no-omission cover.*
- **为什么必补**：**风险最高的一条**。同一套区间 Krawczyk 机器已用于曲面认证（2026 年 7 月，ISSAC）。若不引且不区分，等于把已知工具当新方法，会被直接判为 novelty 不足。

### A5. Li, Wu, Luo, Na, Cheng (2026) — 精度受控 SSI 的直接竞争者

- **作者**：Gangyi Li, Xi Wu, Zhongxuan Luo, Lei Na, Jin-San Cheng
- **标题**：A Precision Controlled Surface-Surface Intersection Algorithm for NURBS
- **期刊**：ACM Transactions on Graphics, **45**, 1–14 (2026-06-29)
- **DOI**：[10.1145/3806045](https://doi.org/10.1145/3806045)
- **勘误**：[10.1145/3839107](https://doi.org/10.1145/3839107)（TOG 45, 1–1, 2026-08-25）
- **定位句（英文草稿）**：*Li et al. achieve user-controlled precision in NURBS surface–surface intersection without interval certification, which isolates our claimed delta as existence-and-uniqueness proof rather than tighter tolerance handling.*
- **备注**：现有引用写作 "Li & Wu 2026 TOG"，作者列表不完整（实为五人），且未提勘误，需修正。

### A6. Yang, Jia (2026) — 小环与相切，区间证书最难的构型

- **作者**：Jieyin Yang, Xiaohong Jia
- **标题**：A Robust and Efficient Intersection Algorithm for NURBS Surfaces: Handling Small Loops and Tangent Intersections
- **期刊**：ACM Transactions on Graphics, **45**, 1–19 (2026-06-29)
- **DOI**：[10.1145/3807948](https://doi.org/10.1145/3807948)
- **勘误**：[10.1145/3839108](https://doi.org/10.1145/3839108)
- **定位句（英文草稿）**：*Small loops and tangential contacts are exactly where interval Krawczyk certification is most likely to fail to contract, so Yang and Jia's handling of these configurations defines the test set our method must survive.*
- **为什么必补**：直接给出"挑战性几何"清单（小环、相切），正好回应"examples on challenging geometries"的指控。

### A7. Wang, Jia, Yang, Wang, Bo, Liu (2025) — watertight 输出 = 无遗漏的竞争者

- **作者**：Yuqing Wang, Xiaohong Jia, Jieyin Yang, Bolun Wang, Pengbo Bo, Yang Liu
- **标题**：Improving the Watertightness of Parametric Surface/Surface Intersection
- **期刊**：Computer Graphics Forum (2025-12-27)
- **DOI**：[10.1111/cgf.70298](https://doi.org/10.1111/cgf.70298)
- **卷期页**：【待核】Crossref 记录尚未分配卷期页（early view）。
- **定位句（英文草稿）**：*Watertightness is the practical symptom of a cover with no missed branches, so this work is the closest empirical yardstick for our no-omission claim and must be compared on the same models.*
- **为什么必补**："无遗漏 cover"是本论文核心卖点，而 watertightness 是该卖点的可测代理指标。

### A8. Botbol, Busé, Chardin, Yildirim (2020) — Busé 近年**最近点**工作

- **作者**：Nicolás Botbol, Laurent Busé, Marc Chardin, Fatmanur Yildirim
- **标题**：Fibers of Multi-Graded Rational Maps and Orthogonal Projection onto Rational Surfaces
- **期刊**：SIAM Journal on Applied Algebra and Geometry, **4**(2), 322–353 (2020)
- **DOI**：[10.1137/19M1289522](https://doi.org/10.1137/19M1289522)
- **定位句（英文草稿）**：*Botbol et al. compute orthogonal projections onto rational surfaces through elimination matrices that depend linearly on the query point, giving a precomputed-algebra alternative to our per-query certified closest-point solve.*
- **为什么必补**：论文声称统一覆盖"最近点"。这是 Busé 组 2020 年后唯一直接对应最近点/正交投影的工作，且思路（预计算消去矩阵）与本方法（逐查询区间收缩）形成清晰对照。**同时它证明 Busé 近年未做求交，只做了投影**——这一点对 (C) 的答复很关键。

### A9. Breiding, Rose, Timme (2023) — 区间认证多项式系统零点的标准引文

- **作者**：Paul Breiding, Kemal Rose, Sascha Timme
- **标题**：Certifying Zeros of Polynomial Systems Using Interval Arithmetic
- **期刊**：ACM Transactions on Mathematical Software, **49**(1), 1–14 (2023-03-21)
- **DOI**：[10.1145/3580277](https://doi.org/10.1145/3580277)
- **定位句（英文草稿）**：*Breiding et al. establish interval-arithmetic certification of polynomial-system zeros as standard practice in numerical algebraic geometry, which is the correctness standard our certificates are held to.*
- **为什么必补**：CAGD 的认证类审稿人会要求引用符号计算社区公认的区间认证基线；缺此条会被认为"certificate"用词不严谨。

### A10. Elber (2022) — "保守/无遗漏"分析的 Elber 先例

- **作者**：Gershon Elber
- **标题**：Detection and computation of conservative kernels of models consisting of freeform curves and surfaces, using inequality constraints
- **期刊**：Computer Aided Geometric Design, **94**, 102075 (2022)
- **DOI**：[10.1016/j.cagd.2022.102075](https://doi.org/10.1016/j.cagd.2022.102075)
- **定位句（英文草稿）**：*Elber derives conservative (guaranteed no-false-negative) geometric answers from inequality constraints on freeform primitives, which is the same epistemic goal our interval cover pursues and the natural precedent to cite for it.*
- **为什么必补**：审稿人点名 Elber；且这是 Elber 近年唯一以"保守性/保证"为核心的自由曲面全局分析工作，正好支撑"无遗漏"论述。

### A11. Liu, Cheng, Liu, Guo, Du, Fu (2026) — 相邻曲面片的求交**检测**

- **作者**：Junyi Liu, Yuan-Yuan Cheng, Shibo Liu, Jia-Peng Guo, Wei Du, Xiao-Ming Fu
- **标题**：Efficient intersection detection of adjacent polynomial parametric surfaces
- **期刊**：Computer Aided Geometric Design, **127**, 102553 (2026-06)
- **DOI**：[10.1016/j.cagd.2026.102553](https://doi.org/10.1016/j.cagd.2026.102553)
- **定位句（英文草稿）**：*Liu et al. target the decision problem for adjacent patches, where shared boundaries make naive interval exclusion fail, and our cover must be shown to remain sound in that same adjacency regime.*
- **为什么必补**：CAGD 同期刊、同年、同子领域，编辑部几乎必然知道；且"相邻片共享边界"是区间排除法的经典失效点，属于必须回应的算例族。

### A12. Duff, Lee (2024) — Krawczyk 认证在符号计算社区的落地

- **作者**：Timothy Duff, Kisun Lee
- **标题**：Certified homotopy tracking using the Krawczyk method
- **出处**：ISSAC '24 — Proc. 2024 Int. Symp. on Symbolic and Algebraic Computation, pp. **274–282** (2024-07-16)
- **DOI**：[10.1145/3666000.3669699](https://doi.org/10.1145/3666000.3669699)
- **定位句（英文草稿）**：*Duff and Lee show the Krawczyk operator certifying entire solution paths rather than isolated roots, which is the path-level guarantee our tracing stage implicitly relies on and should cite explicitly.*
- **为什么必补**：把"certificate"从单点推广到整条路径，正好对应本论文"曲线求交 + 追踪"的路径级保证。与 A4 一起构成认证线的完整引文对。

---

## (B) 可选扩展清单（10 篇）

| # | 文献 | 出处 / DOI | 定位（中文一句话 + 英文草稿） |
|---|---|---|---|
| B1 | Xiao, Chen — *Intersection of two NURBS surfaces with tolerance control* | Computer-Aided Design **203**, 104196；[10.1016/j.cad.2026.104196](https://doi.org/10.1016/j.cad.2026.104196)。**注**：Crossref issued 2027-02（2026 年在线），引用时年份需按最终卷期确认 **【待核】** | 同为"容差受控"路线，但不给存在唯一性证明。*Tolerance-controlled but certificate-free — a direct foil for our guarantee claim.* |
| B2 | Chu, Li, Jia, Yang, Kang — *Computing the intersection of two ellipsoids based on a fast algebraic topology determination strategy* | CAGD **119**, 102442；[10.1016/j.cagd.2025.102442](https://doi.org/10.1016/j.cagd.2025.102442) | 代数拓扑判定 + 快速分类，可作为"拓扑保证"对照算例来源。*Fast algebraic topology determination for a quadric subcase of SSI.* |
| B3 | Zhang, Wu, Cheng, Ding — *Computing the intersection between a rational parametric curve and a rational parametric surface* | CAGD **109**, 102275；[10.1016/j.cagd.2024.102275](https://doi.org/10.1016/j.cagd.2024.102275) | 曲线–曲面求交的最新代数方法，与本论文"曲线求交"分支直接可比。*Curve–surface intersection via algebraic elimination, a benchmark for our curve branch.* |
| B4 | Shao, Chen — *Topological classification of the intersection curves of two quadrics using a set of discriminants* | CAGD **107**, 102244；[10.1016/j.cagd.2023.102244](https://doi.org/10.1016/j.cagd.2023.102244) | 判别式驱动的拓扑分类，提供 ground-truth 拓扑标签用于验证 no-omission。*Discriminant-based topology labels usable as ground truth.* |
| B5 | Spainhour, Weiss — *Robust Containment Queries over Collections of Trimmed NURBS Surfaces via Generalized Winding Numbers* | TOG **45**, 1–21；[10.1145/3797957](https://doi.org/10.1145/3797957) | trimmed NURBS 上的鲁棒查询，覆盖"trimmed surface"这一审稿人常问场景。*Robust queries on trimmed NURBS, the trimming regime our cover must also address.* |
| B6 | Kim 组 — *Precise Hausdorff distance computation for freeform surfaces based on computations with osculating toroidal patches* | CAGD **86**, 101967；[10.1016/j.cagd.2021.101967](https://doi.org/10.1016/j.cagd.2021.101967) | Kim 组近年距离/最近点方向代表，与 Park 2020 的 toroidal patch BVH 同源。*Osculating-toroidal-patch machinery extended to precise Hausdorff distance.* |
| B7 | Amend, Goertzen — *A framework for symmetric self-intersecting surfaces* | Graphical Models **145**, 101329；[10.1016/j.gmod.2026.101329](https://doi.org/10.1016/j.gmod.2026.101329) | 覆盖 Graphical Models 期刊，且自带已知自交结构的构造性算例，可直接用作测试模型。*Constructive families with known self-intersection, useful as generated test cases.* |
| B8 | Dou, Cheng, Wen — *A Maple Implementation for Deterministically Certifying Isolated Simple Zeros of Over-Determined Polynomial Systems with Interval Arithmetic* | AppliedMath **6**(1), 15；[10.3390/appliedmath6010015](https://doi.org/10.3390/appliedmath6010015) | 超定系统的区间认证——SSI 恰是超定（3 方程 / 4 未知）情形，方法学高度对口。*Interval certification for over-determined systems, structurally the SSI setting.* |
| B9 | He, Wang, Zhao — *Self-intersection detection algorithm of sweep surfaces based on geometric features of spine curves* | Computer-Aided Design **192**, 104024；[10.1016/j.cad.2025.104024](https://doi.org/10.1016/j.cad.2025.104024) | 扫掠面自交的专用几何判据，代表"利用结构先验"的另一类思路。*Structure-exploiting self-intersection detection for sweep surfaces.* |
| B10 | Mourrain 等 — *A certified iterative method for isolated singular roots* | J. Symbolic Computation **115**, 223–247；[10.1016/j.jsc.2022.08.006](https://doi.org/10.1016/j.jsc.2022.08.006) | Mourrain 近年的"认证"工作（奇异根），可支撑"certified"术语并为奇异/相切构型提供理论依托。*Certified iteration for singular roots, relevant to tangential degeneracies.* |

---

## (C) 与现有引用的重复 / 冲突检查

### C1. 已引文献核对结果（逐条）

| 现有引用 | 核验到的真实记录 | 结论 |
|---|---|---|
| Cheng et al. 2023 TOG | Cheng, Zhang, Xiao, Li, *Topology driven approximation to rational surface-surface intersection via interval algebraic topology analysis*, TOG **42**, 1–16, [10.1145/3592452](https://doi.org/10.1145/3592452) | ✅ 存在。**但注意**：该文已用 interval 方法做拓扑驱动求交，是本论文"interval"卖点最强的已引先例，正文必须明确区分（详见 C2） |
| Yang et al. 2023 TOG | Yang, Jia, Yan, *Topology Guaranteed B-Spline Surface/Surface Intersection*, TOG **42**, 1–16, [10.1145/3618349](https://doi.org/10.1145/3618349) | ✅ 存在。与上条同为 TOG 42、pp. 1–16，**极易在参考文献表中混淆**，需逐字核对 |
| Li, Jia, Chen 2025 TOG | Li, Jia, Chen, *Fast Determination and Computation of Self-intersections for NURBS Surfaces*, TOG **44**, 1–16, [10.1145/3727620](https://doi.org/10.1145/3727620) | ✅ 存在，作者与卷期正确 |
| Li & Wu 2026 TOG (precision-controlled SSI) | **实为 Li, Wu, Luo, Na, Cheng**，*A Precision Controlled Surface-Surface Intersection Algorithm for NURBS*, TOG **45**, 1–14, [10.1145/3806045](https://doi.org/10.1145/3806045) | ⚠️ **作者不全**（漏 Luo, Na, Cheng），且**存在勘误** [10.1145/3839107](https://doi.org/10.1145/3839107)。需修正 |
| Zhang et al. 2025 CAGD | Zhang, Cheng, Liu, Zhang, CAGD **118**, 102432, [10.1016/j.cagd.2025.102432](https://doi.org/10.1016/j.cagd.2025.102432) | ✅ 存在 |
| Li, Yang, Jia 2026 CAD (SSI survey) | Li, Yang, Jia, *Advances and challenges in surface–surface intersection computation — An overview*, CAD **193**, 104039, [10.1016/j.cad.2026.104039](https://doi.org/10.1016/j.cad.2026.104039) | ✅ 存在。**这篇 survey 是回应"未引用 state of the art"的最佳武器**：应在 Introduction 用它逐条对照自查 |
| Park, Son, Kim, Elber 2020 CAD | CAD **127**, 102866, [10.1016/j.cad.2020.102866](https://doi.org/10.1016/j.cad.2020.102866) | ✅ 存在 |
| Mourrain & Pavone 2009；Thang, Busé, Mourrain 2009；Elber & Cohen 1993；Pekerman et al. 2008；Elber, Grandine, Kim 2009；Grandine & Klein 1997；Hu et al. 1997；Sederberg 1988/1989；Ko & Patrikalakis 2008 | 均为 2009 年及更早，不属于本次 2020–2026 检索范围 | ✅ 未发现重复条目；但**其后续工作缺失**见 C3 |

### C2. 潜在"新颖性冲突"（必须在正文主动化解）

1. **Cheng et al. 2023 TOG（已引）已经用 interval 代数拓扑做 SSI。** 若论文把"interval"本身当作卖点，会被认为增量不足。建议把卖点重新表述为：**统一性（一个算子覆盖求交/最近点/自交/重叠）+ 存在唯一性证书（而非仅拓扑正确）+ 无遗漏 cover 的可证性**，并把 Cheng 2023 明确列为"仅拓扑、无唯一性证书"的对照。
2. **Burr, Hauenstein, Lee 2026 ISSAC（新增 A4）已用 interval Krawczyk 认证曲面。** 这是同一算子的先验使用，必须引用并说明差异（他们认证的是曲面逼近/同伦，不涉及求交与重叠的几何判定）。
3. **A1（Jia 2022 CAD）与 A8（Botbol 2020 SIAA）分别占据"求交"和"最近点"的代数路线制高点**，本论文的对照实验若不含矩阵表示法基线，量化比较会被认为不完整。

### C3. 审稿人点名作者的 2020–2026 现状（**已核实的负面结论，可直接写进 response letter**）

- **Laurent Busé**：2020–2026 共 31 条记录（OpenAlex A5088510335），**无一篇以曲线/曲面求交为主题**。相关者仅 A8（正交投影/最近点，2020）与 *Solving bihomogeneous polynomial systems with a zero-dimensional projection*（ISSAC 2025, pp. 206–214, [10.1145/3747199.3747563](https://doi.org/10.1145/3747199.3747563)）。其余转向 Sylvester form、双有理映射、张量分解。
- **Bernard Mourrain**：2020–2026 共 45 条记录（OpenAlex A5003683965），**同样无求交论文**；转向矩/SOS 多项式优化、几何连续样条、认证数值（B10）。
- **Myung-Soo Kim 组（KAIST/SNU）**：2020–2026 求交/自交相关为 Park–Son–Kim–Elber 2020（已引）、A2（Park–Hong–Kim–Elber 2021）、B6（Hausdorff 距离 2021）；2022 年后该组转向 lattice / 微结构 / 增材制造，**求交线已基本停止**。
- **Gershon Elber**：2020–2026 共 33 条记录（OpenAlex A5055421841），求交/自交相关为 Park 2020（已引）、A2、A10；其余为 V-reps、微结构、晶格。
- **结论**：**不存在被漏掉的 Busé/Mourrain/Kim/Elber 2020–2026 求交论文**。真正的 state of the art 已转移到 Jia/Cheng/Yang 组（KAUST/CAS 系）与 Fu Xiao-Ming 组。这一点既能正面回应审稿人，也说明为何"必补清单"里 A1/A3/A5/A6/A7/A11 权重更高。建议 response letter 明确写出这句。

### C4. 引用格式风险

- ACM TOG 2025–2026 各条页码（1–14 / 1–16 / 1–19 / 1–20 / 1–21）为 ACM **advance pagination**，非最终期号页码。定稿前应重新核对 **【待核】**。
- A7（CGF 2025）尚无卷期页，属 early view **【待核】**。
- B1（CAD 104196）在线 2026、卷期标 2027-02，年份标注需谨慎 **【待核】**。

---

## (D) 检索方法与检索词（可复现）

### D1. 为何不用 Google Scholar / dblp

- `dblp.org` 在本环境返回 HTTP 202 + Anubis 反爬挑战页（"Making sure you're not a bot!"），无法取得结构化记录。
- `api.semanticscholar.org` 返回 HTTP 429（限流）。
- 因此改用 **Crossref REST API**（权威 DOI 元数据：作者、container-title、volume、issue、page、issued）与 **OpenAlex API**（作者消歧 ID、引用数、venue）。二者均为机器可读、可复现、可逐条核验。

### D2. Crossref 查询（基础 URL：`https://api.crossref.org/works`）

```
# 1. 标题精确定位
?query.bibliographic=Topology+guaranteed+and+error+controlled+curve+tracing+parametric+surface-surface+intersection
 &rows=5&select=DOI,title,author,container-title,volume,page,issued

# 2. 期刊定向扫描（CAGD）
?query.container-title=Computer+Aided+Geometric+Design&query.bibliographic=intersection
 &rows=25&filter=from-pub-date:2020-01-01,type:journal-article
 &select=DOI,title,container-title,volume,page,issued,author

# 3. 期刊定向扫描（CAD）
?query.container-title=Computer-Aided+Design&query.bibliographic=intersection
 &rows=25&filter=from-pub-date:2020-01-01,type:journal-article

# 4. 期刊定向扫描（TOG）
?query.container-title=ACM+Transactions+on+Graphics&query.bibliographic=intersection+NURBS+surface
 &rows=20&filter=from-pub-date:2022-01-01,type:journal-article

# 5. 期刊定向扫描（Graphical Models）
?query.container-title=Graphical+Models&query.bibliographic=intersection+self-intersection+curve
 &rows=20&filter=from-pub-date:2020-01-01,type:journal-article

# 6. 主题检索（自交）
?query.bibliographic=self-intersection+detection+NURBS+surface
 &rows=15&filter=from-pub-date:2020-01-01,type:journal-article

# 7. 作者 + 主题（Kim 组）
?query.author=Myung-Soo+Kim&query.bibliographic=surface+intersection
 &rows=20&filter=from-pub-date:2020-01-01

# 8. 认证数值基线
?query.bibliographic=Certifying+zeros+of+polynomial+systems+using+interval+arithmetic&rows=3

# 9. DOI 直查（逐条核验 A4 / A12 / A8）
/10.1145/3815436.3815444
/10.1145/3666000.3669699
/10.1137/19m1289522
```

> 注：Crossref 的 `query.author=Laurent+Buse` 对 "Busé" 的模糊匹配**严重失效**（返回大量以 "Buse" 为名的无关作者，如航空燃料、心理治疗专著）。作者维度必须改用 OpenAlex 的作者 ID，见 D3。

### D3. OpenAlex 查询（基础 URL：`https://api.openalex.org`）

```
# 作者消歧（取得稳定 ID）
/authors?search=Laurent%20Bus%C3%A9&select=id,display_name,works_count
   -> A5088510335
/authors?search=Bernard%20Mourrain&select=id,display_name,works_count
   -> A5003683965
   Myung-Soo Kim -> A5100718911 ;  Gershon Elber -> A5055421841

# 作者全量近年产出（本报告 C3 的依据）
/works?filter=authorships.author.id:A5088510335,from_publication_date:2020-01-01
 &per_page=50&select=doi,title,publication_year,biblio,type,cited_by_count&sort=publication_date:desc
# 同上替换 author.id 为 A5003683965 / A5100718911 / A5055421841

# 主题 + fulltext 检索
/works?search=surface%20intersection&filter=from_publication_date:2020-01-01&per_page=25
/works?filter=title_and_abstract.search:Krawczyk,from_publication_date:2018-01-01&per_page=30
/works?filter=title_and_abstract.search:self-intersection%20surface%20freeform,from_publication_date:2019-01-01
```

### D4. 主题词清单（实际使用，含命中/未命中）

| 检索词 | 结果 |
|---|---|
| `surface-surface intersection`, `SSI NURBS` | ✅ 高产出（TOG/CAGD/CAD 主战场） |
| `self-intersection NURBS surface` | ✅ 产出 A2、B9、B7 |
| `interval arithmetic geometric` / `interval Krawczyk` | ✅ **关键命中 A4**；泛化的 "interval arithmetic intersection curve" 噪声极大（混入椭圆曲线密码学） |
| `certified / validated numerics geometry` | ✅ 产出 A9、A12、B8、B10 |
| `Krawczyk`（直接查） | ✅ 705 条中筛出 A4/A12 及 2019 年若干应用；同名人物噪声多 |
| `closest point NURBS` / `orthogonal projection rational surface` | ⚠️ "closest point" 直接检索几乎全被 IGA 接触力学占据；改以作者维度（Busé/Kim）才命中 A8、B6 |
| `matrix representation intersection` | ⚠️ Crossref 检索弱；改以期刊扫描命中 A1 |
| `subdivision solver polynomial system` | ❌ 0 条（OpenAlex AND 语义过严）；相关线索改由引用网络获得（Mourrain 2009 后续已转移，见 C3） |
| `topology guarantee SSI` | ✅ 产出 Cheng 2023、Yang 2023、A6、B2、B4 |
| `trimmed surface intersection` | ⚠️ 弱；改由 TOG 扫描命中 B5、Antolin/Buffa/Cirillo (CAD 156:103448) |

### D5. 覆盖的期刊/会议

Computer-Aided Design; Computer Aided Geometric Design; ACM TOG; Computer Graphics Forum; Graphical Models; CAD&A; J. Symbolic Computation; ACM TOMS; SIAM J. Appl. Algebra Geom.; AppliedMath; ISSAC 2024/2025/2026; CSIAM Trans. Appl. Math.

> **未覆盖（诚实声明）**：JCAM、JSC 之外的数值分析期刊仅抽样；CAD&A 未做系统性 2020–2026 全量扫描（其 Crossref 元数据质量差，混入大量非几何论文）；中文期刊（如《计算机辅助设计与图形学学报》）仅偶然命中，未系统检索。ISSAC/SC 会议论文集仅通过 ACM DOI 前缀与 OpenAlex 索引间接覆盖，可能仍有遗漏 **【待核】**。

---

## 附：一句话总结

审稿意见说"未引用 Elber / Kim / Busé / Mourrain 的 state of the art"——检索证明这四位在 2020–2026 **确实没有新的求交论文**（C3，已核实），所以正确做法不是硬找他们的新文，而是：(i) 补 A2/A10/A8 交代这四条线的**最新落点**；(ii) 把真正的 state of the art 换成 A1/A3/A5/A6/A7/A11 并做定量对比；(iii) 用 A4/A9/A12 把"认证"锚定到符号计算社区标准。风险最高的是 **A4（interval Krawczyk 已被用于曲面认证）** 与 **已引的 Cheng 2023 TOG（interval 已用于 SSI）**，这两条决定新颖性论证能否成立。
