# 四步加强路线评估与修订（CAD 拒稿 → CAGD）

- 评估日期：2026-11-04
- 评估对象：`C:\Users\15632\Doubao\chats\2026-10-13\new-chat\CAD整理\`（rev19 清理版 + 拒稿通知 + OCCT 对比 + 文献规划）
- 评估对象论文：`P1_CAD_rev19_清理版.docx`（CAGD 投稿版，20 页 / 8 节 / 11 表 / 5 图 / 31 refs）
- 评估对象路线：用户提出的四步（① 真实零件批量测试；② Elber/Kim/Busé/Mourrain 结构性对比；③ 近年文献补引；④ 更多开源实现对比）
- 本文件同时是「每一步留文件」的第 0 号产物（S0），后续步骤产物见 §5 表格与 `docs/experiments/` 下各目录。

---

## 0. 结论速览

**方向正确，但需要修改；不需要推翻。**

四步恰好一一对应主编的三条批评（见 §1），说明路线没有跑偏。问题出在**执行形态**：

| # | 原步骤 | 判断 | 必须的修改 |
|---|---|---|---|
| 1 | 真实零件提取 + 批量测试 | 方向对，**当前实现不能写进论文** | 改成"有界协议 + OCCT 交叉校验 + 可扩展性曲线"三件套；不要报"全零件通过率" |
| 2 | Elber/Kim/Busé/Mourrain 结构性对比 | **不充分**：rev19 的 Table 9 已经是结构性对比 | 必须补 1–2 个**可运行复现**（Busé–Mourrain 矩阵法、Pekerman–Elber–Kim 细分自交），否则主编原话"benchmarked against"仍未回答 |
| 3 | 近年求交/自交/最近点文献补引 | 必要，但只补引会被判"堆引用" | 每篇必须带一句"与本方法的差异"，并落进 Table 9 / 对比节 |
| 4 | 更多开源实现对比 | 方向对，**现实可选对象比想象少** | 把"更多开源实现"改为"更多**可运行求解范式**"：SISL 曲面例程 + 代数求解器 + OCCT 两条 API 路径 |

**路线之外还缺三件事**（这是本次评估最重要的结论，见 §4）：
(a) 性能/可扩展性叙事 —— 真实零件会直接暴露认证 cover 的成本，不处理会反噬；
(b) 证书严格性的一致口径 —— "certified" 在 5 个算子上强度不同，需要一张保证等级表；
(c) 论文定位 —— 5 个算子塞进 20 页，CAGD 审稿人会判"宽而不深"，主线应收到"认证 SSI"。

**做完文献核查后的追加结论（最重要，见 §8）**：四步全部是"补证据"，**没有一步是"补新颖性"**。而真正的风险恰恰在后者——
已引的 Cheng et al. TOG 2023 已用 interval 做 SSI，Burr et al. ISSAC 2026 已用 interval Krawczyk 认证曲面，
Yang & Jia TOG 2025 已做 NURBS 曲面重叠域提取（直接吃掉了本文的"共面重叠"卖点）。
因此修订版路线必须加入 **S0.5「新颖性重定位」**：把贡献从"5 个算子并列"收成"认证自交 + 统一证书机器"，
并把共面重叠降级。否则四步做完仍会因 incremental 被拒。

---

## 1. 主编意见的可检验拆解

拒稿信（`CAD拒稿通知_20261026.jpg`，Manuscript CADJ-26-00624）原文：

> "…your manuscript does not meet the required quality standards … In particular, **the paper is way behind the journal level in terms of examples on challenging geometries and comparison against the state of the art (both qualitative and quantitative)**. This is an old topic, many relevant state of the art papers were formerly published in CAD, but are **not referenced nor the proposed work is benchmarked against them**, see works of Elber, Kim, Buse, Mourrain, etc."

拆成三条可检验要求：

- **R1 例子**：challenging geometries 的**数量与难度**都要够。注意"challenging"不等于"随机大批量"——主编要的是相切、退化、小环、多根这类公认失效构型，且要来自**真实 CAD 数据**或**公认基准**，不能全是自造解析例。
- **R2 对比**：qualitative **and** quantitative。定性表只能算一半。
- **R3 点名**：Elber（IRIT 学派）、Kim M.-S.（KAIST，自交 + BVH）、Busé / Mourrain（代数-矩阵表示、细分求解）。关键词是 **benchmarked against**，即要在**同一输入上比结果**，而不是在 related work 里提一句。

一个容易误判的点：主编写的是 "not referenced"，但 rev19 其实已经引了 [20]–[31]（含 Elber、Kim、Busé、Mourrain 共 7 篇，见 `01_文章正文/...批注版` 的 85 条批注，其中 12 条明确是"回应主编"）。这说明 rev19 是**收到拒稿后已经改过一轮的稿子**。因此再改时真正缺的不是"引用"，而是 **R2 的 benchmark 部分**——把 step 2 理解成"补结构性对比"会重复 rev19 已做过的事。

---

## 2. rev19 已经做了什么、还缺什么

### 2.1 已做的回应（都是加分项，应保留）

| 位置 | 内容 | 评价 |
|---|---|---|
| §1.1 重写 | 补 CAD/CAGD 经典与近期 12 篇 [20]–[31]，按 tracing / subdivision-algebraic / interval / BVH / topology 五条线组织 | 已达标（R3 的 "referenced"） |
| 贡献 (5) | 把"困难几何基准 + 定性/定量双层对比"写进 contribution | 好，但依赖 §7.7/7.8 的真材料 |
| §7.7 + Table 8 | 9 类困难构型（多根、曲线-曲面、有理权重、相切、退化相切、角点接触、折叠自交、共面重叠、秩亏） | 已达 R1 的"难度"，但**全是自造解析例**，无真实数据 |
| §7.8 + Table 9/10/11 | Table 9 定性 16 类方法；Table 10 定量对 OCCT/SISL；Table 11 逐例 | Table 9 已覆盖 Elber/Kim/Busé/Mourrain —— **step 2 的"结构性对比"部分已经存在** |
| §7.9 | 端到端 + ACIS 曲率分类对齐 | 保持 |

### 2.2 缺口（决定四步怎么改）

| 缺口 | 证据 | 对应步骤 |
|---|---|---|
| 例子全部合成，无真实 CAD 零件 | §7.1–7.4 全是随机生成；§7.7 全解析 | Step 1 |
| 定量对比只比"数根"，且只在曲线-曲线一对上做过 SISL 对比 | Table 10 三类基线里，SISL 只有 s1857（曲线），OCCT 是逐例探针 | Step 1/2/4 |
| 点名的方法没有一个是"同输入同容差"跑出来的 | Table 9 无代码方法只给定性；正文 §7.8 自己承认"无法做 runtime 对比" | Step 2 |
| 没有性能/可扩展性证据 | 全文最快数字是曲线单对 0.028–1.203 ms；无曲面/零件级吞吐 | 新增 |
| 证书强度口径不一 | §2.3/2.4 是定理级；§6.2 点包含是证据级；§3.5 条件数阈值是启发式；摘要统一写 "certified" | 新增 |

### 2.3 rev19 数字自检（抽样复核，均成立）

- 6,338 − 6,265 = **73** ✓；1.28 = 6,338/4,950 ✓；4,675/9,000 = 0.519 ≈ 0.52 ✓
- Table 1 逐行合计 4,675 = Total ✓
- 500/500 自交 + 500/500 规则片；500/500 重叠 + 500/500 分离 ✓（与 `experiment_*.md` 一致）
- 结论：**rev19 的数据没有发现硬错误**，问题在"数据覆盖不到主编要求"，不在"数据算错"。

---

## 3. 逐条评估四步

### 第一步：真实零件批量测试 —— 方向对，当前实现会自伤

**已发生的事实**（`docs/experiments/real_parts/`，2026-11-04 11:00–14:22 的运行）：

- 数据集：SimJEB（MIT，382 个 STEP 支架，Open Data Commons Attribution，`stp/README_cad.txt` 有许可），已全部转成 Bézier patch（`data/*.txt`，NSURF/SURF/ENDSURF 纯文本，851 MB）。
- 规模分布（`manifest.csv`，381 件）：faces 中位数 **282**（56–6,146）；patches 中位数 **1,067**（85–74,647）；faces>500 的有 76 件，patches>5,000 的有 53 件。
- 跑了 19 个（`selected.txt`，挑的是面数最少的档），结果：

| 零件 | faces | patches | 结果 | 耗时 |
|---|---|---|---|---|
| 228 | 93 | 103 | rc=0 | 202 s |
| 387 | 82 | 110 | rc=0 | 82 s |
| 411 | 76 | **7,480** | rc=0（3,314,150 个 patch 对，其中 try 18,759） | **542 s** |
| 107 | 132 | 804 | rc=0 | 298 s |
| 10 / 59 / 69 / 70 / 124 / 200 / 205 / 214 / 248 / 506 / 534 / 538 / 546 / 570 | 56–499 | — | **rc=124 超时** | ≥300 s |
| 614 | 149 | 7,946 | rc=124 超时 | ≥600 s |

- 本次复核再加一测：56 面 / 197 patch 的零件 248，`timeout 120` 仍**跑不完**（无任何输出即被 kill）。

**为什么慢（根因，不是 bug）**：
1. `cover(A,B,tau)` 把**整条 1D 交线**用参数盒覆盖到宽度 τ=0.02，盒子数 ∝ (1/τ)^d，且相邻面共享一条完整边界（交集贯穿参数域），必须细分到底；
2. 真实 B-rep 的绝大多数"相交面对"是**相邻面**（共享边），这是**退化接触**而不是横截求交，恰好是最贵的那类；
3. 每个 face 被 `GeomConvert_BSplineSurfaceToBezierSurface` 切成多片（最大 74,647 patch/件），patch 对爆炸。

**结论**：如果照现状"跑完 382 件、报通过率"，审稿人看到的是"每件 1–10 分钟、还有 4%–8% undecided 的认证器"，等于**自证不可用于工程**。这条路不是不能走，是**不能这么走**。

**改法（3 条，缺一不可）**：

- **P1 有界协议**：不做"全零件全 face-pair 认证"，改为**分层抽样 + 每对预算**：每件抽 ≤N 个 bbox 相交 patch 对（另抽 bbox 不相交对作阴性对照），记录每个 patch 对的分类（disjoint / 1D / 2D / coplanar / undecided）与耗时；给出"耗时 vs faces/patches"的可扩展性曲线。
- **P2 真值来源**：真实零件上唯一可信的真值是另一个内核。用 OCCT `GeomAPI_IntSS` 对同一 patch 对求交，构造 2×2 混淆矩阵（认证 1D vs OCCT 有交线；disjoint vs OCCT 无），报告一致率 + 分歧案例清单。**这才是"证明了什么"**。
- **P3 覆盖检验**：对 OCCT 给出的交线采样点，检查是否落在认证 cover 盒内（cover 是"无遗漏"包装器，这是它唯一能直接验的性质）。这条比"数根一致"更有说服力，且与 §6.3 的 Proposition 1 呼应。

**验收标准**：`results.csv`（≥12 件）+ `occt_crosscheck.csv`（≥600 对）+ 一致率数字 + 1 条可扩展性曲线 + 明确写"认证 cover 适合局部认证，不适合全 B-rep 扫描，工业级需要 BVH + trimming + 局部化"。

---

### 第二步：Elber/Kim/Busé/Mourrain 结构性对比 —— 不充分

**为什么不够**：
- rev19 的 Table 9 已经是 16 行结构性对比（guarantee / 退化处理 / 是否开源），step 2 再做一遍是**重复劳动**；
- 主编原话是 "benchmarked against them"，**benchmark 是动词**；
- rev19 §7.8 自己写了两条 caveat："无法获得 [16–19,21,22,24–30] 的可运行实现"、"不做 runtime 对比"。审稿人会直接把这两条当作"没回答主编"。

**可行的分层策略**（按证据强度从高到低）：

| 层级 | 做法 | 覆盖谁 | 可行性 |
|---|---|---|---|
| L1 同输入可运行复现 | ① Busé–Mourrain 矩阵表示法（curve/surface，矩阵特征值给交点）；② Pekerman–Elber–Kim 型细分自交检测；③ Grandine–Klein 型追踪（可用 OCCT 的 IntPatch 近似代表） | Busé/Mourrain、Elber/Kim | 高（① 2–3 天，② 2 天） |
| L2 复算其论文公开算例 | 用原文给出的测试例与数字逐条复算，列表"原文结果 / 我们的复算 / 差异" | Elber & Cohen（IRIT，无源码） | 中 |
| L3 定性表 | rev19 Table 9 已有 | 全部 | 已做 |

**修改建议**：step 2 与 step 4 **合并**成"可执行基线"一条工作线，产出 `docs/experiments/baselines/`，论文里替换/扩写 Table 10，并新增一张"方法 × 保证类型 × 同输入结果 × 运行时间"的表；Table 9 保留但降为背景。

---

### 第三步：近年文献补引 —— 必要，但要"引用 + 区分 + 落表"

- 检索范围（建议）：2023–2026，Computer-Aided Design / CAGD / ACM TOG / Graphical Models / CAD&A / SGP-CAGD 特刊 / ISSAC-SNC（certified 方向）。
- 重点补：**自交**（Li–Jia–Chen TOG 2025 已引，需查其后续与对比数据）、**closest point / distance query**（目前几乎没引，是"最近点"这条卖点的明显缺口）、**trimmed surface / topology**（CAGD 主场）、**certified / validated numerics**（ISSAC 方向，编辑部里做 interval 的人会看）。
- **反模式警告**：rev19 已有 31 条，再盲目加到 45+ 会被判 citation stuffing。规则：每篇新引用必须在正文有**一句差异化说明**，或在 Table 9 占一行；只出现在 related work 列表里的不加。
- 产物：`docs/experiments/literature_2026/recent_lit_review.md`（必补 ≤12 篇 + 定位句草稿 + DOI 核对 + 检索词），已在进行中。

---

### 第四步：更多开源实现对比 —— 方向对，可选对象比想象少

**现实核验**（重要，避免做无用功）：

| 候选 | 是否有 NURBS/ Bézier 求交 | 可用性 | 结论 |
|---|---|---|---|
| OCCT（IntPatch / GeomAPI_IntSS / IntCS / BRepAlgoAPI_Section） | 有 | 已用（8.0.1，conda pythonocc） | 保留，可再补 IntPatch 原生路径 |
| SINTEF SISL | 曲线求交已用（s1857）；曲面求交例程存在但接口复杂 | 开源、可编译 | **值得补**（step 4 的主力） |
| CGAL | 有 NURBS 曲线/曲面类，但**没有 NURBS SSI**；只有三角网格/AABB | 开源 | 不能做同任务对比，最多做"AABB 预筛"对照，价值低 |
| IRIT（Elber） | 有（业界最强之一） | **非开源**，学术申请制 | 只能 L2 复算，不能运行对比 |
| ACIS / Parasolid | 有 | 商业 | 只能引述文献 |
| FreeCAD / OpenVSP / IfcOpenShell | — | 都调 OCCT | 不能当独立基线 |
| 代数求解器：PHCpack / Bertini / Singular / Macaulay2 | 可解多项式系统（graph-patch SSI 化成 f=0 的多项式系统） | 开源 | **值得补 1 个**：这是 Busé/Mourrain 路线的"可运行代理"，能回答"代数法在同一 fixture 上给几个根" |

**修改建议**：把 step 4 的表述从"更多开源实现"改为"**更多可运行求解范式**"——细分（SISL/OCCT）、代数（PHCpack 或 Bertini 解 f=0）、追踪（OCCT IntPatch）各一个代表 + 矩阵法（Busé–Mourrain）。这既满足 "quantitative"，又避免"凑开源数量"的尴尬。

---

## 4. 路线之外还缺的三件事（决定成败）

### (a) 性能 / 可扩展性叙事

现状：论文最快的数字是**曲线单对** 0.028–1.203 ms，没有任何曲面级/零件级吞吐；真实零件一跑就是 100–550 s/件。
若 step 1 只报"能跑"，审稿人立刻问"为什么这么慢、和 OCCT 差几个数量级"。
**必须主动写**：认证 cover 的成本模型（∝ 交集维数 × (1/τ)^d）、它与 tracing 方法的定位差异（认证 vs 追踪）、以及"工业级需要 BVH + trimmed patch + 局部化，本文给的是可认证内核"。
→ 这其实能把弱点转成 contribution 的边界声明（rev19 已经在结论里写了三条 limitation，这里再加一条"成本"）。

### (b) 证书严格性的一致口径

论文里 "certified" 至少三种强度：
1. **定理级**：§2.4 Krawczyk（`K(X)⊂int(X)`）—— 真证书；
2. **包含级**：§6.3 cover band（严格区间排除，是 no-missing 包装器）—— 真证书，但只对 graph patch；
3. **证据级/启发式**：§6.2 点包含（3×3 采样）、§3.5 κ>10⁴ 分类、§4.1 退化分类（单侧检验）。

CAGD 审稿人一定会抓 "nextafter 膨胀 + Y 怎么保证区间化" 和 "点包含的 3×3 网格不是证明"。
**建议**：新增一小节 + 一张表，明确每个结论的保证等级（Theorem / Enclosure / Heuristic），并统一摘要措辞（"certified root" 只用于 1、2 类）。这一步不增加实验，只增加可信度，性价比最高。

### (c) 定位：5 个算子 vs 1 条主线

rev19 同时讲：曲线求交、最近点、退化恢复、曲面自交、共面重叠、graph-patch SSI。CAGD 的典型审稿意见会是"broad but shallow"。
**建议二选一**：
- **方案 A（推荐，改投 CAGD）**：主线收窄为"**认证 SSI：无遗漏 cover + 自交参数对证书**"，曲线/最近点/退化/重叠作为"同一原语的其他实例"压缩到 1 节 + 1 表。标题相应改（去掉 "Intersection, Closest-Point, Degeneracy and Surface Self-Intersection" 的并列）。
- **方案 B（保留 5 算子，重投 CAD）**：必须补真实数据 + 可运行基线，篇幅会到 30 页，且要正面回应同一主编。

---

## 5. 修订版路线（6 步，可执行，含产物与验收）

| 步 | 做什么 | 产出文件（仓库内） | 验收标准 | 进论文位置 |
|---|---|---|---|---|
| **S0** | 定位与结构决策（本文件 §4c） | `docs/experiments/route_review/route_assessment_2026-11-04.md`（本文件） | 选定方案 A/B | 全文结构、标题、摘要 |
| **S1** | 文献补强（原 step 3） | `docs/experiments/literature_2026/recent_lit_review.md` + `refs_add.bib` | 必补 8–12 篇，每篇一句差异；无 DOI 未核对项 | §1.1、Table 9、参考文献 |
| **S2** | 真实零件（原 step 1，改协议） | `docs/experiments/real_parts/{extract_simjeb.py, run_real_batch.py, results.csv, occt_crosscheck.csv, summary.md, scalability.png}` | ≥12 件有界跑完；OCCT 交叉 ≥600 对；一致率与分歧清单；可扩展性曲线 | 新增 §7.x "Real-part benchmark"，Table 12 |
| **S3** | 可执行基线（原 step 2+4 合并） | `docs/experiments/baselines/{sisl_ssi/, matrix_method/, comparison.md, results.csv}` | 每个范式 ≥1 个可运行基线；同 fixture 结果表；失败模式如实记录 | §7.8 改写，Table 10/11 扩 |
| **S4** | 证书严格性补强（新增） | `docs/experiments/route_review/certificate_levels.md` + 必要的 `Y` 区间化说明/代码 | 保证等级表；摘要措辞统一；Krawczyk 的 Y 处理写清楚 | §2.4 扩、新增 §2.5 或附录 |
| **S5** | rev20 写作（原"更新文章"） | `01_文章正文/P1_CAGD_rev20_*.docx` | 上述实验全部入表；结构按 S0 决策；图表编号连续；页数可控 | 全文 |

**顺序理由**：S1 先做（定义对比对象）→ S2/S3 并行（实验）→ S4 并行（写作性补强）→ S5 汇总。原路线的 2→3 顺序反了：先有"和谁比"，才知道"补谁的文献"。

---

## 6. 需要拍板的三个点

1. **目标期刊**：CAGD（理论向，契合度最高）还是 CAD 重投（范围更合但要正面回应同一主编、需要更强真实数据）？
   建议：**CAGD**，把 CAD 的批评当作"加强实验"的理由写在 cover letter 里。
2. **真实零件的定位**：主证据，还是"可扩展性 + 真实数据存在性"的支撑证据？
   建议：**后者**（主证据仍用解析困难例 + 可运行基线；真实零件证明"不是只会在玩具例子上跑"）。
3. **是否允许收窄主线**（把曲线/最近点从"并列贡献"降为"同一原语的实例"）？
   这会改动标题、摘要、贡献列表 —— 需要作者确认，因为涉及申博材料里的成果叙述。

---

## 7. 风险与备选

| 风险 | 概率 | 备选 |
|---|---|---|
| 真实零件有界跑仍太慢（>10 min/件） | 中 | 只保留 8–10 件小件 + 全部 bbox 不相交对（阴性对照）；把大件写成"scalability limit" |
| SISL 曲面求交 API 不可用 | 中 | 退回 OCCT IntPatch 两条路径 + 矩阵法，SISL 只保留曲线 |
| Busé–Mourrain 矩阵法复现超期 | 低-中 | 只做 curve/surface 最小版本（3 个 fixture），在正文说明"最小复现，非完整实现" |
| 证书严格性补强被审稿人继续攻击 | 中 | 把"certified"限定为"Krawczyk-verified in floating-point interval arithmetic"，并引 Burr/Duff 的同类做法 |
| 时间线 | — | CAGD 无常设 deadline；若为 2027 申博材料，建议 S0–S5 控制在 3 周内，先投再迭代 |

---

## 8. 文献核查后的重大修订：真正风险是"新颖性"，不是"例子数量"（2026-11-04 补充）

在四步之外另做了一次 2020–2026 结构化文献核查（Crossref + OpenAlex，逐条 DOI 可核验），报告见
`docs/experiments/literature_2026/recent_lit_review.md`。结论**改变了对路线优先级的判断**：

### 8.1 三个必须正面处理的先验工作

| 先验工作 | 与本文的关系 | 风险 | 处置 |
|---|---|---|---|
| **Cheng, Zhang, Xiao, Li, TOG 42 (2023)** —— interval 代数拓扑做 rational SSI（rev19 已引 [16]） | 已经用 **interval** 做 SSI，且给拓扑保证 | **高**：若把 "interval + SSI" 当卖点，被判增量不足 | 卖点必须改为"**统一性 + 存在唯一性证书 + 可证无遗漏**"，并明确 Cheng 只给拓扑、不给唯一性证书 |
| **Burr, Hauenstein, Lee, ISSAC '26** —— interval Krawczyk 认证曲面逼近（rev19 已引 [5]） | **同一套 Krawczyk 机器已用于曲面认证** | **最高**：核心机器被先验使用 | 引文里把它写成"曲面逼近/同伦方向"，并强调本文是**几何判定**（求交/自交/重叠）而非逼近 |
| **Yang & Jia, TOG 44 (2025)** —— *Overlap Region Extraction of Two NURBS Surfaces* | **直接覆盖本文"共面重叠"卖点**，且更一般（非仿射） | **高**：贡献 (4) 会被判"已有更好解" | 建议把 coplanar-overlap **从并列贡献降级**为"同一原语的一个受限实例"，引 A3 并说明边界；不要把它写成第 4 个独立贡献 |

另：**Jia, Li, Cheng, CAD 150 (2022)** 是 Thang–Busé–Mourrain 矩阵表示路线的**唯一实质后续**，也是"代数路线"当前真正该被 benchmark 的对象（step 2/4 的基线和引用都要换成它）。

### 8.2 一个可以直接写进 response letter 的核实结论

用 OpenAlex 作者 ID 全量核对（Busé A5088510335 / Mourrain A5003683965 / Kim A5100718911 / Elber A5055421841）：

> **这四位在 2020–2026 没有新的求交论文。** Busé 转向正交投影（SIAM 2020）与双齐次系统（ISSAC 2025）；Mourrain 转向矩/SOS 与认证数值；Kim 组 2022 年后转向 lattice/增材；Elber 近年只有保守核（CAGD 2022）。真正的 state of the art 已转移到 Jia / Cheng / Yang 组与 Fu Xiao-Ming 组。

这条对 cover letter 极有价值：**回应"没有 benchmark Elber/Kim/Busé/Mourrain"时，可以说"我们引用了这四条线的全部相关落点，并进一步与它们路线的最新代表（矩阵法 Jia 2022、自交 Park 2021、重叠 Yang & Jia 2025）做了同输入对比"**——这比硬凑旧文献强得多。

### 8.3 由此产生的路线修订（在 §5 的 S0–S5 之上）

- **S0 增加"新颖性重定位"**：把贡献列从 5 条收成 2+3：
  - 主贡献 A：**认证自交**——用 2×2 方系统 Krawczyk 证明"自交参数对"的存在唯一性（Park 2021 / Li 2025 都只做检测/计算，不做参数对证书）；
  - 主贡献 B：**统一证书机器 + 保证等级表**（Theorem / Enclosure / Heuristic 三级），并明确 no-missing cover 只对 graph patch 成立；
  - 支撑：曲线求交/最近点/退化恢复 = 同一原语的实例（压缩篇幅）；共面重叠 = 受限实例（引 A3，降级）。
- **S3 的基线对象改成**：Jia–Li–Cheng 2022（矩阵法）、Park–Hong–Kim–Elber 2021（自交）、Yang & Jia 2025（重叠）、OCCT IntPatch、SISL。
- **S5 增加"response letter / cover letter"**：单独成文，逐条回应 R1/R2/R3，并附"点名四家 2020–2026 无新求交论文"的核实结论。

### 8.4 需要作者拍板的第 4 个点

**是否拆成两篇**：
- 方案 1（一篇，推荐给 CAGD）：收窄为"认证自交 + 统一证书机器"，共面重叠降级，曲线/最近点压缩。风险最小。
- 方案 2（两篇）：认证自交 → CAGD；统一框架 + 工程实现 → CAD&A / GMP。工作量分散，但每篇都更"纯"。
- 方案 3（保持 5 贡献不减）：必须正面比赢 Yang & Jia 2025 的重叠——以现有仿射实现不可能，不建议。

---

## 附录 A：本次评估用到的证据文件

### 执行回执（本评估之后的实际执行结果）

| 步骤 | 回执文件 | 一句话结论 |
|---|---|---|
| ① 真实零件 | `execution_feedback_step1_realparts.md` + `../real_parts/results/summary.md` | soundness 0 反例（1,404 个 disjoint 对）、cover 覆盖率 100%、但 20 件仅 3 件在 300 s 内跑完 → 只能定位为 local certification kernel；同时发现并修复 harness 的 AABB 负坐标 bug |
| ②④ 可执行基线 | `execution_feedback_step2_4.md` + `../baselines/comparison.md` | SISL 并非全面失效（C 强于 OCCT，A/E/F 相当）；差异在 D 的拓扑碎片化（1e-6 时 1,628 条微曲线）与"无证书"；矩阵法给"快但需重数分析"的 trade-off |
| ③ 文献 | `../literature_2026/recent_lit_review.md` | 点名四家 2020–2026 无新求交论文；真正的 SOTA 在 Jia/Cheng/Yang 组；两个新颖性冲突（Cheng 2023 interval SSI、Burr 2026 Krawczyk 曲面）必须正面处理 |
| 证书审计 | `certificate_levels.md` | 定理级 / 包含级 / 证据级三级；`Y` 的处理可正面回答；4D Krawczyk 只给存在性 |


- 拒稿信：`C:\Users\15632\Doubao\chats\2026-10-13\new-chat\CAD整理\02_审稿与拒稿材料\CAD拒稿通知_20261026.jpg`
- rev19 正文：`01_文章正文\P1_CAD_rev19_清理版.docx`（189 段 / 11 表）；批注版含 85 条批注，其中 12 条是"回应主编"的修改说明
- 既往真实零件工作：`D:\qqdownload\nurbs-certified-framework\docs\experiments\real_parts\`（manifest.csv / selected.txt / data / logs）
- 既往 OCCT 对比：`04_OCCT对比与复现\仓库_occt-comparison\OCCT_COMPARISON.md`（Table 9–11 的实验底稿）
- 期刊规划：`05_文献与规划\刘洋申博成果与期刊规划_2026-10.html`（P1→CAD，备选 CAGD）
- 代码基线：`cpp/include/nurbs/cert_ssi.hpp`（classify_pair / cover / Krawczyk）、`cpp/tests/bench_real_batch.cpp`

## 附录 B：本文件生成方式

- 由 DSH Agent Team 于 2026-11-04 生成；数据来自对上述文件的实际读取与对 WSL 环境的实际运行（g++ 11.4 / OCCT 8.0.1 / SimJEB 382 件）。
- 可复现命令见 `docs/experiments/real_parts/` 与 `docs/experiments/baselines/` 下各 README。
