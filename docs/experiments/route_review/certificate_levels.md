# 证书等级审计（Certificate Levels）— rev19 的 "certified" 口径核对

- 日期：2026-11-04
- 目的：CAGD 审稿人（尤其 interval / symbolic-numeric 方向）会逐条问"这个 certified 到底证了什么"。本文件把 rev19 全文的保证性说法按**可机检强度**分三级，逐条溯源到代码，并给出改稿时的措辞建议。
- 审计对象：`cpp/include/nurbs/cert_ssi.hpp`（2,129 行）+ `cpp/include/nurbs/interval_arith.hpp` + 论文 rev19 §2.3/2.4/3.2/3.4/6.2/6.3。

---

## 0. 结论

- **定理级（Level T）**：曲线求交 2×2 Krawczyk、自交参数对 2×2 Krawczyk —— 这两个是真证书，代码实现与论文陈述一致（含 `Y` 的处理，见 §2.1）。
- **包含级（Level E）**：graph-patch SSI 的 cover band（严格排除 + 资源触顶时保守入 cover）、2D 重合的 cover 交集 —— 是"无遗漏包装器"，但不是存在性证书；且只对 graph patch 成立，rev19 已正确限定。
- **证据级/启发式（Level H）**：共面重叠的点包含（3×3 采样）、κ>10⁴ 相切分类、退化分类单侧检验、`exact_coplanar_overlap` 的阈值判定 —— 这些**不是证明**，rev19 部分已加 caveat，但摘要与贡献列表把它们与 Level T 混在一起说 "certified"。
- **建议**：新增一节（§2.5 "Guarantee levels"）+ 一张表，摘要里 "certified" 只保留给 Level T/E，Level H 一律改 "evidence" / "heuristic" / "detected at the stated tolerance"。这是零实验成本的加分项。

---

## 1. 逐条审计

| # | rev19 的声称 | 代码机制（文件:行） | 实际强度 | 判定 |
|---|---|---|---|---|
| 1 | 曲线求交的每个解由 interval Krawczyk 认证存在且唯一（§2.4/§3.2） | `krawczyk_slice_code` / `probe_box`；`cover` 的 2×2 版本在 `cert_ssi.hpp:681` 附近；`Interval` 每步 nextafter 外扩（`interval_arith.hpp:50+`） | **T**：`K(X)⊂X`（代码 `cert_ssi.hpp:751-762`）+ 收缩判据 | 成立。写作时补一句 `Y` 的处理（§2.1） |
| 2 | 自交参数对 `(s1,s2)` 由 2×2 Krawczyk 证明唯一（§5.2） | 同上，系统 `F(s1,s2)=0` | **T** | 成立；这是全文最硬的一块 |
| 3 | graph-patch SSI 的 cover band 是"可证无遗漏"（§6.3, Proposition 1） | `cover()`：`hull_aabb` 不相交即排除（`cert_ssi.hpp:352`）；资源触顶时**保守地把当前盒并入 cover**（`:342-348`） | **E**：排除可靠 ⇒ 无漏；但 cover 不是存在性证明 | 成立，但措辞要改："no-missing **enclosure**"，不是 "certified intersection" |
| 4 | 2D 重合（type=2）由点包含 + 法向判定（§6.1） | `coincidence_ratio`：在 [0.15,0.85]² 采 2×2 点比对残差 | **H**：有限采样，非证明 | 必须降级 |
| 5 | 共面重叠由 3×3 点包含证据支持（§6.2） | 3×3 内缩网格 + 点包含谓词 | **H**（rev19 已限定 affine） | 降级；且注意 Yang & Jia TOG 2025 已做一般 NURBS 重叠域提取（见 `literature_2026/recent_lit_review.md` A3） |
| 6 | 退化分类 sound（§4.1） | `classify_degeneracy`：控制点/权重阈值单侧判定 | **H**：触发即可靠（恰好重合时），不触发不能证非退化 | rev19 已写 caveat，保留 |
| 7 | κ>10⁴ 判相切（§3.5） | 根处 Jacobian 条件数 | **H**：标签而已 | rev19 已写 caveat，保留 |
| 8 | "zero false positives / zero false certificates" 全regression成立 | 各 regression 脚本 | **E**（在测试集上） | 保留，但注明是"在报告的全部 regression 上"，不是"对所有输入" |

---

## 2. 两个最容易被审稿人抓的技术点（提前处理）

### 2.1 `Y` 是怎么处理的？—— 这是可以正面回答的

代码里 `Y` 由**中点 Jacobian 的普通 double 求逆**得到，再以 `Interval::point(Y[a][c])` 参与运算（`cert_ssi.hpp:717-720`、`:822-834`；`krawczyk_slice_code` 的同名逻辑在 `:717-720`）。

**这并不破坏证书**：经典 Krawczyk 定理对**任意实矩阵** `Y` 成立，只要求它是实矩阵、可计算。`Y` 取"近似的数值逆"正是该定理的设计意图（不需要区间求逆）。因此：
- 证书成立性：不受影响；
- 需要说明的是：`J3inv` 的浮点误差只影响 `Y` 这个实数的取值，不影响 `K(X)` 是 `g(x)=x-YF(x)` 的区间像这一事实。
- **建议**：正文明确写 "Y is an arbitrary real matrix computed from the midpoint Jacobian in double precision; the Krawczyk theorem does not require Y to be an interval enclosure of the inverse." 一句话就能堵住最常见质疑。

参考：Neumaier (1990) §5；Moore–Kearfott–Cloud (2009) Ch. 8；Breiding–Rose–Timme (TOMS 2023) 对"认证"术语的用法。

### 2.2 4D 欠定 Krawczyk 只给存在性

`krawczyk_4d_code`（`cert_ssi.hpp:798+`）对 3×4 满秩系统用 Brouwer 不动点给"存在一个零点"，**没有唯一性**（沿 1D 解流形本来就不唯一）。rev19 已声明"只用 2×2 方系统"，与代码一致。若在附录/补充材料里用到 4D 版本，必须写明"existence only"。

---

## 3. 建议加入论文的表格（可直接用）

**Table X. Guarantee levels of the reported statements.**

| Statement | Mechanism | Level | Scope |
|---|---|---|---|
| Curve intersection root | 2×2 interval Krawczyk, `K(X)⊂int X` | Theorem (existence + uniqueness, floating-point interval arithmetic) | NURBS/Bézier curve pairs |
| Shortest-distance candidate | stationary-point certification + min over certified candidates | Theorem (each stationary point certified); global min by enumeration | curve pairs |
| Self-intersection parameter pair `(s1,s2)` | 2×2 interval Krawczyk | Theorem (existence + uniqueness) | lemniscate-cylinder family (random reparameterization + affine transforms) |
| SSI cover band (graph patches) | strict interval/hull exclusion, no box discarded without a proof of absence | Enclosure (no-missing superset of the intersection set) | graph patches only |
| Regular patch is fold-free | all boxes excluded | Enclosure (proof of absence of a fold) | lemniscate family |
| Coplanar overlap region | convex-hull AABB + 3×3 point-inclusion evidence | Evidence (affine/coplanar only) | affine planar patches |
| Tangency vs. transversal label | condition number threshold 10⁴ | Heuristic label (does not affect certificates) | all certified roots |
| Degenerate entity class | one-sided control-point/weight test | Sound when triggered; not a non-degeneracy proof | — |

---

## 4. 改稿措辞清单（逐句替换建议）

| 位置 | rev19 现文 | 建议改为 |
|---|---|---|
| Abstract | "every solution is certified by the classical interval Krawczyk test" | 保留，但后接 "…; the exclusion-based SSI cover and the coplanar-overlap module provide enclosure-level and evidence-level statements respectively (Section 2.5)." |
| Contributions (3) | "Certified detection of genuine surface self-intersection …" | "**Theorem-level** certification of a genuine self-intersection parameter pair" |
| Contributions (4) | "Certified detection of coplanar patch overlap …" | 降级并入贡献 (2) 的实例；或改为 "Evidence-level detection … restricted to affine patches" |
| §6.2 | "certifies a common area" | "supports a common area"（rev19 已有 caveat，统一到全文） |
| §6.3 | "no-missing-point certificate" | "no-missing **enclosure** certificate"（避免被读成存在性） |

---

## 5. 若要进一步加固（可选，2–3 天）

1. 把 `hull_aabb` 排除换成**区间残差排除**（`f` 的区间不含 0 才排除）——更紧，且和 Proposition 1 的文字完全对齐。当前实现是控制点 AABB，不是残差区间，文字上容易被抓。
2. 给 `cover` 加一个"排除盒数 vs 覆盖盒数 vs 触顶盒数"的输出（现在 `cap_hit` 有标志但没在论文里报），并把"触顶盒一律入 cover"写进 Proposition 1 的假设（保证 no-missing 仍然成立）。
3. 对 §6.2 的点包含，换成一个真正的包含证书（如区间 Newton 证明参数存在），而不是 3×3 采样。

---

## 附录：审计所读代码位置

- `cpp/include/nurbs/interval_arith.hpp:1-60`（nextafter 外扩策略）
- `cpp/include/nurbs/cert_ssi.hpp:331-398`（cover / CoverResult / node_cap 保守处理）
- `cpp/include/nurbs/cert_ssi.hpp:681-777`（切片 2×2 Krawczyk）
- `cpp/include/nurbs/cert_ssi.hpp:798-900+`（4D 欠定 Krawczyk，existence only）
- `cpp/include/nurbs/cert_ssi.hpp:1928-2074`（exact_coplanar_overlap / coincidence_ratio / coplanar_candidate）
- `cpp/include/nurbs/cert_ssi.hpp:2075-2129`（classify_pair：类型分派与预算）
