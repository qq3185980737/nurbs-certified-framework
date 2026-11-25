# 矩阵表示法基线（minimal reproduction）

**定位声明（必须先读）**

本目录的矩阵方法实现是 **minimal reproduction（最小复现）**，**不是**
Thang–Busé–Mourrain 2009 或 Jia–Li–Cheng, *Computer-Aided Design* 150:103303
(2022) 的完整实现。完整构造需要从两张有理 Bézier 曲面的 moving surfaces /
Dixon resultant 构造多项式矩阵表示（bidegree (2,2) 的隐式化与消元机器），
本会话内未实现。这里复现的是矩阵表示法的**最后一步——"矩阵 → 广义特征值"**
在 graph-patch 特化下的可运行版本；所有数字只能标注为
"**minimal reproduction**"，不能当作原论文方法的性能数字引用。

- 脚本：`scripts/matrix_baseline.py`（numpy 2.4.6，无 scipy 依赖）
- 结果：`results/matrix_baseline_raw.txt`、`results/matrix_curve_surface.csv`、
  `results/matrix_certified_1d.csv`、`results/matrix_surface_surface.csv`

## 1. 两个子实验

### (1) 3D Bézier 曲线 × 双三次 graph 曲面（TB–M 2009 的 curve/surface 场景，特化）

曲面 `S(u,v) = (u, v, f(u,v))`，`f(u,v) = 0.5u³ − 0.4v³ + 0.15uv + 0.05`（双三次）。
曲线 `C(t)` 写成齐次形式 `(X(t),Y(t),Z(t),W(t))`。交点条件：

```
Z(t)·W(t)^{p+q−1} − Σ a_ij·X(t)^i·Y(t)^j·W(t)^{p+q−i−j} = 0
```

得到一个单变量多项式 `g(t)`——**这就是该特化下的"matrix representation"的系数**；
其根用**伴随矩阵的特征值**（`numpy.linalg.eig`，脚本里显式构造 companion matrix）
一次性求出，再对落在 [0,1] 的实特征值做 Newton 精化并计算
`|C(t) − S(x(t),y(t))|`。

4 个 fixture（`results/matrix_curve_surface.csv`）：

| case | 真值根数 | 实根数 | 去重后 | 全部特征值数 | 最大残差 | 耗时 ms | 根 |
|---|---|---|---|---|---|---|---|
| CS1_transverse | 1 | 1 | 1 | 3 | 0.0 | 0.23 | 0.005983 |
| CS3_three_roots | 3 | 3 | 3 | 3 | 4.2e-17 | 0.05 | 0.2 / 0.5 / 0.8 |
| CS1_tangential_double | 1（二重根） | **2** | 1 | 2 | 2.8e-17 | 0.03 | 0.5 / 0.5 |
| CS1_rational_arc（有理） | 1 | 1 | 1 | **12** | 2.2e-16 | 0.09 | 0.229367 |

观察：
- 横截/多根情形：特征值法一次给出全部根，残差 ~1e-17，耗时 <0.25 ms。
- **重根（相切）情形：二重根在数值上变成两个几乎相同的特征值（0.5, 0.5），
  必须做 cluster 去重才知道"1 个根"**——矩阵法本身不直接给重数语义。
- **有理曲线：清分母（乘 `W^{p+q}`）引入了多余特征值**（12 个里只有 1 个是真根，
  其余是 `W(t)=1+t²` 的复根/无穷远根），需要按虚部与残差过滤。原论文用更精细的
  matrix representation 避免这类 extraneous factors；这是本最小复现的已知缺陷。

### (1b) 同一多项式上的 certified 区间计数（对照）

用与 `cpp/tests/bench_ssi.cpp` 相同的严格区间算术（`nextafter` 向外）+ 二分排除，
对同一个 `g(t)` 做**certified 计数**（覆盖盒的连通分量数 = 根簇数）：

| case | 覆盖盒 | 排除盒 | certified 根簇 | 耗时 ms |
|---|---|---|---|---|
| CS1_transverse | 1 | 34 | 1 | 0.7 |
| CS3_three_roots | 21 | 622 | 3 | 10.6 |
| CS1_tangential_double | **185364** | 262142 | 1 | **7014.6** |
| CS1_rational_arc | 1 | 34 | 1 | 2.5 |

对照结论（可直接写进论文）：**矩阵特征值法在相切/重根上极快（0.03 ms）但需要
重数分析；certified 区间法给出认证根簇数，但在二重根上代价爆炸（7 s、18 万覆盖盒）
——这正是"认证 vs 一次求全部根"的 trade-off，二者互补。**

### (2) 2×2 / 3×3 graph-patch 曲面×曲面（Jia 2022 的"2×2"最小版本）

两张 graph patch 的交线是 `{g(s,t)=0}`，`g = f1−f2`。对每个固定的 `t`，把
`g(·,t)` 看成 `s` 的多项式（deg = 2 对应 2×2 伴随矩阵，deg = 3 对应 3×3），
其特征值就是 `s` 的根。201 个 `t` 采样（`results/matrix_surface_surface.csv`）：

| fixture | 真值 | 实根总数 | 曲线点 | 复数特征值 | 最大残差 | 耗时 ms | 每 t 根数 |
|---|---|---|---|---|---|---|---|
| A | 1 条对角线 | 202 | 202 | 0 | 0.0 | 18.3 | 1–2 |
| B | 0 | 0 | 0 | 0 | — | 0.7 | 0 |
| C | 1 个角点 | 1 | 1 | 0 | 0.0 | 5.1 | 0–1 |
| D | 1 条切触对角线 | **402** | 402 | 0 | 1.1e-16 | 25.6 | **2** |
| E | 1 条内部曲线 | 93 | 93 | 402 | 2.8e-17 | 15.5 | 0–1 |

观察：
- A/B/C/E 的零点集合被正确找出（A 在 t=0 处 `s²=0` 是二重根，多算 1 个点）。
- **D（退化切触）每 t 都返回 2 个相同特征值**（`(s−t)²` 的二重根），
  点数直接翻倍（402 vs 201 个参数位置）；**必须做重数/聚类分析才能得到
  "一条曲线"**。矩阵法对重根的处理是"特征值重数"，与几何上的切触语义不一致。
- E 中约 2/3 的 `t` 只有复特征值（该 t 处无实交点），被正确过滤。

## 2. 为什么没有做完整的 Jia 2022 / TB–M 2009 实现（如实记录）

1. **需要隐式化/消元机器**：Jia 2022 的核心是为一对有理曲面构造 matrix
   representation（moving surfaces / Dixon resultant），其大小与 bidegree 有关，
   需要多元结式与符号消元；用 numpy 从零实现的工作量远超本任务预算。
2. **graph-patch 特化已经覆盖了论文的判断需求**：论文 fixture A–E 的求交条件就是
   `f1−f2=0`，其"矩阵表示"在特化下退化为单变量/单参数多项式，正好用来验证
   "特征值一次求全部根"这一机制及其**重数缺陷**。
3. 因此本基线只用于**机制对照与定性讨论**，不进入论文的运行时对比表；
   若审稿要求，可在后续工作中补完整 Jia 2022 实现。

## 3. 可比性边界

- 本实现 **exploits graph-patch 结构**（直接得到高度场 `f`），而 SISL/OCCT 是
  通用曲面求交器（不知道 graph 结构）；因此**耗时不可与 SISL/OCCT 直接比较**，
  只能比较"根数/重数/认证能力"。
- 无任何存在性/唯一性/完备性证明；有理清分母的多余根需要经验过滤。
- 重根（切触）情形需要额外的聚类与重数分析，数值上病态（二重根的
  特征值分裂 ~√ε）。
