# SISL 4.7 曲面求交基线（可运行，实测）

本文件记录 **SISL 4.7** 在论文共享 fixture 上的曲面-曲面 / 曲面-曲线求交基线。
所有数字均来自本目录下的程序实际运行，原始输出见 `results/`。

## 1. 环境与来源（可复现）

| 项 | 值 |
|---|---|
| OS | WSL2 Ubuntu 22.04（8 核，15 GB） |
| 编译器 | g++ (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0 |
| CMake | 3.22.1 |
| SISL 版本 | **4.7**（`ChangeLog` 首条：2021-03-16 "Released version 4.7 of SISL"） |
| 源码位置 | `/mnt/d/sisl_dl/src/SISL-master`（SINTEF/SISL GitHub master 的 zip 解压） |
| 构建产物 | `sisl_build/libsisl.a`（2.75 MB，`-O2 -fPIC`） |

说明：本机 WSL **未安装 git**（`git: command not found`），因此没有执行
`git clone`，而是复用已存在的 SISL 4.7 源码（与 GitHub master 一致，版本由
`ChangeLog` 证实）。网络本身可用（`curl https://github.com` → HTTP 200）。

构建（见 `scripts/build_sisl.sh`）：

```bash
cmake -S /mnt/d/sisl_dl/src/SISL-master -B docs/experiments/baselines/sisl_build \
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="-O2 -fPIC" \
      -Dsisl_COMPILE_APPS=OFF -Dsisl_COMPILE_EXAMPLES=OFF -Dsisl_COMPILE_VIEWER=OFF
cmake --build docs/experiments/baselines/sisl_build -j 8
```

## 2. SISL 的求交 API（读 `include/sisl.h`、`doc/manual/func/*.tex`、`app/surf_surf_int.C` 实证）

| 函数 | 作用 | 关键输出 |
|---|---|---|
| `s1859(surf1, surf2, epsco, epsge, &npt, &par1, &par2, &ncr, &intcurves, &stat)` | 曲面×曲面交集**拓扑** | 孤立点数 `npt`、交线条数 `ncr`、每条交线的参数域引导点（`epar1/epar2`） |
| `s1310(surf1, surf2, intcurve, epsge, maxstep, makecurv=2, graphic, &stat)` | 把交线**追踪**成 3D B-spline + 两个参数域曲线 | `pgeom/ppar1/ppar2`；`stat=3` 表示遇到奇异点/退化曲面 |
| `s1858(surf, curve, epsco, epsge, &npt, &par1, &par2, &ncr, &intcurves, &stat)` | 曲面×曲线 | 交点参数（曲面 2 个、曲线 1 个） |
| `s1857(curve1, curve2, epsco, epsge, ...)` | 曲线×曲线 | 交点参数 |
| `s1850 / s1851 / s1852 / s1853 / s1854` | 曲线×平面/直线、曲面×平面/球/柱/锥 | 同类接口 |
| `newSurf(n1,n2,k1,k2,knot1,knot2,coef,kind,dim,copy)` | 构造 Bézier/B-spline 曲面（`kind=3` 多项式 Bézier，`kind=4` 有理 Bézier） | — |

**结论：SISL 4.7 有完整的曲面-曲面求交高层例程（`s1859` + `s1310`），可以
直接作为基线，无需退化为"不可行"。** 仓库自带 `app/surf_surf_int.C` 是官方调用示例。

Graph patch（论文 fixture）用 degree 2×2 Bézier 曲面精确表示：
`S(u,v) = (u, v, f(u,v))`，控制网 x 坐标取 `i/(n-1)`、y 取 `j/(n-1)`（Bézier 基下
`Σ (i/n) B_i(u) = u` 精确成立），z 为 `f` 的 Bernstein 系数。

## 3. 主表：fixture A–F（同输入、同容差）

程序：`src/sisl_ssi_harness.cpp`；构建/运行：`scripts/build_run_sisl.sh`；
结果：`results/sisl_ssi.csv`（摘要）、`results/sisl_ssi_detail.txt`（逐条引导点/追踪曲线端点）。

| fixture | 真值 | epsge | s1859 状态 | 交点数 | 交线数 | 引导点 | s1310 错误 | 残差 | 耗时 ms |
|---|---|---|---|---|---|---|---|---|---|
| A `f=s²−t²` | 1 条对角线 s=t | 1e-6 | 0 | 0 | 1 | 3 | 0 | 1.1e-16 | 0.77 |
| B `f=0.3` | 0 交 | 1e-6 | 0 | 0 | 0 | 0 | 0 | — | 0.002 |
| C `f=s+0.5t` | 1 个角点 (0,0) 切触 | 1e-6 | 0 | **1** | 0 | 0 | 0 | 0 | 0.16 |
| D `f=(s−t)²` | 1 条退化切触对角线 | 1e-6 | 0 | **346** | **1628** | 3256 | **1536** | 9.5e-7 | **2923.6** |
| E 双三次 `u³+v³=0.1` | 1 条内部交线 | 1e-6 | 0 | 0 | 1 | 2 | 0 | 4.1e-7 | 0.76 |
| F 有理四分之一圆柱 × 水平直线 | 1 交点 (0.866025,0.5,0.5) | 1e-6 | 0 | **1** | 0 | 0 | 0 | 1.1e-16 | 0.06 |

残差定义：
- 单点：`|S1(p1) − S2(p2)|`（在 s1859 返回的参数处求值）；
- 交线：对 s1310 追踪出的 3D 曲线密集采样，直接在 `(x,y)` 处求值
  `max(|P−S1(x,y)|, |P−S2(x,y)|)`（graph patch 的 x,y 即两张曲面的参数）。

关键观察：

1. **横截构型（A、B、E、F）SISL 全部正确且极快**（0.06–0.8 ms，残差 ≤ 4e-7）。
2. **C（角点切触）SISL 正确返回 1 个孤立交点**——这是它比 OCCT 的
   `GeomAPI_IntSS` 强的地方（OCCT 该 API 的公开结果只有曲线）。
3. **D（退化切触，`f=(s−t)²`）SISL 不是"返回 0"，而是把一条解析上光滑的
   对角线炸成 1628 条微曲线 + 346 个孤立点，其中 1536 条追踪失败
   （`s1310` 返回错误），耗时 2.9 s。** 在 epsge=1e-4 时为 70 条曲线 + 54 点、
   116 ms。几何上每个输出点都落在对角线上（`max|u−v| ≈ 2e-7`、残差 ~1e-14），
   但**拓扑完全不可用，且没有任何"这条对角线是一个整体、存在且唯一"的证书**。
4. epsge 从 1e-6 放宽到 1e-4，D 的碎片数从 1628 降到 70——**输出拓扑随容差
   剧烈变化**，说明它不是一个稳定的几何算子。

`results/sisl_ssi.csv` 同时给出 epsge=1e-4 的整行数据，可对照。

## 4. 失效构型证据（论文中"无证书/不完整"的直接证据）

程序：`src/sisl_failure_probe.cpp`；运行：`scripts/run_sisl_failure.sh`；
结果：`results/sisl_failure_probe.txt`。

### 4.1 多根：5 次 Bézier 曲线 × 直线（`s1857`）

`C(t) = (t, 30·∏_{k=1..5}(t−k/6))` 与 x 轴在 (0,1) 内有 5 个横截交点。

```
stat=0  num_points=5  num_curves=0   (ground truth: 5 simple roots)
root[0] t=0.500000000  C(t).y=-6.6e-17
...
root[4] t=0.833333333  C(t).y=-2.8e-17
max |C(t).y| = 6.592e-17
```

SISL 找齐 5 个根、残差 ~1e-17；但**它只返回几何交点，不返回根的重数/存在性，
也不给根数上界或完备性证明**。论文 §7.7(i) 的对照可写"根数一致，但 SISL 无证书"。

### 4.2 曲面自交

SISL **没有曲面自交算子**。对 `s1859(S,S)`（同一张曲面）探测：`stat=0，
npt=0, ncr=4`——重合代码路径返回参数方形的 4 条边界当作"交线"，与自交无关。
→ 论文可写"SISL 不提供曲面自交（self-intersection）算子，本文的
self-intersection 定位是能力差异而非速度差异"。

### 4.3 二维重合（coplanar overlap）

程序：`src/sisl_overlap_probe.cpp`；结果：`results/sisl_overlap_probe.txt`。

| 构型 | SISL `s1859` 输出 | 解读 |
|---|---|---|
| G 两张完全相同 patch（全重合） | `stat=0, ncr=4`，四条参数域边界 | 把"2D 重合区域"的**边界**当成交线返回，**不报告"交集是 2 维区域"** |
| G2 两张相同的 z=0 平面 | 同上（4 条边界） | 同上 |
| H 共面平板部分重叠（x∈[0,1] 与 [0.5,1.5]） | `ncr=4`，恰好是重叠区域的边界 | 同上；用户无法从输出契约区分"两条交线"与"二维重叠" |
| H2 共面平板仅共边（x=1 接触） | `ncr=1`，x=1 的一条线段 | 1D 交集，正确 |
| J 近共面 `z=0` vs `z=1e-7u`（低于 epsge） | `ncr=4` | 容差内视为重合，退化为边界环 |

**这是 SISL 基线最重要的可比性边界：SISL 没有"二维重合"这一输出类别，
只能退回边界环（且 `stat=0`，不报警）。** 论文的 coplanar-overlap 分类结果
（仿射片成立）是 SISL 输出契约里根本不存在的信息。

### 4.4 未做的构型（如实记录）

- **小环/近切触（CAGD 2026 Yang & Jia 那类小环构型）**：本会话未构造。原因：
  可靠的 3D 小环 fixture 需要额外设计（两曲面在很小参数区域内二次相切形成闭合环），
  且 SISL 在切触/退化构型上的行为已被 D（§4 主表）与重合探测（§4.3）覆盖；
  时间优先级低于 OCCT 重跑与矩阵基线。若审稿需要，可作为后续补充。

## 5. 可比性边界（写进论文的表述）

1. SISL 是**几何求解器**，返回近似曲线/点，**不对根的存在性、唯一性、完备性作证明**；
   容差 `epsge` 只影响细分/追踪，不构成证书。
2. 输出拓扑是**不稳定的**：退化构型 D 的"交线数/点数"随 `epsge` 从 70/54 变到
   1628/346，因此"交线数"不能直接当作方法的正确性指标。
3. SISL 不对 2D 重合给出分类（见 §4.3），也不提供曲面自交算子（§4.2）。
4. 在简单横截构型上 SISL 非常快（亚毫秒级），**本文不应声称速度优势**；
   本文的差异是**证书与无遗漏**（cover 盒带 ⊇ 解集）。
