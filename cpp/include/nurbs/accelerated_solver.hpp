// ============================================================================
// accelerated_solver.hpp — BVH 加速求解（空间索引减少 LM 初值）
// ----------------------------------------------------------------------------
// 独立于 unifiedSolve 的加速路径：
//   对每条曲线建立"参数分段 BVH"（叶节点 AABB = 该参数段控制点凸包盒），
//   对大量曲线对做求交时：每个曲线索引只需构建一次，随后每次查询仅做
//   两树遍历 + 对相交叶对跑少量 LM 初值（替代全局 6×6 网格 + 候选区）。
//   可复用索引使建树开销被摊薄，从而显著提升批量求解吞吐。
// 语义与 unifiedSolve 一致（相交/相切/分离），解集合在容差内一致；
// 初值不同导致收敛轨迹不同，不保证逐位一致。
// 原 unifiedSolve 保持不动，以维持与 JS 的逐位等价（见 tests/test_equivalence）。
// ============================================================================
#ifndef NURBS_ACCELERATED_SOLVER_HPP
#define NURBS_ACCELERATED_SOLVER_HPP

#include "unified_solver.hpp"
#include "interval_arith.hpp"

namespace nurbs {

// 曲线空间索引（BVH）。构建后可在任意多次求交查询中复用。
class BVHIndex {
public:
    BVHIndex() = default;
    void build(const NURBSCurve& curve, int maxDepth = 7);
    // 根节点 AABB（整条曲线的空间包围盒；构建后有效），供 sweep-and-prune 等场景使用
    const AABB& rootBox() const { return nodes[0].box; }
    // 导出叶节点 14-DOP 与参数区间（供 GPU / 外部并行遍历使用）。
    void exportLeaves(std::vector<double>& dopOut, std::vector<std::pair<double, double>>& params) const {
        dopOut.clear(); params.clear();
        for (const auto& nd : nodes) if (nd.leaf) {
            for (int d = 0; d < 14; ++d) dopOut.push_back(nd.dop[d]);
            params.emplace_back(nd.uStart, nd.uEnd);
        }
    }

private:
    friend std::vector<SolveResult> solveWithIndex(
        const BVHIndex&, const NURBSCurve&, const BVHIndex&, const NURBSCurve&,
        const SolverOptions&, int);
    friend BVHIndex buildIndex(const NURBSCurve&, int);
    friend void collectLeafPairs(const BVHIndex&, int, const BVHIndex&, int,
                                 std::vector<std::pair<int, int>>&);
    friend std::pair<size_t, size_t> adaptiveInitCounts(
        const BVHIndex&, const NURBSCurve&, const BVHIndex&, const NURBSCurve&);
    friend size_t dopLeafPairCount(const BVHIndex&, const NURBSCurve&, const BVHIndex&, const NURBSCurve&);
    friend size_t aabbLeafPairCount(const BVHIndex&, const NURBSCurve&, const BVHIndex&, const NURBSCurve&);
    int buildRec(const NURBSCurve& curve, double a, double b, int depth, int maxDepth);
    // box: AABB；dop: 14-DOP（7 个方向 min/max，比 AABB 更紧，用于相交精剪）
    struct Node { double uStart = 0, uEnd = 1; AABB box{}; double dop[14]; bool leaf = false; int child0 = -1, child1 = -1; };
    std::vector<Node> nodes;
};

// 统计 14-DOP 精剪后的相交叶对数（对照纯 AABB 叶对 aabbLeafPairCount，量化包围体收紧）。
size_t dopLeafPairCount(const BVHIndex& ia, const NURBSCurve& c1, const BVHIndex& ib, const NURBSCurve& c2);
// 纯 AABB 相交叶对数（AABB 基线）。
size_t aabbLeafPairCount(const BVHIndex& ia, const NURBSCurve& c1, const BVHIndex& ib, const NURBSCurve& c2);

// 构建索引（等价于 BVHIndex::build 后返回）
BVHIndex buildIndex(const NURBSCurve& curve, int maxDepth = 7);

// 用预建索引求解两条曲线的全部解（曲线与索引需对应）。
// numThreads：0=自动（硬件并发），1=强制串行，>1=指定线程数。
std::vector<SolveResult> solveWithIndex(const BVHIndex& ia, const NURBSCurve& c1,
                                        const BVHIndex& ib, const NURBSCurve& c2,
                                        const SolverOptions& options = {},
                                        int numThreads = 0);

// 便捷入口：每次调用内部建索引再求解（适合单次/偶发查询）
std::vector<SolveResult> unifiedSolveAccelerated(const NURBSCurve& c1, const NURBSCurve& c2,
                                                 const SolverOptions& options = {});

// 自适应初值策略统计（供论文/基准量化初值削减）：返回 {相交叶对数, 自适应初值总数}。
// 对照"固定 9 初值"时为 9×叶对数；二者差值即自适应削减的初值（若不丢解则为净收益）。
std::pair<size_t, size_t> adaptiveInitCounts(const BVHIndex& ia, const NURBSCurve& c1,
                                             const BVHIndex& ib, const NURBSCurve& c2);

// ---- 批量多曲线求交（任务2：BVH 森林 + sweep-and-prune 候选对筛选 + 并行） ----
// 单次求交结果：曲线对 (i,j) + 该对全部解
struct BatchIntersection { int i, j; std::vector<SolveResult> r; };

// 用根 AABB 做 sweep-and-prune：按 minX 排序、滑动窗口筛出 AABB 相交的候选曲线对，
// 避免全 N² 两两求交。boxes[i] 与 curves[i] 一一对应。
std::vector<std::pair<int, int>> sweepAndPrunePairs(const std::vector<AABB>& boxes);

// 集合内所有曲线两两求交（BVH 森林复用 + 候选对级别并行）。
// numThreads：0=自动，1=串行，>1=指定。返回每条相交曲线对的最佳解。
std::vector<BatchIntersection> batchSolveAll(const std::vector<NURBSCurve>& curves,
                                             const SolverOptions& options = {},
                                             int numThreads = 0);

// ---- 保证性求交（区间/Krawczyk 解存在性证明）----
// 对收敛解 (u,v) 的邻域盒 [ub]×[vb]，用区间算术扩展 + Krawczyk 算子证明该盒内
// 存在唯一解（严格舍入区间，NURBS 凸包区间求值）。返回 true 表示已证明。
// 仅对横截（良态）解适用；相切/退化（Jacobian 奇异）返回 false（正确行为）。
bool proveUniqueInBox(const NURBSCurve& c1, const NURBSCurve& c2,
                      const Interval& ub, const Interval& vb);

// GS-Newton 区间证明（Hansen-Sengupta Gauss-Seidel 变体）。把线性化系统 A(x−m)=−Y f(m)
// 逐变量用区间消元更新（比 Krawczyk 一步算子更紧），允许在更大的盒上证明存在唯一，
// 或同等覆盖率下更少区间运算。返回：1=证明存在唯一；0=未定（继续细分）；-1=证明无解（排除）。
// 成功时 ub/vb 被内缩到证明盒。
int proveUniqueGS(const NURBSCurve& c1, const NURBSCurve& c2, Interval& ub, Interval& vb);

// 相切交点精化（对标《计算机辅助设计与图形学学报》2025 有理插值处理相切的思想）：
// 对相切/近奇异交点（Jacobian 奇异，LM 线性收敛、Krawczyk 不适用），用阻尼
// Gauss-Newton 最小化 ||c1(u)−c2(v)||²，将切点定位精化到机器精度。
// 返回精化后的残差距离（越小越精）。
double refineTangentPoint(const NURBSCurve& c1, const NURBSCurve& c2,
                          double& u, double& v, int iters = 80);

// ---- B/C: 解分类 + 条件数/误差界 ----
enum class IntersectionKind { TRANSVERSE, TANGENT, OVERLAP, UNKNOWN };
struct ClassifiedSolution {
    double u, v;                 // 交点参数
    double distance;             // 几何距离（backward error ||f||）
    double det;                  // 交点处 Jacobian 行列式
    double cond;                 // 条件数估计 κ(J)（F 范数 / |det|）
    IntersectionKind kind;       // 交点类型
};
// 对求解结果做后置分类：|det(J)| ≥ ε → 横截；|det(J)| < ε → 相切（近奇异）。
// 同时给每个交点 Jacobian 条件数与 backward error，量化结果可信度。
std::vector<ClassifiedSolution> classifySolutions(const NURBSCurve& c1, const NURBSCurve& c2,
                                                  const std::vector<SolveResult>& sol);
// 检测两曲线是否存在重合段（退化重叠）。基于控制点一致性 + 公共参数区采样距离。
bool detectOverlap(const NURBSCurve& c1, const NURBSCurve& c2, double tol = 1e-6);

// ---- ④ 全局保证性求交：区间细分全局枚举 ----
// 曲线在参数区间 [u] 上的区间扩展（两套，各司其职）：
//   curveInterval      : 一致收敛（中心+Lipschitz，δ→0 收敛到点值）→ 叶盒 Krawczyk 证明用。
//   curveIntervalHull  : 凸包严格区间（控制点 min/max，严格包含真值，不低估）→ 全局剪枝/排除用，
//                        保证"未被排除的盒覆盖所有真实交点"（无遗漏）。
struct CurveIntervalEval { Interval fx, fy, dfx, dfy; };
CurveIntervalEval curveInterval(const NURBSCurve& c, const Interval& u);
CurveIntervalEval curveIntervalHull(const NURBSCurve& c, const Interval& u);
// 严格一致收敛区间（全局排除用）：f(m) ± L·δ，L 用导数控制点凸包范数上界构造
// （对有理 NURBS 含分母界），保证 c(u) ⊂ c(m) ± L·δ 严格成立（不漏解），
// 又比凸包区间紧（随 δ→0 收敛）→ 排除测试快且严格。
CurveIntervalEval curveIntervalStrict(const NURBSCurve& c, const Interval& u);

// 全局枚举解：叶盒解（含证明盒、分类）
struct GlobalIntersection {
    double u, v;          // 证明盒中心（交点近似）
    Interval ub, vb;      // 证明盒（参数域）
    bool proven;          // 是否已被 Krawczyk 证明存在唯一
    IntersectionKind kind;// 分类（横截已证明 / 疑似相切退化）
};
struct GlobalEnumerationResult {
    std::vector<GlobalIntersection> solutions;
    int provenCount = 0;       // 已被 Krawczyk 证明存在唯一的横截解数
    int tangentCount = 0;      // 疑似相切/退化（Krawczyk 失败且近奇异）解数
    int evaluatedBoxes = 0;    // 细分评估的盒总数
    bool complete = false;     // 全参数域已覆盖
};
// 从全参数域根盒出发，区间细分全局枚举两 NURBS 曲线交点：
//   凸包严格区间做排除测试剪枝（不漏解），叶盒用 Krawczyk 证明存在唯一（无漏 + 唯一），
//   相切/退化（证明失败）按条件数分类保留。结果决定性与启发式不同，可直接作为真值上界。
GlobalEnumerationResult enumerateIntersections(const NURBSCurve& c1, const NURBSCurve& c2,
                                               double leafTol = 1e-4);

// ---- ⑤ 退化/奇异曲线处理（结合 Monterde NURBS 奇异点结构定理、CASC2023 尖点分析）----
// 退化检测：控制点全重合→POINT（退化点）；控制点共线→COLLINEAR（退化到线段）；
//           权重 wmin 过小→RATIONAL（有理投影退化，Lipschitz wmin² 分母界爆炸）。
// 这些退化情形使常规 Krawczyk（Jacobian 奇异）失效或求解器漏解/噪声解，
// 需专门处理（退化点→点包含证明；共线→降维；RATIONAL→Box 收缩定位）。
enum class CurveDegeneracy { NONE, POINT, COLLINEAR, RATIONAL };
CurveDegeneracy detectDegeneracy(const NURBSCurve& c, Point2* outPoint = nullptr, double* outWmin = nullptr);
// 退化点 P × 曲线：点包含证明（广义 Moore-Penrose Krawczyk，超定 R→R²，A=G⁻¹Jᵀ）。
// 证明存在唯一 v∈vb 使 c(v)=P（区间含 0 + 残差验证，防超定不相容误证）。
// 对退化点曲线×普通曲线/退化到曲线×曲线交点恢复证明；vb 成功时内缩。
bool provePointOnCurve(const NURBSCurve& c, Point2 p, Interval& vb);

} // namespace nurbs

#endif // NURBS_ACCELERATED_SOLVER_HPP
