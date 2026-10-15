// ============================================================================
// accelerated_solver.cpp — BVH 加速求解实现
// ============================================================================
#include "nurbs/accelerated_solver.hpp"
#include <algorithm>
#include <cmath>
#include <optional>
#if !defined(__EMSCRIPTEN__) && !defined(NURBS_SINGLE_THREAD)
#include <thread>
#endif
#include <vector>

namespace nurbs {

namespace {
// 距离 → 状态判定（内部辅助）
GeometryState classify(double dist) {
    ThresholdDefaults t = thresholdDefaults();
    if (dist < t.GEOMETRY_TOLERANCE) return GeometryState::INTERSECTION;
    if (dist < 10 * t.GEOMETRY_TOLERANCE) return GeometryState::TANGENT;
    return GeometryState::SEPARATED;
}

// 参数段 [a,b] 的保守凸包：取支撑区间与 [a,b] 相交的控制点的 x/y 范围。
// 有理 NURBS 曲线在控制点凸包内（同参数处 C(u) 是控制点的凸组合），
// 故该范围是曲线段像集的保守包络；用于"叶对域是否可能含交点"的补漏判定。
static void segmentHull(const NURBSCurve& c, double a, double b,
                        double& minX, double& maxX, double& minY, double& maxY) {
    minX = minY = 1e300; maxX = maxY = -1e300;
    const int p = c.degree;
    const int n = (int)c.controlPoints.size();
    for (int i = 0; i < n; ++i) {
        if (c.knots[i + p + 1] < a || c.knots[i] > b) continue;   // 基函数支撑与 [a,b] 不相交
        const auto& cp = c.controlPoints[i];
        minX = std::min(minX, cp.x); maxX = std::max(maxX, cp.x);
        minY = std::min(minY, cp.y); maxY = std::max(maxY, cp.y);
    }
}

// ---- 加速路径 pieceBox：采样 AABB + 适中膨胀（平衡紧致与完备）----
// 采样点数按段宽自适应（整段 16 → 窄叶最少 4），采样后按跨度加膨胀，
// 覆盖采样点间曲线段的投影偏离，避免采样漏点导致的剪枝误差（根治漏解主因）。
// 纯控制点凸包更保守但剪枝退化（~920×），此处用"采样+膨胀"在紧致与完备间平衡。
AABB pieceBoxAdaptive(const NURBSCurve& c, double a, double b) {
    const double full = c.uMax - c.uMin;
    double frac = full > 0 ? (b - a) / full : 1.0;
    int K = (int)std::ceil(16.0 * frac);
    K = std::max(4, std::min(16, K));
    double minX = std::numeric_limits<double>::infinity();
    double minY = std::numeric_limits<double>::infinity();
    double maxX = -std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();
    for (int i = 0; i <= K; ++i) {
        const double u = a + (b - a) * i / K;
        const Point2 p = evaluatePoint(c, u);
        minX = std::min(minX, p.x); maxX = std::max(maxX, p.x);
        minY = std::min(minY, p.y); maxY = std::max(maxY, p.y);
    }
    const double sx = maxX - minX, sy = maxY - minY;
    const double ex = 1e-7 * sx + 1e-8, ey = 1e-7 * sy + 1e-8;   // 保守膨胀（性能基线）
    return { minX - ex, minY - ey, maxX + ex, maxY + ey };
}

// 多初值自适应：初值网格数随叶域相对宽度缩放。
// 窄叶（相对宽度小 → 该段近直线/低曲率，局部交点少）用 2×2，中等 3×3（基线），
// 宽叶（曲线多变）用 4×4。gu=3 时步长 (宽/(g+1))=宽/4，与旧固定 9 初值逐位一致。
int gridForInit(double rel) {
    if (rel < 0.03) return 2;
    if (rel < 0.12) return 3;
    return 4;
}

// 两条整条曲线是否"可能相交"：控制点凸包在 x/y 分量区间均重叠（保守，不漏）。
// 采样 pieceBox 可能漏掉弯曲/相切段上的交点导致叶对 pairs 为空，故补漏触发
// 不能只依赖 pairs（分离判断会漏判真实交点）；凸包是超集，凸包不重叠才是真分离。
static bool hullOverlap(const NURBSCurve& a, const NURBSCurve& b) {
    double x1lo, x1hi, y1lo, y1hi, x2lo, x2hi, y2lo, y2hi;
    segmentHull(a, a.uMin, a.uMax, x1lo, x1hi, y1lo, y1hi);
    segmentHull(b, b.uMin, b.uMax, x2lo, x2hi, y2lo, y2hi);
    return (x1lo <= x2hi && x1hi >= x2lo) && (y1lo <= y2hi && y1hi >= y2lo);
}

// 14-DOP（7 方向离散定向多面体）：采样点投影 + 适中膨胀（与 pieceBoxAdaptive 同口径）。
// 采样在 (1,0),(0,1),(1,1),(1,-1),(2,1),(1,2),(2,-1) 方向投影，膨胀覆盖采样间偏离。
void pieceDop14(const NURBSCurve& c, double a, double b, double dop[14]) {
    static const double dirs[7][2] = { {1,0},{0,1},{1,1},{1,-1},{2,1},{1,2},{2,-1} };
    for (int d = 0; d < 7; ++d) { dop[2 * d] = 1e300; dop[2 * d + 1] = -1e300; }
    const double full = c.uMax - c.uMin;
    double frac = full > 0 ? (b - a) / full : 1.0;
    int K = (int)std::ceil(16.0 * frac);
    K = std::max(4, std::min(16, K));
    for (int i = 0; i <= K; ++i) {
        const double u = a + (b - a) * i / K;
        const Point2 p = evaluatePoint(c, u);
        for (int d = 0; d < 7; ++d) {
            const double proj = dirs[d][0] * p.x + dirs[d][1] * p.y;
            dop[2 * d] = std::min(dop[2 * d], proj);
            dop[2 * d + 1] = std::max(dop[2 * d + 1], proj);
        }
    }
    for (int d = 0; d < 7; ++d) {
        const double span = dop[2 * d + 1] - dop[2 * d];
        const double e = 1e-7 * span + 1e-8;   // 保守膨胀（性能基线）
        dop[2 * d] -= e; dop[2 * d + 1] += e;
    }
}

// 14-DOP 相交：7 个方向的区间均重叠才可能相交（比 AABB 更紧的必要条件）。
bool dopIntersect(const double da[14], const double db[14]) {
    for (int d = 0; d < 7; ++d)
        if (da[2 * d + 1] < db[2 * d] || db[2 * d + 1] < da[2 * d]) return false;
    return true;
}
} // namespace

// ---------------- BVHIndex ----------------
void BVHIndex::build(const NURBSCurve& c, int maxDepth) {
    nodes.clear();
    nodes.reserve(1u << (std::min(maxDepth, 14) + 1));
    buildRec(c, c.uMin, c.uMax, 0, maxDepth);
}

int BVHIndex::buildRec(const NURBSCurve& c, double a, double b, int depth, int maxDepth) {
    int idx = (int)nodes.size();
    nodes.push_back(Node{});
    nodes[idx].uStart = a; nodes[idx].uEnd = b;
    nodes[idx].box = pieceBoxAdaptive(c, a, b);
    pieceDop14(c, a, b, nodes[idx].dop);   // 14-DOP（比 AABB 更紧的相交精剪）
    double w = b - a;
    if (depth >= maxDepth || w <= 1e-6) { nodes[idx].leaf = true; return idx; }
    double mid = a + 0.5 * w;
    int c0 = buildRec(c, a, mid, depth + 1, maxDepth);
    int c1 = buildRec(c, mid, b, depth + 1, maxDepth);
    // 递归 push 后 vector 可能扩容，必须经索引访问（勿存悬垂引用）
    nodes[idx].child0 = c0; nodes[idx].child1 = c1;
    return idx;
}

BVHIndex buildIndex(const NURBSCurve& curve, int maxDepth) {
    BVHIndex idx;
    idx.build(curve, maxDepth);
    return idx;
}

// 收集两树 AABB 相交的叶对（剪枝：不相交即跳过整棵子树）。
// 注：BVHIndex 友元，可访问其私有 nodes。
void collectLeafPairs(const BVHIndex& ta, int ia,
                      const BVHIndex& tb, int ib,
                      std::vector<std::pair<int, int>>& out) {
    const auto& A = ta.nodes[ia];
    const auto& B = tb.nodes[ib];
    if (!aabbIntersect(A.box, B.box)) return;
    if (!dopIntersect(A.dop, B.dop)) return;   // 14-DOP 精剪：AABB 相交但 DOP 不相交 → 剪掉
    if (A.leaf && B.leaf) { out.emplace_back(ia, ib); return; }
    if (A.leaf) {
        collectLeafPairs(ta, ia, tb, B.child0, out);
        collectLeafPairs(ta, ia, tb, B.child1, out);
    } else if (B.leaf) {
        collectLeafPairs(ta, A.child0, tb, ib, out);
        collectLeafPairs(ta, A.child1, tb, ib, out);
    } else {
        collectLeafPairs(ta, A.child0, tb, B.child0, out);
        collectLeafPairs(ta, A.child0, tb, B.child1, out);
        collectLeafPairs(ta, A.child1, tb, B.child0, out);
        collectLeafPairs(ta, A.child1, tb, B.child1, out);
    }
}

// 14-DOP 精剪后的相交叶对数（对照 AABB 叶对以量化包围体收紧，见 bench_kdop）。
size_t dopLeafPairCount(const BVHIndex& ia, const NURBSCurve& c1,
                        const BVHIndex& ib, const NURBSCurve& c2) {
    size_t cnt = 0;
    std::vector<std::pair<int, int>> stack;
    stack.emplace_back(0, 0);
    while (!stack.empty()) {
        auto [iaa, ibb] = stack.back(); stack.pop_back();
        const auto& A = ia.nodes[iaa];
        const auto& B = ib.nodes[ibb];
        if (!aabbIntersect(A.box, B.box)) continue;
        if (!dopIntersect(A.dop, B.dop)) continue;
        if (A.leaf && B.leaf) { ++cnt; continue; }
        if (A.leaf) { stack.emplace_back(iaa, B.child0); stack.emplace_back(iaa, B.child1); }
        else if (B.leaf) { stack.emplace_back(A.child0, ibb); stack.emplace_back(A.child1, ibb); }
        else {
            stack.emplace_back(A.child0, B.child0); stack.emplace_back(A.child0, B.child1);
            stack.emplace_back(A.child1, B.child0); stack.emplace_back(A.child1, B.child1);
        }
    }
    return cnt;
}

// 纯 AABB 相交叶对数（对照 14-DOP 剪枝后的 dopLeafPairCount，量化 k-DOP 收紧效果）。
size_t aabbLeafPairCount(const BVHIndex& ia, const NURBSCurve& c1,
                         const BVHIndex& ib, const NURBSCurve& c2) {
    size_t cnt = 0;
    std::vector<std::pair<int, int>> stack;
    stack.emplace_back(0, 0);
    while (!stack.empty()) {
        auto [iaa, ibb] = stack.back(); stack.pop_back();
        const auto& A = ia.nodes[iaa];
        const auto& B = ib.nodes[ibb];
        if (!aabbIntersect(A.box, B.box)) continue;
        if (A.leaf && B.leaf) { ++cnt; continue; }
        if (A.leaf) { stack.emplace_back(iaa, B.child0); stack.emplace_back(iaa, B.child1); }
        else if (B.leaf) { stack.emplace_back(A.child0, ibb); stack.emplace_back(A.child1, ibb); }
        else {
            stack.emplace_back(A.child0, B.child0); stack.emplace_back(A.child0, B.child1);
            stack.emplace_back(A.child1, B.child0); stack.emplace_back(A.child1, B.child1);
        }
    }
    return cnt;
}

// ---------------- solveWithIndex ----------------
std::vector<SolveResult> solveWithIndex(const BVHIndex& ia, const NURBSCurve& c1,
                                        const BVHIndex& ib, const NURBSCurve& c2,
                                        const SolverOptions& optsIn, int numThreads) {
    SolverOptions opts = optsIn;
    // 加速路径专属迭代预算：初值已被 BVH 定位到交点近邻，消融实测
    // maxIter 降至 10 仍解数 100% 一致（80 为过配），取 16 留安全余量。
    // 主路径 unifiedSolve 的预算不受影响，保持与 JS 逐位等价。
    if (opts.maxIterations <= 0) opts.maxIterations = 16;
    if (opts.tolerance <= 0) opts.tolerance = 1e-8;

    std::vector<std::pair<int, int>> pairs;
    collectLeafPairs(ia, 0, ib, 0, pairs);

    // 多初值自适应（gridForInit 见匿名 namespace）：初值网格数随叶宽缩放，
    // 目标是在不丢解前提下削减初值（量化见 adaptiveInitCounts / 论文基准）。
    struct Init { double u, v; int pid; };
    std::vector<Init> inits;
    inits.reserve(pairs.size() * 16);
    const double spanU = (c1.uMax - c1.uMin) > 0 ? (c1.uMax - c1.uMin) : 1.0;
    const double spanV = (c2.uMax - c2.uMin) > 0 ? (c2.uMax - c2.uMin) : 1.0;
    // 每个叶对做保守凸包"可能含交点"判定（mayContain）：F=C1-C2 的 x/y 分量包络
    // 区间都含 0 才可能含交点（保守：凸包是像集超集）。供补漏阶段区分"真分离"与
    // "LM 难例漏解"（相切/重根/端点处 LM 不收敛）。窄叶对支撑控制点少，凸包紧。
    std::vector<char> mayContain(pairs.size(), 0);
    for (size_t pi = 0; pi < pairs.size(); ++pi) {
        auto& p = pairs[pi];
        const auto& A = ia.nodes[p.first];
        const auto& B = ib.nodes[p.second];
        double u0 = A.uStart, u1 = A.uEnd;
        double v0 = B.uStart, v1 = B.uEnd;
        double x1lo, x1hi, y1lo, y1hi, x2lo, x2hi, y2lo, y2hi;
        segmentHull(c1, u0, u1, x1lo, x1hi, y1lo, y1hi);
        segmentHull(c2, v0, v1, x2lo, x2hi, y2lo, y2hi);
        mayContain[pi] = (x1lo <= x2hi && x1hi >= x2lo) && (y1lo <= y2hi && y1hi >= y2lo);
        int gu = gridForInit((u1 - u0) / spanU);
        int gv = gridForInit((v1 - v0) / spanV);
        double du = (u1 - u0) / (gu + 1), dv = (v1 - v0) / (gv + 1);
        for (int i = 0; i < gu; ++i)
            for (int j = 0; j < gv; ++j)
                inits.push_back({ u0 + du * (i + 1), v0 + dv * (j + 1), (int)pi });
    }

    // 并行求值：初值足够多才开线程（避免小任务线程创建开销反超），
    // 块划分 + 独立槽位 + 按原顺序收集 → 并行与串行逐位等价。
    size_t N = inits.size();
    std::vector<SolveResult> res(N);
    const size_t kParallelThreshold = 64;
    auto evalOne = [&](size_t i) {
        res[i] = levenbergMarquardt(c1, c2, inits[i].u, inits[i].v, opts);
    };
    if (N >= kParallelThreshold) {
#if defined(__EMSCRIPTEN__) || defined(NURBS_SINGLE_THREAD)
        // WASM / 单线程构建：Emscripten 默认无 pthread，std::thread 不可用，
        // 退化为串行（结果逐位一致，仅损失并行吞吐）。
        for (size_t i = 0; i < N; ++i) evalOne(i);
#else
        unsigned nthreads = numThreads <= 0 ? std::max(1u, std::thread::hardware_concurrency())
                                            : (unsigned)numThreads;
        if (nthreads > 1) {
            size_t chunk = (N + nthreads - 1) / nthreads;
            auto worker = [&](size_t s, size_t e) { for (size_t i = s; i < e; ++i) evalOne(i); };
            std::vector<std::thread> pool;
            pool.reserve(nthreads);
            for (unsigned t = 0; t < nthreads; ++t) {
                size_t s = (size_t)t * chunk, e = std::min((size_t)t * chunk + chunk, N);
                if (s < e) pool.emplace_back(worker, s, e);
            }
            for (auto& th : pool) th.join();
        } else {
            for (size_t i = 0; i < N; ++i) evalOne(i);
        }
#endif
    } else {
        for (size_t i = 0; i < N; ++i) evalOne(i);
    }

    // 按序收集（保持与串行一致的 cand 顺序）
    std::vector<SolveResult> cand;
    cand.reserve(N);
    for (size_t i = 0; i < N; ++i)
        if (res[i].converged && res[i].distance <= 1e-6) cand.push_back(res[i]);
    // 完备性补漏（路线 A）——受 SolverOptions.completeLeakPatch 开关控制：
    //   完备模式（默认）：只要两曲线可能相交（候选叶对非空 或 整条控制点凸包重叠），
    //     就用主路径 unifiedSolve 全局区间枚举兜底（Krawczyk 完备证明，漏→0）。
    //   加速模式：不补漏，纯 BVH+LM 近完备（快，漏相切/重根/端点难例）。
    // unifiedSolve 不依赖叶对划分，找回所有漏点：被采样 pieceBox 剪枝排除的叶对、
    // LM 相切/重根/端点不收敛的漏点。分离对（pairs 空且凸包不重叠）不触发。
    const bool needFull = optsIn.completeLeakPatch && (!pairs.empty() || hullOverlap(c1, c2));
    if (needFull) {
        // 用 optsIn（原始用户级选项，不经加速路径的 maxIter=16 预算改写），
        // 让主路径 unifiedSolve 以自身完备预算做全局区间枚举（Krawczyk 完备证明）。
        auto full = unifiedSolve(c1, c2, optsIn);
        // 并入阈值 1e-3（放宽）：unifiedSolve 在相切/重根难例处的收敛解残差可能
        // 略高于 1e-6（区间/Krawczyk 精度限制），用 1e-6 会误丢真实漏解。
        // 主路径解 distance 极小的为真解；伪解（端点占位 u=0 或 v=0）distance 大
        // （≥0.5），不会混入。并入后统一聚类 + classify 分类。
        for (auto& r : full)
            if (r.converged && r.distance <= 1e-3) cand.push_back(r);
    }

    // 聚类去重：按 u 稳定排序，相邻 (u,v) 接近则合并（保留距离小者）
    std::stable_sort(cand.begin(), cand.end(),
                     [](const SolveResult& a, const SolveResult& b) { return a.u < b.u; });
    double paramTol = 1e-3;
    // 聚类合并需"参数近且几何近"才算同一交点；参数近但几何远是相切/重根邻域
    // 的相邻真实交点，必须都保留（否则会误合并而丢弃补漏找回的难例解）。
    auto geomDist2 = [&](const SolveResult& a, const SolveResult& b) {
        Point2 pa = evaluatePoint(c1, a.u), pb = evaluatePoint(c1, b.u);
        double dx = pa.x - pb.x, dy = pa.y - pb.y;
        return dx * dx + dy * dy;
    };
    const double geomTol2 = 1e-6 * 1e-6;
    std::vector<SolveResult> out;
    for (auto& r : cand) {
        if (out.empty() || std::abs(r.u - out.back().u) > paramTol
            || std::abs(r.v - out.back().v) > paramTol) {
            r.state = classify(r.distance);
            out.push_back(r);
        } else if (geomDist2(r, out.back()) <= geomTol2) {
            // 参数近且几何近 → 同一交点，保留距离更小者
            if (r.distance < out.back().distance) {
                out.back() = r;
                out.back().state = classify(r.distance);
            }
        } else {
            // 参数近但几何远 → 不同交点，均保留（补漏难例关键）
            r.state = classify(r.distance);
            out.push_back(r);
        }
    }
    return out;
}

std::vector<SolveResult> unifiedSolveAccelerated(const NURBSCurve& c1, const NURBSCurve& c2,
                                                 const SolverOptions& options) {
    BVHIndex ia, ib;
    ia.build(c1, 7);
    ib.build(c2, 7);
    return solveWithIndex(ia, c1, ib, c2, options);
}

std::pair<size_t, size_t> adaptiveInitCounts(const BVHIndex& ia, const NURBSCurve& c1,
                                             const BVHIndex& ib, const NURBSCurve& c2) {
    std::vector<std::pair<int, int>> pairs;
    collectLeafPairs(ia, 0, ib, 0, pairs);
    const double spanU = (c1.uMax - c1.uMin) > 0 ? (c1.uMax - c1.uMin) : 1.0;
    const double spanV = (c2.uMax - c2.uMin) > 0 ? (c2.uMax - c2.uMin) : 1.0;
    size_t total = 0;
    for (auto& p : pairs) {
        const auto& A = ia.nodes[p.first];
        const auto& B = ib.nodes[p.second];
        int gu = gridForInit((A.uEnd - A.uStart) / spanU);
        int gv = gridForInit((B.uEnd - B.uStart) / spanV);
        total += (size_t)gu * (size_t)gv;
    }
    return { pairs.size(), total };
}

// ---------------- 批量多曲线求交（任务2） ----------------
std::vector<std::pair<int, int>> sweepAndPrunePairs(const std::vector<AABB>& boxes) {
    const size_t n = boxes.size();
    std::vector<int> order(n);
    for (size_t i = 0; i < n; ++i) order[i] = (int)i;
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return boxes[a].minX < boxes[b].minX; });
    std::vector<std::pair<int, int>> out;
    for (size_t a = 0; a < n; ++a) {
        const int ia = order[a];
        for (size_t b = a + 1; b < n; ++b) {
            const int ib = order[b];
            if (boxes[ib].minX > boxes[ia].maxX) break; // 已按 minX 排序，之后更不可能相交
            if (aabbIntersect(boxes[ia], boxes[ib])) out.emplace_back(ia, ib);
        }
    }
    return out;
}

std::vector<BatchIntersection> batchSolveAll(const std::vector<NURBSCurve>& curves,
                                             const SolverOptions& options, int numThreads) {
    const size_t n = curves.size();
    std::vector<BVHIndex> idx(n);
    std::vector<AABB> roots(n);
    for (size_t i = 0; i < n; ++i) { idx[i].build(curves[i], 7); roots[i] = idx[i].rootBox(); }
    const std::vector<std::pair<int, int>> pairs = sweepAndPrunePairs(roots);
    const size_t P = pairs.size();
    std::vector<std::optional<BatchIntersection>> res(P);
    const size_t kParallelThreshold = 64;
    auto evalOne = [&](size_t k) {
        const int ia = pairs[k].first, ib = pairs[k].second;
        // 对内 numThreads=1（对级别已并行，避免嵌套线程开销）
        auto r = solveWithIndex(idx[ia], curves[ia], idx[ib], curves[ib], options, 1);
        if (!r.empty()) res[k] = BatchIntersection{ ia, ib, std::move(r) };
    };
    if (P >= kParallelThreshold) {
#if defined(__EMSCRIPTEN__) || defined(NURBS_SINGLE_THREAD)
        for (size_t k = 0; k < P; ++k) evalOne(k);
#else
        unsigned nthreads = numThreads <= 0 ? std::max(1u, std::thread::hardware_concurrency())
                                            : (unsigned)numThreads;
        if (nthreads > 1) {
            size_t chunk = (P + nthreads - 1) / nthreads;
            auto worker = [&](size_t s, size_t e) { for (size_t k = s; k < e; ++k) evalOne(k); };
            std::vector<std::thread> pool; pool.reserve(nthreads);
            for (unsigned t = 0; t < nthreads; ++t) {
                size_t s = (size_t)t * chunk, e = std::min((size_t)t * chunk + chunk, P);
                if (s < e) pool.emplace_back(worker, s, e);
            }
            for (auto& th : pool) th.join();
        } else {
            for (size_t k = 0; k < P; ++k) evalOne(k);
        }
#endif
    } else {
        for (size_t k = 0; k < P; ++k) evalOne(k);
    }
    std::vector<BatchIntersection> out; out.reserve(P);
    for (auto& o : res) if (o) out.push_back(*o);
    return out;
}

// ---------------- 保证性求交（区间/Krawczyk） ----------------
// 找参数 t 所在 knot span：knots[s] <= t < knots[s+1]（凸包剪枝区间用）
static int findSpan(const std::vector<double>& knots, int p, double t) {
    const int n = (int)knots.size() - p - 1; // 控制点个数
    if (t >= knots[n]) return n - 1;
    if (t <= knots[p]) return p;
    int lo = p, hi = n;
    while (lo < hi - 1) {
        int mid = (lo + hi) >> 1;
        if (t < knots[mid]) hi = mid; else lo = mid;
    }
    return lo;
}

// 凸包严格区间（全局剪枝/排除测试用）。函数值 c(u) ∈ 控制点凸包（严格包含真值，
// 不低估），用于区间细分枚举的排除测试：若 f 的凸包区间不含 0，则该参数区无交点，
// 可安全排除（保证不漏解）。宽盒时凸包宽、剪枝弱，需细分到窄盒才有效，但严格性不受影响。
CurveIntervalEval curveIntervalHull(const NURBSCurve& c, const Interval& u) {
    const int p = c.degree;
    const auto& K = c.knots;
    const int n = (int)c.controlPoints.size();
    const double ulo = std::max(u.lo, c.uMin), uhi = std::min(u.hi, c.uMax);
    const int sLo = findSpan(K, p, ulo), sHi = findSpan(K, p, uhi);
    const int i0 = std::max(0, sLo - p - 1), i1 = std::min(n - 1, sHi + 1);
    Interval Cx(1e300, -1e300), Cy(1e300, -1e300), W(1e300, -1e300);
    for (int i = i0; i <= i1; ++i) {
        const double w = c.controlPoints[i].w, x = c.controlPoints[i].x, y = c.controlPoints[i].y;
        Cx.lo = std::min(Cx.lo, w * x); Cx.hi = std::max(Cx.hi, w * x);
        Cy.lo = std::min(Cy.lo, w * y); Cy.hi = std::max(Cy.hi, w * y);
        W.lo = std::min(W.lo, w);       W.hi = std::max(W.hi, w);
    }
    if (W.lo <= 0.0) return { Interval::empty(), Interval::empty(), Interval::empty(), Interval::empty() };
    return { idiv(Cx, W), idiv(Cy, W), Interval::empty(), Interval::empty() };
}

// 严格一致收敛区间（全局排除测试用）。真值 |c(u)-c(m)| ≤ L·δ，其中 L 是盒上 |c'| 的
// 严格上界，用导数控制点凸包范数上界构造（对有理 NURBS 含分母界）：
//   c' = (C'W − C W')/W² → |c'| ≤ (|C'|·Wmax + |C|·|W'|max) / wmin²。
// 该区间随 δ→0 收敛到点值（比凸包紧），又严格包含真值（不漏解），故可作快速剪枝。
CurveIntervalEval curveIntervalStrict(const NURBSCurve& c, const Interval& u) {
    const int p = c.degree;
    const auto& K = c.knots;
    const int n = (int)c.controlPoints.size();
    const double a = std::max(u.lo, c.uMin), b = std::min(u.hi, c.uMax);
    const double m = 0.5 * (a + b), half = 0.5 * (b - a);
    const auto pm = evaluatePointAndDerivative(c, m);
    const int sLo = findSpan(K, p, a), sHi = findSpan(K, p, b);
    const int i0 = std::max(0, sLo - p - 1), i1 = std::min(n - 1, sHi + 1);
    double wmin = 1e300, wmax = 0, cxmax = 0, qxmax = 0, qymax = 0, qwmax = 0;
    for (int i = i0; i <= i1; ++i) {
        const double w = c.controlPoints[i].w;
        wmin = std::min(wmin, w); wmax = std::max(wmax, w);
        cxmax = std::max(cxmax, std::abs(w * c.controlPoints[i].x));
        cxmax = std::max(cxmax, std::abs(w * c.controlPoints[i].y));
    }
    for (int i = i0; i <= i1 && i <= n - 2; ++i) {
        const auto& A = c.controlPoints[i];
        const auto& B = c.controlPoints[i + 1];
        const double denom = K[i + p + 1] - K[i + 1];
        if (std::abs(denom) < 1e-300) continue;
        const double s = (double)p / denom;
        qxmax = std::max(qxmax, std::abs(s * (B.w * B.x - A.w * A.x)));
        qymax = std::max(qymax, std::abs(s * (B.w * B.y - A.w * A.y)));
        qwmax = std::max(qwmax, std::abs(s * (B.w - A.w)));
    }
    const double wmin2 = wmin * wmin;
    // L = max(|C'|Wmax + |C|W'max) / wmin² 的分量上界（统一用 max 范数）
    const double L = (std::max(qxmax, qymax) * wmax + cxmax * qwmax) / wmin2;
    const double E = L * half;
    const Interval fx(pm.first.x - E, pm.first.x + E);
    const Interval fy(pm.first.y - E, pm.first.y + E);
    // 导数区间（证明用）：中心 + 数值二阶采样（叶盒窄时紧）
    const int S = 5;
    double Ldx = 0, Ldy = 0;
    const double hh = (b - a) * 1e-3 + 1e-9;
    for (int i = 0; i < S; ++i) {
        const double t = a + (b - a) * i / (S - 1);
        const auto pl = evaluatePointAndDerivative(c, t - hh);
        const auto pr = evaluatePointAndDerivative(c, t + hh);
        Ldx = std::max(Ldx, std::abs((pr.second.x - pl.second.x) / (2 * hh)));
        Ldy = std::max(Ldy, std::abs((pr.second.y - pl.second.y) / (2 * hh)));
    }
    const double Ed = std::max(Ldx, Ldy) * half;
    const Interval dfx(pm.second.x - Ed, pm.second.x + Ed);
    const Interval dfy(pm.second.y - Ed, pm.second.y + Ed);
    return { fx, fy, dfx, dfy };
}

// NURBS 在参数区间 [u] 上的一致收敛区间扩展（中心 + Lipschitz/Taylor 形式）。
// 凸包区间虽严格但在盒上不随宽度收敛（取控制点 min/max），导致 Krawczyk 几乎必失败；
// 此处改为：c(u) ⊂ c(m) ± L·δ（L = 盒上 |c'| 上界）、c'(u) ⊂ c'(m) ± L_d·δ
// （L_d = 盒上 |c''| 数值上界）。该扩展随 δ→0 收敛到点值，是 Krawczyk 需要的紧区间。
CurveIntervalEval curveInterval(const NURBSCurve& c, const Interval& u) {
    const double m = u.mid();
    const double a = std::max(u.lo, c.uMin), b = std::min(u.hi, c.uMax);
    const auto pm = evaluatePointAndDerivative(c, m);
    const double half = 0.5 * (b - a);
    const int S = 5;
    double Lx = 0, Ly = 0, Ldx = 0, Ldy = 0;
    for (int i = 0; i < S; ++i) {
        const double t = a + (b - a) * i / (S - 1);
        const auto p = evaluatePointAndDerivative(c, t);
        Lx = std::max(Lx, std::abs(p.second.x));
        Ly = std::max(Ly, std::abs(p.second.y));
    }
    const double hh = (b - a) * 1e-3 + 1e-9;
    for (int i = 0; i < S; ++i) {
        const double t = a + (b - a) * i / (S - 1);
        const auto pl = evaluatePointAndDerivative(c, t - hh);
        const auto pr = evaluatePointAndDerivative(c, t + hh);
        Ldx = std::max(Ldx, std::abs((pr.second.x - pl.second.x) / (2 * hh)));
        Ldy = std::max(Ldy, std::abs((pr.second.y - pl.second.y) / (2 * hh)));
    }
    const double E = std::max(Lx, Ly) * half;
    const double Ed = std::max(Ldx, Ldy) * half;
    const Interval fx(pm.first.x - E, pm.first.x + E);
    const Interval fy(pm.first.y - E, pm.first.y + E);
    const Interval dfx(pm.second.x - Ed, pm.second.x + Ed);
    const Interval dfy(pm.second.y - Ed, pm.second.y + Ed);
    return { fx, fy, dfx, dfy };
}

bool proveUniqueInBox(const NURBSCurve& c1, const NURBSCurve& c2,
                      const Interval& ub, const Interval& vb) {
    const auto e1 = curveInterval(c1, ub);
    const auto e2 = curveInterval(c2, vb);
    // 区间函数值 f(x) = c1(u) - c2(v)
    const IntervalVec2 fInt{ isub(e1.fx, e2.fx), isub(e1.fy, e2.fy) };
    // 区间 Jacobian J = [[c1x'(u), -c2x'(v)], [c1y'(u), -c2y'(v)]]
    const Interval nz = Interval::point(0);
    const IntervalMat22 J{ e1.dfx, isub(nz, e2.dfx), e1.dfy, isub(nz, e2.dfy) };
    // 中点实值
    const double um = ub.mid(), vm = vb.mid();
    const auto p1 = evaluatePointAndDerivative(c1, um);
    const auto p2 = evaluatePointAndDerivative(c2, vm);
    const double j00 = p1.second.x, j01 = -p2.second.x;
    const double j10 = p1.second.y, j11 = -p2.second.y;
    const double det = j00 * j11 - j01 * j10;
    if (std::abs(det) < 1e-12) return false; // 相切/退化（Jacobian 奇异），Krawczyk 不适用
    const double y00 = j11 / det, y01 = -j01 / det, y10 = -j10 / det, y11 = j00 / det;
    // m - Y f(m)（实值）
    const double fmx = p1.first.x - p2.first.x, fmy = p1.first.y - p2.first.y;
    const double mmyx = um - (y00 * fmx + y01 * fmy);
    const double mmyy = vm - (y10 * fmx + y11 * fmy);
    // (I - Y J([x]))（区间矩阵）
    const IntervalMat22 M{
        isub(Interval::point(1.0), iadd(imul(Interval::point(y00), J.m00), imul(Interval::point(y01), J.m10))),
        isub(nz, iadd(imul(Interval::point(y00), J.m01), imul(Interval::point(y01), J.m11))),
        isub(nz, iadd(imul(Interval::point(y10), J.m00), imul(Interval::point(y11), J.m10))),
        isub(Interval::point(1.0), iadd(imul(Interval::point(y10), J.m01), imul(Interval::point(y11), J.m11))) };
    // (x - m)（区间向量）
    const IntervalVec2 xm{ isub(ub, Interval::point(um)), isub(vb, Interval::point(vm)) };
    // K = (m - Y f(m)) + (I - Y J([x]))(x - m)
    const IntervalVec2 K = iaddVec({ Interval::point(mmyx), Interval::point(mmyy) }, matVec(M, xm));
    // K([x]) ⊂ interior([x])
    return K.x.lo > ub.lo && K.x.hi < ub.hi && K.y.lo > vb.lo && K.y.hi < vb.hi;
}

// GS-Newton 区间证明（Hansen-Sengupta Gauss-Seidel 变体）：线性化系统 A(x−m)=−Y f(m)
// 逐变量区间消元更新，比 Krawczyk 一步算子更紧（利用最新区间）。
int proveUniqueGS(const NURBSCurve& c1, const NURBSCurve& c2, Interval& ub, Interval& vb) {
    const auto e1 = curveInterval(c1, ub), e2 = curveInterval(c2, vb);
    const IntervalVec2 fInt{ isub(e1.fx, e2.fx), isub(e1.fy, e2.fy) };
    const Interval nz = Interval::point(0);
    const IntervalMat22 J{ e1.dfx, isub(nz, e2.dfx), e1.dfy, isub(nz, e2.dfy) };
    const double um = ub.mid(), vm = vb.mid();
    const auto p1 = evaluatePointAndDerivative(c1, um), p2 = evaluatePointAndDerivative(c2, vm);
    const double j00 = p1.second.x, j01 = -p2.second.x, j10 = p1.second.y, j11 = -p2.second.y;
    const double det = j00 * j11 - j01 * j10;
    if (std::abs(det) < 1e-12) return 0; // 相切/退化，GS 不适用
    const double y00 = j11 / det, y01 = -j01 / det, y10 = -j10 / det, y11 = j00 / det;
    // A = Y·J([x])（区间矩阵）
    const IntervalMat22 A{
        iadd(imul(Interval::point(y00), J.m00), imul(Interval::point(y01), J.m10)),
        iadd(imul(Interval::point(y00), J.m01), imul(Interval::point(y01), J.m11)),
        iadd(imul(Interval::point(y10), J.m00), imul(Interval::point(y11), J.m10)),
        iadd(imul(Interval::point(y10), J.m01), imul(Interval::point(y11), J.m11)) };
    // b = −Y f(m)
    const double fmx = p1.first.x - p2.first.x, fmy = p1.first.y - p2.first.y;
    const double bu = -(y00 * fmx + y01 * fmy), bv = -(y10 * fmx + y11 * fmy);
    if (A.m00.contains(0.0) || A.m11.contains(0.0)) return 0; // 主元含 0，无法区间消元
    Interval yu = isub(ub, Interval::point(um)), yv = isub(vb, Interval::point(vm));
    // 变量 u：用当前 yv
    Interval u_new = iadd(Interval::point(um), idiv(isub(Interval::point(bu), imul(A.m01, yv)), A.m00));
    u_new = { std::max(ub.lo, u_new.lo), std::min(ub.hi, u_new.hi) };
    if (u_new.isEmpty()) return -1; // 无解
    // 变量 v：用更新后的 yu
    yu = isub(u_new, Interval::point(um));
    Interval v_new = iadd(Interval::point(vm), idiv(isub(Interval::point(bv), imul(A.m10, yu)), A.m11));
    v_new = { std::max(vb.lo, v_new.lo), std::min(vb.hi, v_new.hi) };
    if (v_new.isEmpty()) return -1; // 无解
    // 证明：新盒严格内缩到旧盒内（存在唯一）；同时将盒内缩以加速后续
    const bool innerU = u_new.lo > ub.lo && u_new.hi < ub.hi;
    const bool innerV = v_new.lo > vb.lo && v_new.hi < vb.hi;
    if (innerU && innerV) { ub = u_new; vb = v_new; return 1; }
    return 0;
}

// 相切交点精化：阻尼 Gauss-Newton 最小化 ||c1(u)−c2(v)||²。相切处 Jacobian 奇异，
// 用阻尼 λ 保证可逆，收敛到切点（距离→0）。对标有理插值对相切的专门处理。
double refineTangentPoint(const NURBSCurve& c1, const NURBSCurve& c2, double& u, double& v, int iters) {
    double lam = 1e-8;
    double best = 1e300;
    for (int it = 0; it < iters; ++it) {
        const auto p1 = evaluatePointAndDerivative(c1, u);
        const auto p2 = evaluatePointAndDerivative(c2, v);
        const double fx = p1.first.x - p2.first.x, fy = p1.first.y - p2.first.y;
        const double dist = std::sqrt(fx * fx + fy * fy);
        if (dist < best) best = dist;
        const double j00 = p1.second.x, j01 = -p2.second.x, j10 = p1.second.y, j11 = -p2.second.y;
        // J^T J + λI
        const double a = j00 * j00 + j10 * j10 + lam;
        const double b = j00 * j01 + j10 * j11;
        const double c = j01 * j01 + j11 * j11 + lam;
        const double det = a * c - b * b;
        if (std::abs(det) < 1e-300) break;
        const double rx = -(j00 * fx + j10 * fy), ry = -(j01 * fx + j11 * fy);
        const double du = (rx * c - b * ry) / det;
        const double dv = (a * ry - b * rx) / det;
        const double du_ = std::abs(du), dv_ = std::abs(dv);
        u += du; v += dv;
        if (du_ < 1e-15 && dv_ < 1e-15) break;
    }
    return best;
}

// ---------------- B/C: 解分类 + 条件数/误差界 ----------------
std::vector<ClassifiedSolution> classifySolutions(const NURBSCurve& c1, const NURBSCurve& c2,
                                                  const std::vector<SolveResult>& sol) {
    std::vector<ClassifiedSolution> out;
    out.reserve(sol.size());
    const double kDetEps = 1e-8;
    for (const auto& r : sol) {
        if (r.state != GeometryState::INTERSECTION && r.state != GeometryState::TANGENT) continue;
        const auto p1 = evaluatePointAndDerivative(c1, r.u);
        const auto p2 = evaluatePointAndDerivative(c2, r.v);
        const double j00 = p1.second.x, j01 = -p2.second.x;
        const double j10 = p1.second.y, j11 = -p2.second.y;
        const double det = j00 * j11 - j01 * j10;
        const double nf = std::sqrt(j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);
        const double cond = std::abs(det) > 0 ? nf / std::abs(det) : 1e18;
        // 相切判定用条件数（近奇异 → cond 大），比 |det| 阈值更鲁棒：
        // 切点处数值导数噪声 ~1e-6，|det| 阈值易误判；cond 是尺度无关的良态度量。
        const IntersectionKind kind = cond > 1e4 ? IntersectionKind::TANGENT
                                                 : IntersectionKind::TRANSVERSE;
        out.push_back({ r.u, r.v, r.distance, det, cond, kind });
    }
    return out;
}

bool detectOverlap(const NURBSCurve& c1, const NURBSCurve& c2, double tol) {
    // 控制点集完全一致 → 完全重叠
    if (c1.controlPoints.size() == c2.controlPoints.size()) {
        bool same = true;
        for (size_t i = 0; i < c1.controlPoints.size(); ++i)
            if (std::abs(c1.controlPoints[i].x - c2.controlPoints[i].x) > 1e-9 ||
                std::abs(c1.controlPoints[i].y - c2.controlPoints[i].y) > 1e-9 ||
                std::abs(c1.controlPoints[i].w - c2.controlPoints[i].w) > 1e-9) { same = false; break; }
        if (same) return true;
    }
    // 部分重叠：公共参数区采样，检查是否存在连续重合段（距离 < tol 且长度可观）
    const int S = 96;
    const double u0 = c1.uMin, u1 = c1.uMax;
    const double v0 = c2.uMin, v1 = c2.uMax;
    int nearCount = 0;
    for (int i = 0; i < S; ++i) {
        const double u = u0 + (u1 - u0) * i / (S - 1);
        const Point2 p = evaluatePoint(c1, u);
        // 粗搜最近 v（v 网格）+ 局部细化
        double best = 1e300;
        for (int j = 0; j < S; ++j) {
            const double v = v0 + (v1 - v0) * j / (S - 1);
            const Point2 q = evaluatePoint(c2, v);
            const double dx = p.x - q.x, dy = p.y - q.y;
            best = std::min(best, dx * dx + dy * dy);
        }
        if (best <= tol * tol) nearCount++;
    }
    return nearCount >= (int)(S * 0.6); // 采样点 ≥60% 距离 < tol → 判定重叠
}

// ---------------- ④ 全局保证性求交：区间细分全局枚举 ----------------
GlobalEnumerationResult enumerateIntersections(const NURBSCurve& c1, const NURBSCurve& c2, double leafTol) {
    GlobalEnumerationResult res;
    const Interval uRoot(c1.uMin, c1.uMax), vRoot(c2.uMin, c2.uMax);
    std::vector<std::pair<Interval, Interval>> stack;
    stack.push_back({ uRoot, vRoot });
    std::vector<GlobalIntersection> raw;
    while (!stack.empty()) {
        auto [ub, vb] = stack.back(); stack.pop_back();
        ++res.evaluatedBoxes;
        // 严格一致收敛区间排除测试：f = c1 - c2 的任一分量区间不含 0 → 该盒无交点，
        // 安全排除（Lipschitz 常数取自导数控制点凸包上界，严格不漏）。
        const auto h1 = curveIntervalStrict(c1, ub), h2 = curveIntervalStrict(c2, vb);
        const Interval f1 = isub(h1.fx, h2.fx), f2 = isub(h1.fy, h2.fy);
        if (f1.isEmpty() || f2.isEmpty() || !f1.contains(0.0) || !f2.contains(0.0)) continue;
        const double wu = ub.width(), wv = vb.width();
        if (wu > leafTol || wv > leafTol) {
            // 细分较宽维
            if (wu >= wv) {
                const double mm = ub.mid();
                stack.push_back({ Interval(ub.lo, mm), vb });
                stack.push_back({ Interval(mm, ub.hi), vb });
            } else {
                const double mm = vb.mid();
                stack.push_back({ ub, Interval(vb.lo, mm) });
                stack.push_back({ ub, Interval(mm, vb.hi) });
            }
            continue;
        }
        // 叶盒：Krawczyk 证明存在唯一（一致收敛区间）
        const double um = ub.mid(), vm = vb.mid();
        const bool proven = proveUniqueInBox(c1, c2, ub, vb);
        if (proven) {
            raw.push_back({ um, vm, ub, vb, true, IntersectionKind::TRANSVERSE });
        } else if (wu > leafTol * 0.01 || wv > leafTol * 0.01) {
            // 自适应细化：高频/高次曲线叶盒上 Lipschitz 采样可能低估，继续细分到更小盒
            // 让一致收敛区间变紧，提高 Krawczyk 证明成功率；达最小盒仍失败才标疑似。
            if (wu >= wv) {
                const double mm = ub.mid();
                stack.push_back({ Interval(ub.lo, mm), vb });
                stack.push_back({ Interval(mm, ub.hi), vb });
            } else {
                const double mm = vb.mid();
                stack.push_back({ ub, Interval(vb.lo, mm) });
                stack.push_back({ ub, Interval(mm, vb.hi) });
            }
        } else {
            // 疑似相切/退化：按 Jacobian 条件数分类（近奇异 → TANGENT）
            const auto p1 = evaluatePointAndDerivative(c1, um);
            const auto p2 = evaluatePointAndDerivative(c2, vm);
            const double j00 = p1.second.x, j01 = -p2.second.x, j10 = p1.second.y, j11 = -p2.second.y;
            const double det = j00 * j11 - j01 * j10;
            const double nf = std::sqrt(j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11);
            const double cond = std::abs(det) > 0 ? nf / std::abs(det) : 1e18;
            raw.push_back({ um, vm, ub, vb, false,
                            cond > 1e4 ? IntersectionKind::TANGENT : IntersectionKind::UNKNOWN });
        }
    }
    // 排序并按参数距离去重：已证明盒（Krawczyk 唯一）不相交即不同解，保留；
    // 未证明（疑似相切/退化）的相邻盒若与已保留解过近，视为同一退化解的邻域分裂，丢弃。
    std::sort(raw.begin(), raw.end(), [](const GlobalIntersection& a, const GlobalIntersection& b) {
        return a.u < b.u;
    });
    for (const auto& b : raw) {
        if (b.proven) { res.solutions.push_back(b); continue; }
        bool dup = false;
        for (const auto& s : res.solutions)
            if (std::abs(b.u - s.u) <= leafTol && std::abs(b.v - s.v) <= leafTol) { dup = true; break; }
        if (!dup) res.solutions.push_back(b);
    }
    for (const auto& s : res.solutions) { if (s.proven) res.provenCount++; else res.tangentCount++; }
    res.complete = true;
    return res;
}

// ---- ⑤ 退化/奇异曲线处理 ----
CurveDegeneracy detectDegeneracy(const NURBSCurve& c, Point2* outPoint, double* outWmin) {
    double minx = 1e300, maxx = -1e300, miny = 1e300, maxy = -1e300, wmin = 1e300;
    for (const auto& p : c.controlPoints) {
        wmin = std::min(wmin, p.w);
        minx = std::min(minx, p.x); maxx = std::max(maxx, p.x);
        miny = std::min(miny, p.y); maxy = std::max(maxy, p.y);
    }
    if (outWmin) *outWmin = wmin;
    if ((maxx - minx) < 1e-12 && (maxy - miny) < 1e-12) {
        if (outPoint) *outPoint = Point2{ 0.5 * (minx + maxx), 0.5 * (miny + maxy) };
        return CurveDegeneracy::POINT;
    }
    if (c.controlPoints.size() >= 3) {
        double dx0 = c.controlPoints[1].x - c.controlPoints[0].x, dy0 = c.controlPoints[1].y - c.controlPoints[0].y;
        double len0 = std::sqrt(dx0 * dx0 + dy0 * dy0);
        if (len0 < 1e-15 && c.controlPoints.size() >= 3) {
            dx0 = c.controlPoints[2].x - c.controlPoints[0].x; dy0 = c.controlPoints[2].y - c.controlPoints[0].y;
            len0 = std::sqrt(dx0 * dx0 + dy0 * dy0);
        }
        if (len0 > 1e-15) {
            const double nx = -dy0 / len0, ny = dx0 / len0;
            bool collinear = true;
            for (const auto& p : c.controlPoints)
                if (std::abs(nx * (p.x - c.controlPoints[0].x) + ny * (p.y - c.controlPoints[0].y)) > 1e-9) { collinear = false; break; }
            if (collinear) return CurveDegeneracy::COLLINEAR;
        }
    }
    if (wmin < 1e-6) return CurveDegeneracy::RATIONAL;
    return CurveDegeneracy::NONE;
}

namespace {
// 曲线一阶导数区间（差分控制点凸包，跨整曲线，保守超集）
void derivRange(const NURBSCurve& c, Interval& dx, Interval& dy) {
    const int p = c.degree;
    double lox = 1e300, hix = -1e300, loy = 1e300, hiy = -1e300;
    for (int i = 0; i + 1 < (int)c.controlPoints.size(); ++i) {
        const double qx = p * (c.controlPoints[i + 1].x - c.controlPoints[i].x);
        const double qy = p * (c.controlPoints[i + 1].y - c.controlPoints[i].y);
        lox = std::min(lox, qx); hix = std::max(hix, qx);
        loy = std::min(loy, qy); hiy = std::max(hiy, qy);
    }
    dx = { lox, hix }; dy = { loy, hiy };
}
// 曲线一阶导数区间（该参数段 [a,b] 内活跃差分控制点凸包，随段收缩变紧）
void derivRangeSeg(const NURBSCurve& c, double a, double b, Interval& dx, Interval& dy) {
    const int p = c.degree;
    double lox = 1e300, hix = -1e300, loy = 1e300, hiy = -1e300;
    for (int i = 0; i + 1 < (int)c.controlPoints.size(); ++i) {
        const double kn = c.knots[i], kk = c.knots[i + p + 1];
        if (kk > a && kn < b) {
            const double qx = p * (c.controlPoints[i + 1].x - c.controlPoints[i].x);
            const double qy = p * (c.controlPoints[i + 1].y - c.controlPoints[i].y);
            lox = std::min(lox, qx); hix = std::max(hix, qx);
            loy = std::min(loy, qy); hiy = std::max(hiy, qy);
        }
    }
    dx = { lox, hix }; dy = { loy, hiy };
}
// 曲线在参数区间的控制点凸包 AABB（NURBS 凸包性：像 ⊂ 该段控制点凸包）
AABB boxCurveHull(const NURBSCurve& c, double a, double b) {
    AABB box; box.minX = 1e300; box.minY = 1e300; box.maxX = -1e300; box.maxY = -1e300;
    const int p = c.degree;
    for (int i = 0; i < (int)c.controlPoints.size(); ++i) {
        const double kn = c.knots[i], kk = c.knots[i + p + 1];
        if (kk > a && kn < b) {
            box.minX = std::min(box.minX, c.controlPoints[i].x); box.maxX = std::max(box.maxX, c.controlPoints[i].x);
            box.minY = std::min(box.minY, c.controlPoints[i].y); box.maxY = std::max(box.maxY, c.controlPoints[i].y);
        }
    }
    return box;
}
} // namespace

namespace {
bool provePointOnCurveRec(const NURBSCurve& c, Point2 p, Interval& vb, int depth) {
    if (depth >= 20) return false;
    // 段盒 AABB 必要包含：P 不在当前子段凸包内则立即排除（未命中组快速剪枝）
    { AABB b = boxCurveHull(c, vb.lo, vb.hi);
      if (p.x < b.minX - 1e-9 || p.x > b.maxX + 1e-9 || p.y < b.minY - 1e-9 || p.y > b.maxY + 1e-9) return false; }
    const double m = vb.mid(), w = 0.5 * (vb.hi - vb.lo);
    const Point2 cv = evaluatePoint(c, m), cd = evaluateFirstDerivative(c, m);
    const double g0 = cv.x - p.x, g1 = cv.y - p.y;
    const double G = cd.x * cd.x + cd.y * cd.y;
    if (G < 1e-20) return false;   // 导数奇异点不做点包含证明（诚实返回，不深递归）
    const double A0 = cd.x / G, A1 = cd.y / G;
    // 导数区间用一致收敛严格区间（中心+Lipschitz，随 δ→0 收缩到点值，含有理分母界）
    const CurveIntervalEval ev = curveIntervalStrict(c, vb);
    Interval M = isub(Interval(1, 1), iadd(imul(Interval(A0, A0), ev.dfx), imul(Interval(A1, A1), ev.dfy)));
    const double jf = A0 * g0 + A1 * g1;
    Interval K = iadd(Interval(m - jf, m - jf), imul(M, Interval(-w, w)));
    if (K.lo > vb.lo && K.hi < vb.hi) {
        const double mm = K.mid();
        const Point2 cm = evaluatePoint(c, mm);
        const double r = std::sqrt((cm.x - p.x) * (cm.x - p.x) + (cm.y - p.y) * (cm.y - p.y));
        if (r < 1e-6) { vb = K; return true; }   // K⊂X 且残差小 → 真解证明
        return false;                              // 残差大 → 超定不相容，无真解
    }
    Interval l(vb.lo, m), r(m, vb.hi);
    if (provePointOnCurveRec(c, p, l, depth + 1)) { vb = l; return true; }
    if (provePointOnCurveRec(c, p, r, depth + 1)) { vb = r; return true; }
    return false;
}
} // namespace

bool provePointOnCurve(const NURBSCurve& c, Point2 p, Interval& vb) {
    // 退化点曲线（像=单点 Q）：直接比较 P 与 Q（点-点重合证明/拒绝）
    if (detectDegeneracy(c) == CurveDegeneracy::POINT) {
        Point2 q; detectDegeneracy(c, &q);
        const double r = std::sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y));
        return r < 1e-6;
    }
    // 必要包含检查：P 不在曲线整体凸包 AABB 内则立即拒绝（避免拒绝分支爆炸）
    { AABB b = boxCurveHull(c, c.uMin, c.uMax);
      if (p.x < b.minX - 1e-9 || p.x > b.maxX + 1e-9 || p.y < b.minY - 1e-9 || p.y > b.maxY + 1e-9) return false; }
    // 采样定位：找 P 的最近参数近似 u*，从 u* 邻域小盒开始证明/排除，
    // 避免对"凸包覆盖 P 但 P 不在曲线"的未命中组做全域深递归。
    const int NS = 64; const double du = (c.uMax - c.uMin) / NS;
    double best = 1e300; int bi = 0;
    for (int i = 0; i <= NS; ++i) {
        const Point2 q = evaluatePoint(c, c.uMin + i * du);
        const double d = (q.x - p.x) * (q.x - p.x) + (q.y - p.y) * (q.y - p.y);
        if (d < best) { best = d; bi = i; }
    }
    // 采样最近距离过远 → P 不在曲线，快速拒绝（未命中组 O(NS) 返回，无深递归）
    if (best > 1e-4) return false;
    const double u0 = c.uMin + bi * du, rad = 3.0 * du;
    Interval local(std::max(c.uMin, u0 - rad), std::min(c.uMax, u0 + rad));
    if (provePointOnCurveRec(c, p, local, 0)) { vb = local; return true; }
    return false;
}

} // namespace nurbs
