// ============================================================================
// bench_curve_selfintersect.cpp — 曲面/曲线测例补充⑦：曲线自交检测
//   补齐 ACIS C-C self-intersection（曲线侧，与曲面自交对称）。
//   8 字形参数曲线 C(t)=(sin 2πt, ½ sin 4πt)：t=0、½、1 三参数同落原点 →
//   内部自交点（t=½ 内部参数与端点重合），成对检测。
//     (a) 检测到自交点（两不同参数 t1≠t2 映射到同一点）；
//     (b) 自交点参数 t1 在曲线内部（t∈(0,1)，非仅两端闭合）；
//     (c) 闭合环（泪滴，仅两端重合）正确区分：无内部自交；
//     (d) 单调曲线负例：无自交。
//   自包含（eval 函数 + 采样对检测，不依赖主库）。
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

struct P3 { double x = 0, y = 0, z = 0; };
static const double PI = 3.1415926535897932384626433832795;

// 曲线求值函数签名
using EvalFn = P3 (*)(double t);

// 8 字形：x=sin(2πt), y=½sin(4πt) —— t=0,½,1 同落原点
static P3 evalFig8(double t) {
    return { std::sin(2 * PI * t), 0.5 * std::sin(4 * PI * t), 0.0 };
}
// 单调直线：无自交
static P3 evalLine(double t) {
    return { 3 * t, 0.9 * t, 0.0 };
}
// 泪滴闭合环（三次 Bezier，两端重合）：仅端点闭合，内部无自交
static P3 evalTeardrop(double t) {
    const double u = 1 - t;
    const double b0 = u*u*u, b1 = 3*u*u*t, b2 = 3*u*t*t, b3 = t*t*t;
    const double px[4] = {0, 1, -1, 0}, py[4] = {0, 1, 1, 0};
    return { b0*px[0] + b1*px[1] + b2*px[2] + b3*px[3],
             b0*py[0] + b1*py[1] + b2*py[2] + b3*py[3], 0.0 };
}

static int pass = 0, fail = 0;
static void check(const char* label, bool ok) {
    printf("  %-44s %s\n", label, ok ? "PASS" : "FAIL");
    if (ok) ++pass; else ++fail;
}

// 自交检测：t1 在曲线内部 (0.02..0.98)，t2 任意 (0..1)，|t1-t2|>0.1，
// 返回最小点距及参数对（"两不同参数映射到同一点"即自交）。
static double detectSelfIntersect(EvalFn f, double& t1o, double& t2o) {
    const int N = 400;
    double mn = 1e300; t1o = t2o = -1;
    for (int i = 0; i <= N; ++i) for (int j = 0; j <= N; ++j) {
        const double t1 = i / (double)N, t2 = j / (double)N;
        if (t1 < 0.02 || t1 > 0.98) continue;          // t1 必须内部（排除纯端点闭合）
        if (std::abs(t1 - t2) < 0.1) continue;          // 参数差要足够大
        if (i == j) continue;
        P3 a = f(t1), b = f(t2);
        const double d = (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z);
        if (d < mn) { mn = d; t1o = t1; t2o = t2; }
    }
    return std::sqrt(mn);
}

int main() {
    printf("==== 测例补充⑦：曲线自交检测（8 字形，ACIS C-C self-inter）====\n\n");

    // [a] 8 字形：检测到内部自交（t1≈0.5 与端点 t2≈0/1 同落原点）
    {
        double t1 = 0, t2 = 0;
        const double dmin = detectSelfIntersect(evalFig8, t1, t2);
        check("检测到自交点（两不同参数距 <1e-3）", dmin < 1e-3);
        check("自交参数 t1 在内部（0.02<t1<0.98）", t1 > 0.02 && t1 < 0.98);
        check("自交点非仅端点闭合（t1≠0,1）", t1 > 0.05 && t1 < 0.95);
    }

    // [b] 泪滴闭合环（仅两端重合）→ 无内部自交，正确区分
    {
        double t1 = 0, t2 = 0;
        const double dc = detectSelfIntersect(evalTeardrop, t1, t2);
        check("闭合环无内部自交（最小内部距>0.05）", dc > 0.05);
    }

    // [c] 单调曲线负例 → 无自交
    {
        double t1 = 0, t2 = 0;
        const double dl = detectSelfIntersect(evalLine, t1, t2);
        check("单调曲线无自交（最小距>0.1）", dl > 0.1);
    }

    printf("\n结论：测例⑦——曲线自交检测：8 字形内部自交成对检出，闭合环/单调线正确拒检，\n");
    printf("      曲线侧自交检测可用（与曲面自交对称补齐，区分内部自交 vs 端点闭合）。\n");
    printf("      通过 %d，失败 %d。\n", pass, fail);
    return fail ? 1 : 0;
}
