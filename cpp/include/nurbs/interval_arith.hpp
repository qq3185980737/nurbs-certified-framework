// ============================================================================
// interval_arith.hpp — 严格舍入的双精度区间算术（IEEE-754，向下/向上舍入）
// ----------------------------------------------------------------------------
// 供"保证性求交"使用：对 BVH 定位到的收敛解构造参数盒，用区间扩展的
// 函数值/一阶导做 Krawczyk 存在性验证。
// 舍入策略（跨平台一致保证"结果严格包含真值"）：
//   * MSVC 原生（论文实验平台）：_controlfp_s 切换 _RC_DOWN/_RC_UP，严格向外舍入；
//   * Emscripten/WASM 及其他非 MSVC 平台：不切换硬件舍入模式（WASM 不支持），
//     改用 std::nextafter 对区间端点各外扩 1 ulp，同样保证严格包含真值。
// ============================================================================
#ifndef NURBS_INTERVAL_ARITH_HPP
#define NURBS_INTERVAL_ARITH_HPP

#include <float.h>
#include <cmath>
#include <algorithm>
#include <limits>

namespace nurbs {

#ifdef _MSC_VER
#  define NURBS_NATIVE_ROUNDING 1
#endif

struct Interval {
    double lo, hi;
    Interval(double l = 0.0, double h = 0.0) : lo(l), hi(h) {}
    static Interval point(double x) { return Interval(x, x); }
    static Interval empty() { return Interval(1e300, -1e300); }
    bool isEmpty() const { return lo > hi; }
    bool contains(double x) const { return lo <= x && x <= hi; }
    double mid() const { return 0.5 * (lo + hi); }
    double width() const { return hi - lo; }
};

#ifdef NURBS_NATIVE_ROUNDING
#  define NURBS_RC_DOWN _RC_DOWN
#  define NURBS_RC_UP   _RC_UP
#  define NURBS_OUTWARD(lo, hi) ((void)0)
struct RoundingGuard {
    unsigned old_;
    explicit RoundingGuard(unsigned rc) { _controlfp_s(&old_, rc, _MCW_RC); }
    ~RoundingGuard() { unsigned t; _controlfp_s(&t, old_, _MCW_RC); }
};
#else
#  define NURBS_RC_DOWN 0
#  define NURBS_RC_UP   1
#  define NURBS_OUTWARD(lo, hi) \
     do { (lo) = std::nextafter((lo), -std::numeric_limits<double>::infinity()); \
          (hi) = std::nextafter((hi),  std::numeric_limits<double>::infinity()); } while (0)
struct RoundingGuard { explicit RoundingGuard(unsigned) {} ~RoundingGuard() {} };
#endif

inline Interval iadd(const Interval& a, const Interval& b) {
    Interval r;
    { RoundingGuard g(NURBS_RC_DOWN); r.lo = a.lo + b.lo; }
    { RoundingGuard g(NURBS_RC_UP);   r.hi = a.hi + b.hi; }
    NURBS_OUTWARD(r.lo, r.hi);
    return r;
}
inline Interval isub(const Interval& a, const Interval& b) {
    Interval r;
    { RoundingGuard g(NURBS_RC_DOWN); r.lo = a.lo - b.hi; }
    { RoundingGuard g(NURBS_RC_UP);   r.hi = a.hi - b.lo; }
    NURBS_OUTWARD(r.lo, r.hi);
    return r;
}
inline Interval imul(const Interval& a, const Interval& b) {
    Interval r;
    { RoundingGuard g(NURBS_RC_DOWN);
      double m1 = a.lo * b.lo, m2 = a.lo * b.hi, m3 = a.hi * b.lo, m4 = a.hi * b.hi;
      r.lo = std::fmin(std::fmin(m1, m2), std::fmin(m3, m4)); }
    { RoundingGuard g(NURBS_RC_UP);
      double m1 = a.lo * b.lo, m2 = a.lo * b.hi, m3 = a.hi * b.lo, m4 = a.hi * b.hi;
      r.hi = std::fmax(std::fmax(m1, m2), std::fmax(m3, m4)); }
    NURBS_OUTWARD(r.lo, r.hi);
    return r;
}
inline Interval idiv(const Interval& a, const Interval& b) {
    if (b.lo <= 0.0 && b.hi >= 0.0) return Interval::empty();
    Interval r;
    { RoundingGuard g(NURBS_RC_DOWN);
      r.lo = std::fmin(std::fmin(a.lo / b.lo, a.lo / b.hi),
                       std::fmin(a.hi / b.lo, a.hi / b.hi)); }
    { RoundingGuard g(NURBS_RC_UP);
      r.hi = std::fmax(std::fmax(a.lo / b.lo, a.lo / b.hi),
                       std::fmax(a.hi / b.lo, a.hi / b.hi)); }
    NURBS_OUTWARD(r.lo, r.hi);
    return r;
}
inline Interval isqrt(const Interval& a) {
    if (a.lo < 0.0) return Interval::empty();
    Interval r;
    { RoundingGuard g(NURBS_RC_DOWN); r.lo = std::sqrt(a.lo); }
    { RoundingGuard g(NURBS_RC_UP);   r.hi = std::sqrt(a.hi); }
    NURBS_OUTWARD(r.lo, r.hi);
    return r;
}
inline Interval isqr(const Interval& a) {
    Interval r;
    { RoundingGuard g(NURBS_RC_DOWN);
      r.lo = std::fmin(a.lo * a.lo, a.hi * a.hi);
      if (a.lo < 0.0 && a.hi > 0.0) r.lo = 0.0; }
    { RoundingGuard g(NURBS_RC_UP);
      r.hi = std::fmax(a.lo * a.lo, a.hi * a.hi); }
    NURBS_OUTWARD(r.lo, r.hi);
    return r;
}

// 2×2 区间矩阵（Krawczyk 算子用）
struct IntervalMat22 {
    Interval m00, m01, m10, m11;
};

// 2 维区间向量
struct IntervalVec2 { Interval x, y; };

inline IntervalVec2 iaddVec(const IntervalVec2& a, const IntervalVec2& b) {
    return { iadd(a.x, b.x), iadd(a.y, b.y) };
}
inline IntervalVec2 matVec(const IntervalMat22& M, const IntervalVec2& v) {
    // (M00*v.x + M01*v.y, M10*v.x + M11*v.y)
    return { iadd(imul(M.m00, v.x), imul(M.m01, v.y)),
             iadd(imul(M.m10, v.x), imul(M.m11, v.y)) };
}

} // namespace nurbs
#endif
