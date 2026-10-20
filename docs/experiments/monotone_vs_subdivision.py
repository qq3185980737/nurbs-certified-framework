# -*- coding: utf-8 -*-
"""1D 求根策略对比（分段三次 Bézier，20 个跨度各含 1 个根）。

A) 区间细分 + 排除（enumerateIntersections 的 1D 类比，逐跨度凸包区间求值）
B) 求导定单调区间 + 二分
B') 单调段 + 牛顿（二分兜底）
C) 区间细分 + 单调性剪枝 —— 注意：剪枝判据必须按跨度边界处理，否则会把跨段盒误剪；
   本原型未做该处理，实测漏 1 个根（19/20），**仅作对照，不用于任何结论**。

结论（A vs B）：根越密，极值点越多，先切出"导数定号"的单调段本身就比直接细分更贵，
因此本框架采用区间细分 + 凸包排除；单调性只在导数区间定号时作为廉价剪枝候选。
统计口径：区间求值次数（≈ 盒评估成本）。
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')

NP = 20                      # 20 个跨度，每段 1 个根
POLES_Y = []                 # 每个跨度的 4 个 y 控制值（三次 Bézier）
SPANS = []
for k in range(NP):
    x0, x1 = k / NP, (k + 1) / NP
    # y 值在段内单调穿越 0：端点异号
    y0 = 1.0 if k % 2 == 0 else -1.0
    y1 = -y0
    # 控制点取单调插值（保证段内单调、端点符号相反）
    p = [y0, y0 * 0.35 + y1 * 0.0, y1 * 0.65 + y0 * 0.0, y1]
    p = [y0, y0 * 0.6, y1 * 0.6, y1]
    SPANS.append((x0, x1, p))

def hull_span(p, a, b, x0, x1):
    """在 [a,b] ⊂ [x0,x1] 上做 de Casteljau（Bernstein 凸包）"""
    t0 = (a - x0) / (x1 - x0)
    t1 = (b - x0) / (x1 - x0)
    def split(c, t):
        cur = list(c); left = [cur[0]]; right = [cur[-1]]
        while len(cur) > 1:
            cur = [(1 - t) * cur[i] + t * cur[i + 1] for i in range(len(cur) - 1)]
            left.append(cur[0]); right.append(cur[-1])
        return left, right[::-1]
    c = p
    if t0 > 0:
        _, c = split(c, t0)
        t1 = (t1 - t0) / (1 - t0)
    if t1 < 1:
        c, _ = split(c, t1)
    return min(c), max(c)

def deriv_bernstein(p):
    n = len(p) - 1
    return [n * (p[i + 1] - p[i]) for i in range(n)]

def f_at(x):
    k = min(int(x * NP), NP - 1)
    x0, x1, p = SPANS[k]
    t = (x - x0) / (x1 - x0)
    return ((1 - t) ** 3) * p[0] + 3 * ((1 - t) ** 2) * t * p[1] + 3 * (1 - t) * t * t * p[2] + t ** 3 * p[3]

def eval_box(lo, hi):
    """跨度为单个 Bézier 段时直接求包络；否则按跨度分解"""
    k0 = min(int(lo * NP), NP - 1)
    k1 = min(int(hi * NP), NP - 1)
    if k0 != k1:
        return (-1e9, 1e9)                      # 跨段：保守不排除
    x0, x1, p = SPANS[k0]
    return hull_span(p, lo, hi, x0, x1)

def eval_dbox(lo, hi):
    k0 = min(int(lo * NP), NP - 1)
    k1 = min(int(hi * NP), NP - 1)
    if k0 != k1:
        return (-1e9, 1e9)
    x0, x1, p = SPANS[k0]
    return hull_span(deriv_bernstein(p), lo, hi, x0, x1)

# ---------- A) 细分 ----------
def method_A(tol=1e-4):
    evals = 0; kept = 0; boxes = []
    stack = [(0.0, 1.0)]
    while stack:
        lo, hi = stack.pop()
        evals += 1
        a, b = eval_box(lo, hi)
        if not (a <= 0.0 <= b):
            continue
        if hi - lo <= tol:
            kept += 1; boxes.append((lo, hi)); continue
        m = 0.5 * (lo + hi)
        stack.append((lo, m)); stack.append((m, hi))
    return evals, kept, boxes

# ---------- B) 单调 + 二分 ----------
def method_B(tol=1e-14):
    evals = 0; slabs = []; stack = [(0.0, 1.0)]
    while stack:
        lo, hi = stack.pop()
        evals += 1
        a, b = eval_dbox(lo, hi)
        if a > 0 or b < 0 or (hi - lo) < 1e-9:
            slabs.append((lo, hi)); continue
        m = 0.5 * (lo + hi)
        stack.append((lo, m)); stack.append((m, hi))
    roots = []
    for lo, hi in slabs:
        evals += 2
        fa, fb = f_at(lo), f_at(hi)
        if (fa > 0) == (fb > 0):
            continue
        a, b, fa = lo, hi, fa
        while b - a > tol:
            evals += 1
            m = 0.5 * (a + b); fm = f_at(m)
            if (fm > 0) == (fa > 0):
                a, fa = m, fm
            else:
                b = m
        roots.append(0.5 * (a + b))
    return evals, len(slabs), roots

# ---------- C) 细分 + 单调剪枝 ----------
def method_C(tol=1e-4):
    evals = 0; kept = 0; pruned = 0; boxes = []
    stack = [(0.0, 1.0)]
    while stack:
        lo, hi = stack.pop()
        evals += 1
        a, b = eval_box(lo, hi)
        if not (a <= 0.0 <= b):
            continue
        da, db = eval_dbox(lo, hi)
        if da > 0 or db < 0:
            evals += 2
            if (f_at(lo) > 0) == (f_at(hi) > 0):
                pruned += 1; continue
        if hi - lo <= tol:
            kept += 1; boxes.append((lo, hi)); continue
        m = 0.5 * (lo + hi)
        stack.append((lo, m)); stack.append((m, hi))
    return evals, kept, boxes, pruned

TRUE_ROOTS = [((k + 0.5) / NP) for k in range(NP)]

def covered(boxes):
    return sum(1 for r in TRUE_ROOTS if any(lo - 1e-9 <= r <= hi + 1e-9 for lo, hi in boxes))

A = method_A(); B = method_B(); C = method_C()
print("fixture: 20 段三次 Bézier，每段 1 个零点（共 20 个）")
print("A 细分          : 求值 %5d, 保留盒 %4d, 覆盖 %d/20" % (A[0], A[1], covered(A[2])))
print("B 单调+二分     : 求值 %5d, 单调段 %4d, 找到 %d/20" % (B[0], B[1], len(B[2])))
print("C 细分+单调剪枝 : 求值 %5d, 保留盒 %4d, 剪枝 %3d, 覆盖 %d/20" % (C[0], C[1], C[3], covered(C[2])))

# ---------- B') 单调段 + 牛顿(带二分兜底) ----------
def fprime_at(x):
    k = min(int(x * NP), NP - 1)
    x0, x1, p = SPANS[k]
    t = (x - x0) / (x1 - x0)
    d = deriv_bernstein(p)
    return ((1 - t) ** 2) * d[0] + 2 * (1 - t) * t * d[1] + t * t * d[2]

def method_Bp(tol=1e-14):
    evals = 0; slabs = []; stack = [(0.0, 1.0)]
    while stack:
        lo, hi = stack.pop()
        evals += 1
        a, b = eval_dbox(lo, hi)
        if a > 0 or b < 0 or (hi - lo) < 1e-9:
            slabs.append((lo, hi)); continue
        m = 0.5 * (lo + hi)
        stack.append((lo, m)); stack.append((m, hi))
    roots = []
    for lo, hi in slabs:
        evals += 2
        fa, fb = f_at(lo), f_at(hi)
        if (fa > 0) == (fb > 0):
            continue
        a, b, fa = lo, hi, fa
        for _ in range(60):
            evals += 1
            x = b - fb * (b - a) / (fb - fa) if fb != fa else 0.5 * (a + b)   # 割线/牛顿式
            fp = fprime_at(x); evals += 1
            xn = x - f_at(x) / fp if fp != 0 else 0.5 * (a + b)
            if not (a < xn < b):
                xn = 0.5 * (a + b)
            fn = f_at(xn); evals += 1
            if fn == 0.0 or abs(b - a) < tol:
                a = b = xn; break
            if (fn > 0) == (fa > 0):
                a, fa = xn, fn
            else:
                b, fb = xn, fn
        roots.append(0.5 * (a + b))
    return evals, len(slabs), roots

Bp = method_Bp()
print("B' 单调+牛顿二分 : 求值 %5d, 单调段 %4d, 找到 %d/20" % (Bp[0], Bp[1], len(Bp[2])))
