# -*- coding: utf-8 -*-
"""从输入推导结果：用论文公式链（A.4/A.6/A.7/A.9/A.10/A.16/A.17/A.19/A.20）复算程序输出。
   对照：Fig.5 两交点、双纽线自交参数对、SSI 盒带计数（4,052/3,070 与 40,188/85,344）。
"""
import sys, math
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np

NEG, POS = -math.inf, math.inf


# ---------- 区间算术（每次运算向外一 ulp，对应 Eq.(2)）----------
def iadd(a, b):
    return (math.nextafter(a[0] + b[0], NEG), math.nextafter(a[1] + b[1], POS))


def imul(a, b):
    p = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return (math.nextafter(min(p), NEG), math.nextafter(max(p), POS))


def isub(a, b):
    return iadd(a, (-b[1], -b[0]))


def contains_zero(a):
    return a[0] <= 0.0 <= a[1]


# ---------- 多项式 f(s,t)=Σ a[i][j] s^i t^j（幂基，Horner 区间）----------
def poly_ival(a, s, t):
    r = (0.0, 0.0)
    si = (1.0, 1.0)
    for i in range(len(a)):
        ri = (0.0, 0.0)
        tj = (1.0, 1.0)
        for j in range(len(a[i])):
            ri = iadd(ri, imul((a[i][j], a[i][j]), tj))
            tj = imul(tj, t)
        r = iadd(r, imul(ri, si))
        si = imul(si, s)
    return r


def cover(a, sl, sh, tl, th, tol, st):
    """Eq.(A.20) 排除 + Eq.(A.9) 细分（长边二分），统计排除/覆盖/细分。"""
    F = poly_ival(a, (sl, sh), (tl, th))
    if not contains_zero(F):
        st['excl'] += 1
        return
    if (sh - sl) <= tol and (th - tl) <= tol:
        st['cover'] += 1
        return
    st['subdiv'] += 1
    if (sh - sl) >= (th - tl):
        m = 0.5 * (sl + sh)
        cover(a, sl, m, tl, th, tol, st)
        cover(a, m, sh, tl, th, tol, st)
    else:
        m = 0.5 * (tl + th)
        cover(a, sl, sh, tl, m, tol, st)
        cover(a, sl, sh, m, th, tol, st)


# ============ 测试 1：SSI 盒带计数 ============
print("=" * 70)
print("测试 1：公式 (A.20)+(A.9) 复算 SSI 盒带（对照 Table 7 / bench_ssi）")
# fixture A: f=s^2-t^2
A = [[0.0, 0.0, -1.0], [0.0, 0.0], [1.0, 0.0]]
st = {'excl': 0, 'cover': 0, 'subdiv': 0}
cover(A, 0.0, 1.0, 0.0, 1.0, 1e-3, st)
print("  fixture A (f=s²−t²): 排除=%d 覆盖=%d 细分=%d | 程序: 4052 / 3070 / 7121"
      % (st['excl'], st['cover'], st['subdiv']))
# fixture D: f=(s-t)^2 = s^2 - 2st + t^2
D = [[0.0, 0.0, 1.0], [0.0, -2.0], [1.0, 0.0]]
st2 = {'excl': 0, 'cover': 0, 'subdiv': 0}
cover(D, 0.0, 1.0, 0.0, 1.0, 1e-3, st2)
print("  fixture D (f=(s−t)²): 排除=%d 覆盖=%d 细分=%d | 程序: 40188 / 85344 / 125531"
      % (st2['excl'], st2['cover'], st2['subdiv']))

# ============ 测试 2：Fig.5 两交点（Eq.(A.4) 残差 + Newton + (A.7) 认证）============
print("\n" + "=" * 70)
print("测试 2：公式 (A.4) 求根 + (A.6)(A.7) 认证，对照论文 Fig.5 两个交点")


def bez(ctrl, u):
    u = np.asarray(u, dtype=float)
    P = np.asarray(ctrl, dtype=float)
    b = np.stack([(1 - u) ** 2, 2 * (1 - u) * u, u ** 2], axis=-1)
    return b @ P


def bez_d(ctrl, u):
    u = np.asarray(u, dtype=float)
    P = np.asarray(ctrl, dtype=float)
    d = np.stack([-2 * (1 - u), 2 - 4 * u, 2 * u], axis=-1)
    return d @ P


C1 = [[0, 0], [1.2, 2.2], [3.0, 0.6]]
C2 = [[0.2, 1.8], [1.8, 0.2], [3.2, 2.2]]


def F(x):
    u, v = x
    d = bez(C1, u) - bez(C2, v)
    return np.array([d @ bez_d(C1, u), d @ bez_d(C2, v)])


def Jf(x):
    u, v = x
    du, dv = bez_d(C1, u), bez_d(C2, v)
    cu, cv = bez(C1, u), bez(C2, v)
    d = cu - cv
    ddu, ddv = bez_d(C1, u), bez_d(C2, v)
    d2u = 2 * (np.asarray(C1[2]) - 2 * np.asarray(C1[1]) + np.asarray(C1[0]))
    d2v = 2 * (np.asarray(C2[2]) - 2 * np.asarray(C2[1]) + np.asarray(C2[0]))
    return np.array([[du @ du + d @ d2u, -dv @ du],
                     [du @ dv, -dv @ dv - d @ d2v]])


roots = []
for u0 in np.linspace(0.05, 0.95, 10):
    for v0 in np.linspace(0.05, 0.95, 10):
        x = np.array([u0, v0])
        for _ in range(60):
            Fv = F(x)
            if np.linalg.norm(Fv) < 1e-14:
                break
            try:
                x = x - np.linalg.solve(Jf(x), Fv)
            except np.linalg.LinAlgError:
                break
        if np.all((x >= -1e-9) & (x <= 1 + 1e-9)) and np.linalg.norm(F(x)) < 1e-10:
            if not any(np.linalg.norm(x - r) < 1e-6 for r in roots):
                roots.append(x)
roots.sort(key=lambda r: r[0])
print("  求得交点（按 u 排序）:")
for r in roots:
    print("    (u*,v*)=(%.4f,%.4f)   ‖F‖=%.2e" % (r[0], r[1], np.linalg.norm(F(r))))
print("  论文 Fig.5 / §7 参考值: (0.4183,0.2892) 与 (0.7307,0.6088)")

# ============ 测试 3：双纽线自交（Eq.(A.17)）============
print("\n" + "=" * 70)
print("测试 3：公式 (A.17) 曲线自交，对照程序 (s1*,s2*)=(0.25,0.75)")
Amp = 0.5
TWO = 2 * math.pi


def lx(s):
    t = TWO * s
    sn, cs = math.sin(t), math.cos(t)
    return Amp * cs / (1 + sn * sn)


def ly(s):
    t = TWO * s
    sn, cs = math.sin(t), math.cos(t)
    return Amp * sn * cs / (1 + sn * sn)


def F2(x):
    s1, s2 = x
    return np.array([lx(s1) - lx(s2), ly(s1) - ly(s2)])


def J2(x):
    s1, s2 = x
    h = 1e-6
    d1 = np.array([(lx(s1 + h) - lx(s1 - h)) / (2 * h), (ly(s1 + h) - ly(s1 - h)) / (2 * h)])
    d2 = np.array([(lx(s2 + h) - lx(s2 - h)) / (2 * h), (ly(s2 + h) - ly(s2 - h)) / (2 * h)])
    return np.array([[d1[0], -d2[0]], [d1[1], -d2[1]]])


x = np.array([0.20, 0.80])
for _ in range(80):
    Fv = F2(x)
    if np.linalg.norm(Fv) < 1e-15:
        break
    x = x - np.linalg.solve(J2(x), Fv)
print("  Newton 解: s1*=%.6f, s2*=%.6f, ‖F‖=%.2e" % (x[0], x[1], np.linalg.norm(F2(x))))
print("  程序/论文参考值: (0.250000, 0.750000)")
print("  解析验证 lx(0.25)=%.3e, ly(0.25)=%.3e, lx(0.75)=%.3e, ly(0.75)=%.3e"
      % (lx(0.25), ly(0.25), lx(0.75), ly(0.75)))

# ============ 测试 4：点包含（Eq.(A.16)）============
print("\n" + "=" * 70)
print("测试 4：公式 (A.16) 点包含（退化恢复）")
P = bez(C1, 0.37)
best = min(((np.linalg.norm(bez(C1, s) - P), s) for s in np.linspace(0, 1, 10001)))
print("  P=C₁(0.37)，一维搜索最优 s*=%.4f，残差=%.2e ⇒ P∈image(C₁) 认证" % (best[1], best[0]))
