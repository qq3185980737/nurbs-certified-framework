#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ssi_verifier.py — 路线 C 原型 v2：独立精确有理数验证器（支持有理权重 w≠1）

证书 JSON：
  patches: 两个有理 Bezier 面片（齐次控制网 [wx,wy,wz,w]，hex float → 精确 Fraction）
  tree:    4D 参数盒细分树（节点 = 按某坐标中点二分；叶子 = excluded 或 output）
  certified: 切片存在性盒（Krawczyk）

验证（全部精确有理数算术，无浮点）：
  1) 树结构：根 [0,1]^4；子盒 = 父盒按中点精确二分；坐标合法
  2) excluded 叶子：验证器自行重建子片，用有理控制点凸包 AABB 精确判定严格不相交
  3) output 叶子：盒在 [0,1]^4 内
  4) certified 盒：验证器自行重建面片，精确有理数切片 Krawczyk 通过（存在性）
用法：python3 ssi_verifier.py cert.json
"""
import hashlib, json
import sys
from fractions import Fraction as F


def parse_rat(s):
    """证书中的有理数：\"n/d\" 或 hex float 字符串"""
    if "/" in s:
        return F(s)
    return F.from_float(float.fromhex(s))


def parse_hex(s):
    # 向后兼容旧证书：hex 与 n/d 均可
    return parse_rat(s)


def _seg_dist2(p, a, b):
    """精确有理数平方距离：3D 点 p 到线段 ab（全部 Fraction）"""
    dx, dy, dz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    L2 = dx * dx + dy * dy + dz * dz
    if L2 == 0:
        return (p[0] - a[0]) ** 2 + (p[1] - a[1]) ** 2 + (p[2] - a[2]) ** 2
    dot = (p[0] - a[0]) * dx + (p[1] - a[1]) * dy + (p[2] - a[2]) * dz
    if dot <= 0:
        return (p[0] - a[0]) ** 2 + (p[1] - a[1]) ** 2 + (p[2] - a[2]) ** 2
    if dot >= L2:
        return (p[0] - b[0]) ** 2 + (p[1] - b[1]) ** 2 + (p[2] - b[2]) ** 2
    base = (p[0] - a[0]) ** 2 + (p[1] - a[1]) ** 2 + (p[2] - a[2]) ** 2
    return base - dot * dot / L2


class FI:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi

    def __add__(self, o):
        return FI(self.lo + o.lo, self.hi + o.hi)

    def __sub__(self, o):
        return FI(self.lo - o.hi, self.hi - o.lo)

    def __mul__(self, o):
        p = [self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi]
        return FI(min(p), max(p))

    def __neg__(self):
        return FI(-self.hi, -self.lo)

    def inv(self):
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError("interval contains 0")
        return FI(1 / self.hi, 1 / self.lo)

    def __truediv__(self, o):
        return self * o.inv()

    def width(self):
        return self.hi - self.lo


ZERO = FI(F(0), F(0))
ONE = FI(F(1), F(1))


def load_patch(obj):
    nu, nv = obj["nu"], obj["nv"]
    poles = [[[parse_rat(v) for v in p] for p in obj["poles"][i * nv:(i + 1) * nv]] for i in range(nu)]
    for i in range(nu):
        for j in range(nv):
            if poles[i][j][3] <= 0:
                raise ValueError("non-positive weight")
    return {"nu": nu, "nv": nv, "P": poles}


def split_1d(coeffs):
    """coeffs: list of 4-vector Fractions；返回 (左半, 右半) 精确 de Casteljau"""
    n = len(coeffs)
    tri = [coeffs]
    for r in range(1, n):
        row = []
        for i in range(n - r):
            row.append([(tri[r - 1][i][c] + tri[r - 1][i + 1][c]) / 2 for c in range(4)])
        tri.append(row)
    L = [tri[r][0] for r in range(n)]
    R = [tri[n - 1 - r][r] for r in range(n)]
    return L, R


def split_u(P):
    nu, nv = P["nu"], P["nv"]
    Lp = [[None] * nv for _ in range(nu)]
    Rp = [[None] * nv for _ in range(nu)]
    for j in range(nv):
        L, R = split_1d([P["P"][i][j] for i in range(nu)])
        for i in range(nu):
            Lp[i][j] = L[i]
            Rp[i][j] = R[i]
    return {"nu": nu, "nv": nv, "P": Lp}, {"nu": nu, "nv": nv, "P": Rp}


def split_v(P):
    nu, nv = P["nu"], P["nv"]
    Lp = [[None] * nv for _ in range(nu)]
    Rp = [[None] * nv for _ in range(nu)]
    for i in range(nu):
        L, R = split_1d([P["P"][i][j] for j in range(nv)])
        for j in range(nv):
            Lp[i][j] = L[j]
            Rp[i][j] = R[j]
    return {"nu": nu, "nv": nv, "P": Lp}, {"nu": nu, "nv": nv, "P": Rp}


def aabb(P):
    """有理控制点凸包 AABB（正权 ⇒ 曲面含于凸包）"""
    lo = [None, None, None]
    hi = [None, None, None]
    for i in range(P["nu"]):
        for j in range(P["nv"]):
            p = P["P"][i][j]
            w = p[3]
            for c in range(3):
                v = p[c] / w
                lo[c] = v if lo[c] is None else min(lo[c], v)
                hi[c] = v if hi[c] is None else max(hi[c], v)
    return lo, hi


def disjoint(la, ha, lb, hb):
    return any(ha[c] < lb[c] or hb[c] < la[c] for c in range(3))


def fi_net(P):
    return [[[FI(P["P"][i][j][c], P["P"][i][j][c]) for c in range(4)]
             for j in range(P["nv"])] for i in range(P["nu"])]


def eval_net4(net, U, V):
    nu, nv = len(net), len(net[0])
    cols = []
    for j in range(nv):
        row = [net[i][j] for i in range(nu)]
        for r in range(1, nu):
            row = [[row[i][c] + (row[i + 1][c] - row[i][c]) * U for c in range(4)]
                   for i in range(len(row) - 1)]
        cols.append(row[0])
    row = cols
    for r in range(1, nv):
        row = [[row[i][c] + (row[i + 1][c] - row[i][c]) * V for c in range(4)]
               for i in range(len(row) - 1)]
    return row[0]


def deriv_net4(net, deg, axis):
    nu, nv = len(net), len(net[0])
    dd = FI(F(deg), F(deg))
    if axis == 0:
        return [[[(net[i + 1][j][c] - net[i][j][c]) * dd for c in range(4)]
                 for j in range(nv)] for i in range(nu - 1)]
    return [[[(net[i][j + 1][c] - net[i][j][c]) * dd for c in range(4)]
             for j in range(nv - 1)] for i in range(nu)]


def deriv_net4_frac(net, deg, axis):
    """Fraction 控制网的导数控制网（用于精确点导数）"""
    nu, nv = len(net), len(net[0])
    if axis == 0:
        return [[[(net[i + 1][j][c] - net[i][j][c]) * deg for c in range(4)]
                 for j in range(nv)] for i in range(nu - 1)]
    return [[[(net[i][j + 1][c] - net[i][j][c]) * deg for c in range(4)]
             for j in range(nv - 1)] for i in range(nu)]


def point_eval4(net, u, v):
    nu, nv = len(net), len(net[0])
    cols = []
    for j in range(nv):
        row = [net[i][j] for i in range(nu)]
        for r in range(1, nu):
            row = [[row[i][c] + (row[i + 1][c] - row[i][c]) * u for c in range(4)]
                   for i in range(len(row) - 1)]
        cols.append(row[0])
    row = cols
    for r in range(1, nv):
        row = [[row[i][c] + (row[i + 1][c] - row[i][c]) * v for c in range(4)]
               for i in range(len(row) - 1)]
    return row[0]


class Patch:
    def __init__(self, obj):
        self.nu, self.nv = obj["nu"], obj["nv"]
        self.P = [[[parse_rat(x) for x in obj["poles"][i * self.nv + j]] for j in range(self.nv)]
                  for i in range(self.nu)]
        self.positive = all(p[3] > 0 for row in self.P for p in row)
        self.wmin = min(p[3] for row in self.P for p in row)
        self.wmax = max(p[3] for row in self.P for p in row)
        self.H = fi_net({"nu": self.nu, "nv": self.nv, "P": self.P})
        self.Hu = deriv_net4(self.H, self.nu - 1, 0) if self.nu > 1 else None
        self.Hv = deriv_net4(self.H, self.nv - 1, 1) if self.nv > 1 else None

    def point(self, u, v):
        h = point_eval4(self.P, u, v)
        return [h[c] / h[3] for c in range(3)]

    def _clampw(self, h3):
        lo = h3.lo if h3.lo > self.wmin else self.wmin
        hi = h3.hi if h3.hi < self.wmax else self.wmax
        return FI(lo, hi)

    def cart_iv(self, U, V):
        h = eval_net4(self.H, U, V)
        w = self._clampw(h[3])
        return [h[c] / w for c in range(3)]

    def du_iv(self, U, V):
        if self.Hu is None:
            return [ZERO, ZERO, ZERO]
        h = eval_net4(self.H, U, V)
        hu = eval_net4(self.Hu, U, V)
        w = self._clampw(h[3])
        out = []
        for c in range(3):
            num = hu[c] * w - h[c] * hu[3]
            out.append(num / (w * w))
        return out

    def dv_iv(self, U, V):
        if self.Hv is None:
            return [ZERO, ZERO, ZERO]
        h = eval_net4(self.H, U, V)
        hv = eval_net4(self.Hv, U, V)
        w = self._clampw(h[3])
        out = []
        for c in range(3):
            num = hv[c] * w - h[c] * hv[3]
            out.append(num / (w * w))
        return out

    def _deriv_frac(self, u, v, axis):
        if axis == 0:
            net = deriv_net4_frac(self.P, self.nu - 1, 0) if self.nu > 1 else None
        else:
            net = deriv_net4_frac(self.P, self.nv - 1, 1) if self.nv > 1 else None
        if net is None:
            return [F(0), F(0), F(0)]
        h = point_eval4(self.P, u, v)
        dh = point_eval4(net, u, v)
        return [(dh[c] * h[3] - h[c] * dh[3]) / (h[3] * h[3]) for c in range(3)]

    def du_point(self, u, v):
        return self._deriv_frac(u, v, 0)

    def dv_point(self, u, v):
        return self._deriv_frac(u, v, 1)


def inv3(M):
    det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
           - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
           + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    if det == 0:
        return None
    cof = [[M[1][1] * M[2][2] - M[1][2] * M[2][1], -(M[1][0] * M[2][2] - M[1][2] * M[2][0]), M[1][0] * M[2][1] - M[1][1] * M[2][0]],
           [-(M[0][1] * M[2][2] - M[0][2] * M[2][1]), M[0][0] * M[2][2] - M[0][2] * M[2][0], -(M[0][0] * M[2][1] - M[0][1] * M[2][0])],
           [M[0][1] * M[1][2] - M[0][2] * M[1][1], -(M[0][0] * M[1][2] - M[0][2] * M[1][0]), M[0][0] * M[1][1] - M[0][1] * M[1][0]]]
    inv = [[F(0)] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            inv[j][i] = cof[i][j] / det
    return inv


def _krawczyk_ok_single(PA, PB, U, V, S, T):
    SA = PA.cart_iv(U, V)
    SB = PB.cart_iv(S, T)
    Au = PA.du_iv(U, V)
    Av = PA.dv_iv(U, V)
    Bu = PB.du_iv(S, T)
    Bv = PB.dv_iv(S, T)
    J = [[Au[c], Av[c], Bu[c].__neg__(), Bv[c].__neg__()] for c in range(3)]
    Vx = [U, V, S, T]
    m = [(x.lo + x.hi) / 2 for x in Vx]
    pmA = PA.point(m[0], m[1])
    pmB = PB.point(m[2], m[3])
    Fm = [pmA[c] - pmB[c] for c in range(3)]
    for kfix in range(4):
        idx = [c for c in range(4) if c != kfix]
        J3 = [[J[r][idx[c]] for c in range(3)] for r in range(3)]
        J3m = [[(J3[r][c].lo + J3[r][c].hi) / 2 for c in range(3)] for r in range(3)]
        Y = inv3(J3m)
        if Y is None:
            continue
        X3 = [Vx[i] for i in idx]
        m3 = [m[i] for i in idx]
        W = [[ZERO] * 3 for _ in range(3)]
        try:
            for r in range(3):
                for d in range(3):
                    acc = ONE if r == d else ZERO
                    for c in range(3):
                        acc = acc - FI(Y[r][c], Y[r][c]) * J3[c][d]
                    W[r][d] = acc
            K = [ZERO] * 3
            for r in range(3):
                acc = FI(m3[r], m3[r])
                for c in range(3):
                    acc = acc - FI(Y[r][c], Y[r][c]) * FI(Fm[c], Fm[c])
                for d in range(3):
                    acc = acc + W[r][d] * FI(X3[d].lo - m3[d], X3[d].hi - m3[d])
                K[r] = acc
        except ZeroDivisionError:
            continue
        inside = True
        for r in range(3):
            if not (X3[r].lo <= K[r].lo and K[r].hi <= X3[r].hi):
                inside = False
                break
            wd = K[r].width()
            xw = X3[r].width()
            if not (wd <= F(9, 10) * xw or xw < F(1, 10 ** 12)):
                inside = False
                break
        if inside:
            return True
    return False

def _krawczyk_ok_iter(PA, PB, U, V, S, T):
    """与内核 krawczyk_slice 语义一致：逐 kfix 迭代收缩（最多 12 次），精确有理数验证。"""
    X0 = [U, V, S, T]
    for kfix in range(4):
        X = list(X0)
        idx = [c for c in range(4) if c != kfix]
        for _ in range(12):
            try:
                Au = PA.du_iv(X[0], X[1]); Av = PA.dv_iv(X[0], X[1])
                Bu = PB.du_iv(X[2], X[3]); Bv = PB.dv_iv(X[2], X[3])
            except ZeroDivisionError:
                return False
            J = [[Au[c], Av[c], Bu[c].__neg__(), Bv[c].__neg__()] for c in range(3)]
            m = [(x.lo + x.hi) / 2 for x in X]
            pmA = PA.point(m[0], m[1]); pmB = PB.point(m[2], m[3])
            Fm = [pmA[c] - pmB[c] for c in range(3)]
            J3 = [[J[r][idx[c]] for c in range(3)] for r in range(3)]
            J3m = [[(J3[r][c].lo + J3[r][c].hi) / 2 for c in range(3)] for r in range(3)]
            Y = inv3(J3m)
            if Y is None:
                break
            X3 = [X[i] for i in idx]; m3 = [m[i] for i in idx]
            W = [[ZERO] * 3 for _ in range(3)]
            for r in range(3):
                for d in range(3):
                    acc = ONE if r == d else ZERO
                    for c in range(3):
                        acc = acc - FI(Y[r][c], Y[r][c]) * J3[c][d]
                    W[r][d] = acc
            K = [ZERO] * 3
            for r in range(3):
                acc = FI(m3[r], m3[r])
                for c in range(3):
                    acc = acc - FI(Y[r][c], Y[r][c]) * FI(Fm[c], Fm[c])
                for d in range(3):
                    acc = acc + W[r][d] * FI(X3[d].lo - m3[d], X3[d].hi - m3[d])
                K[r] = acc
            inside = True
            for r in range(3):
                if not (X3[r].lo <= K[r].lo and K[r].hi <= X3[r].hi):
                    inside = False; break
                xw = X3[r].width()
                if not (K[r].width() <= F(9, 10) * xw or xw < F(1, 10 ** 12)):
                    inside = False; break
            if inside:
                return True
            nonempty = True
            narrower = False
            newX3 = []
            for r in range(3):
                lo = max(X3[r].lo, K[r].lo); hi = min(X3[r].hi, K[r].hi)
                if lo > hi:
                    nonempty = False; break
                if lo > X3[r].lo or hi < X3[r].hi:
                    narrower = True
                newX3.append(FI(lo, hi))
            if (not nonempty) or (not narrower):
                break
            for i in range(3):
                X[idx[i]] = newX3[i]
    return False


def krawczyk_ok(PA, PB, U, V, S, T):
    """精确单次切片 Krawczyk 检查（快；claimed_slice 用 Newton polish + 小盒验证）。"""
    return _krawczyk_ok_single(PA, PB, U, V, S, T)


def _insert_u(P, U, nu, nv, du, t):
    """精确有理数 Boehm 节点插入（u 方向），与内核 insert_knot_u 同算法"""
    k = du
    for i in range(du, nu):
        if U[i] <= t < U[i + 1]:
            k = i
            break
    Pout = [None] * ((nu + 1) * nv)
    Uout = list(U[:k + 1]) + [t] + list(U[k + 1:])
    for j in range(nv):
        for i in range(0, k - du + 1):
            Pout[i * nv + j] = P[i * nv + j]
        for i in range(k - du + 1, k + 1):
            den = U[i + du] - U[i]
            alpha = (t - U[i]) / den if den != 0 else F(0)
            Pout[i * nv + j] = [alpha * P[i * nv + j][c] + (1 - alpha) * P[(i - 1) * nv + j][c]
                                for c in range(4)]
        for i in range(k + 1, nu + 1):
            Pout[i * nv + j] = P[(i - 1) * nv + j]
    return Pout, Uout, nu + 1


def _insert_v(P, V, nu, nv, dv, t):
    k = dv
    for i in range(dv, nv):
        if V[i] <= t < V[i + 1]:
            k = i
            break
    Pout = [None] * (nu * (nv + 1))
    Vout = list(V[:k + 1]) + [t] + list(V[k + 1:])
    for i in range(nu):
        for j in range(0, k - dv + 1):
            Pout[i * (nv + 1) + j] = P[i * nv + j]
        for j in range(k - dv + 1, k + 1):
            den = V[j + dv] - V[j]
            alpha = (t - V[j]) / den if den != 0 else F(0)
            Pout[i * (nv + 1) + j] = [alpha * P[i * nv + j][c] + (1 - alpha) * P[i * nv + j - 1][c]
                                      for c in range(4)]
        for j in range(k + 1, nv + 1):
            Pout[i * (nv + 1) + j] = P[i * nv + j - 1]
    return Pout, Vout, nv + 1


def p7_convert(nur):
    """精确有理数 NURBS -> Bezier 分块（Boehm 提升内部节点重数到次数），返回块控制网列表"""
    du, dv, nu, nv = nur["du"], nur["dv"], nur["nu"], nur["nv"]
    P = [[parse_rat(x) for x in nur["poles"][i * nv + j]] for i in range(nu) for j in range(nv)]
    U = [parse_rat(x) for x in nur["U"]]
    V = [parse_rat(x) for x in nur["V"]]
    # 检查夹紧且数量一致
    for t in sorted({u for u in U[1:-1]}):
        mult = U.count(t)
        while mult < du:
            P, U, nu = _insert_u(P, U, nu, nv, du, t)
            mult += 1
    for t in sorted({v for v in V[1:-1]}):
        mult = V.count(t)
        while mult < dv:
            P, V, nv = _insert_v(P, V, nu, nv, dv, t)
            mult += 1
    uk = [U[i] for i in range(len(U) - 1) if U[i] < U[i + 1]] + [U[-1]]
    vk = [V[i] for i in range(len(V) - 1) if V[i] < V[i + 1]] + [V[-1]]
    blocks = []
    for su in range(len(uk) - 1):
        for sv in range(len(vk) - 1):
            poles = [P[(su * du + i) * nv + (sv * dv + j)]
                     for i in range(du + 1) for j in range(dv + 1)]
            blocks.append(poles)
    return blocks


def krawczyk4d_ok(PA, PB, U, V, S, T):
    """精确有理数 4D 欠定 Krawczyk：K ⊆ X ⇒ X 内存在 F=0 解（Y 满列秩 3）"""
    SA = PA.cart_iv(U, V)
    SB = PB.cart_iv(S, T)
    Au = PA.du_iv(U, V)
    Av = PA.dv_iv(U, V)
    Bu = PB.du_iv(S, T)
    Bv = PB.dv_iv(S, T)
    J = [[Au[c], Av[c], Bu[c].__neg__(), Bv[c].__neg__()] for c in range(3)]
    X = [U, V, S, T]
    m = [(x.lo + x.hi) / 2 for x in X]
    pmA = PA.point(m[0], m[1])
    pmB = PB.point(m[2], m[3])
    Fm = [pmA[c] - pmB[c] for c in range(3)]
    Jm = [[(J[r][d].lo + J[r][d].hi) / 2 for d in range(4)] for r in range(3)]
    M = [[sum(Jm[r][d] * Jm[s][d] for d in range(4)) for s in range(3)] for r in range(3)]
    Minv = inv3(M)
    if Minv is None:
        return False
    Y = [[sum(Jm[d][a] * Minv[d][c] for d in range(3)) for c in range(3)] for a in range(4)]
    W = [[ZERO] * 4 for _ in range(4)]
    for a in range(4):
        for b in range(4):
            acc = ONE if a == b else ZERO
            for c in range(3):
                acc = acc - FI(Y[a][c], Y[a][c]) * J[c][b]
            W[a][b] = acc
    K = [ZERO] * 4
    for a in range(4):
        acc = FI(m[a], m[a])
        for c in range(3):
            acc = acc - FI(Y[a][c], Y[a][c]) * FI(Fm[c], Fm[c])
        for b in range(4):
            acc = acc + W[a][b] * FI(X[b].lo - m[b], X[b].hi - m[b])
        K[a] = acc
    for a in range(4):
        if not (X[a].lo <= K[a].lo and K[a].hi <= X[a].hi):
            return False
    return True


def _cross3(a, b):
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]]


def check_coincidence(cert, PA, PB):
    """精确验证平面双线性重合证书：仿射参数对应 + 完整齐次控制网 + 盒映射"""
    errors = []
    co = cert["coincidence"]
    al, be, ga = parse_rat(co["alpha"]), parse_rat(co["beta"]), parse_rat(co["gamma"])
    de, ep, ze = parse_rat(co["delta"]), parse_rat(co["epsilon"]), parse_rat(co["zeta"])
    if PA.nu != 2 or PA.nv != 2 or PB.nu != 2 or PB.nv != 2:
        return ["coincidence certificate requires bilinear patches"], 0

    def hom(p):
        w = p[3]
        return [p[0] * w, p[1] * w, p[2] * w, w]

    def sub(a, b):
        return [a[c] - b[c] for c in range(4)]

    A = [[hom(PA.P[i][j]) for j in range(2)] for i in range(2)]
    B = [[hom(PB.P[i][j]) for j in range(2)] for i in range(2)]
    A0 = A[0][0]
    Au = sub(A[1][0], A[0][0])
    Av = sub(A[0][1], A[0][0])
    Auv = sub(sub(A[1][1], A[1][0]), sub(A[0][1], A[0][0]))
    B0 = B[0][0]
    Bs = sub(B[1][0], B[0][0])
    Bt = sub(B[0][1], B[0][0])
    Bst = sub(sub(B[1][1], B[1][0]), sub(B[0][1], B[0][0]))
    c00 = [B0[c] + ga * Bs[c] + ze * Bt[c] + ga * ze * Bst[c] for c in range(4)]
    cu = [al * Bs[c] + de * Bt[c] + (al * ze + ga * de) * Bst[c] for c in range(4)]
    cv = [be * Bs[c] + ep * Bt[c] + (be * ze + ga * ep) * Bst[c] for c in range(4)]
    cuv = [(al * ep + be * de) * Bst[c] for c in range(4)]
    # homogeneous nets may differ by one global nonzero scale
    if A0[3] == 0 or c00[3] == 0:
        errors.append("affine coincidence identity fails")
    else:
        lam = c00[3] / A0[3]
        if lam == 0:
            errors.append("affine coincidence identity fails")
        else:
            ok = True
            for got, want in ((c00, A0), (cu, Au), (cv, Av), (cuv, Auv)):
                for c in range(4):
                    if got[c] != lam * want[c]:
                        ok = False
            if not ok:
                errors.append("affine coincidence identity fails")
    cnt = 0
    for item in co["boxes"]:
        b = [parse_rat(x) for x in item["box"]]
        u0, u1, v0, v1, s0, s1, t0, t1 = b
        ok = True
        for lo, hi in ((u0, u1), (v0, v1), (s0, s1), (t0, t1)):
            if lo < 0 or hi > 1 or lo > hi:
                ok = False
        sm = [al * uu + be * vv + ga for uu in (u0, u1) for vv in (v0, v1)]
        tm = [de * uu + ep * vv + ze for uu in (u0, u1) for vv in (v0, v1)]
        if min(sm) < s0 or max(sm) > s1 or min(tm) < t0 or max(tm) > t1:
            ok = False
        if not ok:
            errors.append("coincidence box mapping fails")
        else:
            cnt += 1
    return errors, cnt


def check_tangency(cert, PA, PB):
    """精确验证相切点证书：F(p)=0 且两法向平行（叉积为零、法向非零）"""
    errors = []
    cnt = 0
    for item in cert["tangent_points"]:
        u, v, s, t = [parse_rat(x) for x in item["point"]]
        ok = True
        for x in (u, v, s, t):
            if x < 0 or x > 1:
                ok = False
        A = PA.point(u, v)
        B = PB.point(s, t)
        if A != B:
            ok = False
        nA = _cross3(PA.du_point(u, v), PA.dv_point(u, v))
        nB = _cross3(PB.du_point(s, t), PB.dv_point(s, t))
        if all(x == 0 for x in nA) or all(x == 0 for x in nB):
            ok = False
        if any(x != 0 for x in _cross3(nA, nB)):
            ok = False
        if not ok:
            errors.append("tangent point certificate fails (F or normals)")
        else:
            cnt += 1
    return errors, cnt


def _slice_claim_ok(PA, PB, U, V, S, T, iters=40):
    """验证 claimed_slice：浮点 Newton polish 到解，再对解的小盒做精确 slice Krawczyk。
    认证的小盒必须完全落在内核声称的 4D 盒内 => 内核声称盒内确有解。"""
    X0 = [U, V, S, T]
    def fp(x):
        return FI(F.from_float(float(x)), F.from_float(float(x)))
    for kfix in range(4):
        idx = [c for c in range(4) if c != kfix]
        p = [float((x.lo + x.hi) / 2) for x in X0]
        try:
            for _ in range(iters):
                u, v, s, t = p
                Au = PA.du_iv(fp(u), fp(v)); Av = PA.dv_iv(fp(u), fp(v))
                Bu = PB.du_iv(fp(s), fp(t)); Bv = PB.dv_iv(fp(s), fp(t))
                A = PA.point(F.from_float(u), F.from_float(v))
                B = PB.point(F.from_float(s), F.from_float(t))
                Fv = [float(A[c] - B[c]) for c in range(3)]
                if max(abs(x) for x in Fv) < 1e-14:
                    break
                J = [[Au[c].lo, Av[c].lo, -Bu[c].hi, -Bv[c].hi] for c in range(3)]
                J3 = [[J[r][idx[c]] for c in range(3)] for r in range(3)]
                M = [[F(J3[r][c]) for c in range(3)] for r in range(3)]
                Jinv = inv3(M)
                if Jinv is None:
                    raise ValueError("singular")
                delta = [-sum(float(Jinv[r][c]) * Fv[c] for c in range(3)) for r in range(3)]
                for k2 in range(3):
                    p[idx[k2]] += delta[k2]
                if any(p[k2] < 0.0 or p[k2] > 1.0 for k2 in range(4)):
                    raise ValueError("out of domain")
        except Exception:
            continue
        h = 1e-9 * max(1.0, max(abs(x) for x in p))
        box = []
        good = True
        for k2 in range(4):
            lo = max(F.from_float(p[k2] - h), X0[k2].lo)
            hi = min(F.from_float(p[k2] + h), X0[k2].hi)
            if lo > hi or lo < 0 or hi > 1:
                good = False
                break
            box.append(FI(lo, hi))
        if not good:
            continue
        if _krawczyk_ok_single(PA, PB, box[0], box[1], box[2], box[3]):
            return True
    return False


def _polyline_dist2(p, polylines):
    best = None
    for pl in polylines:
        for i in range(len(pl) - 1):
            d2 = _seg_dist2(p, pl[i], pl[i + 1])
            if best is None or d2 < best:
                best = d2
    return best


def _exact_solution_points(cert, PA, PB):
    pts, errors = [], []
    for item in cert.get("solution_points", []):
        u, v, s, t = [parse_rat(x) for x in item["point"]]
        if any(x < 0 or x > 1 for x in (u, v, s, t)):
            errors.append("solution point outside [0,1]^4")
            continue
        A = PA.point(u, v)
        B = PB.point(s, t)
        if A != B:
            errors.append("solution point fails exact F=0")
        else:
            pts.append(A)
    return pts, errors


def check_miss(cert, PA, PB):
    """严格遗漏证书：精确交点 S(p) 到 OCCT 折线的平方距离 > delta²（全部有理数精确）"""
    errors = []
    polylines = [[[parse_rat(c) for c in q] for q in pl] for pl in cert["occt_polyline"]]
    delta2 = parse_rat(cert["delta2"])
    if delta2 < 0:
        errors.append("delta2 must be non-negative")
    pts, errs = _exact_solution_points(cert, PA, PB)
    errors.extend(errs)
    cnt = 0
    for A in pts:
        d2 = _polyline_dist2(A, polylines)
        if d2 is None:
            errors.append("empty OCCT polyline")
        elif not (d2 > delta2):
            errors.append("miss claim fails: dist^2 <= delta^2")
        else:
            cnt += 1
    return errors, cnt


def check_coverage(cert, PA, PB):
    """严格覆盖证书：精确交点样本到 OCCT 折线的平方距离 <= eps²"""
    errors = []
    polylines = [[[parse_rat(c) for c in q] for q in pl] for pl in cert["occt_polyline"]]
    eps2 = parse_rat(cert["eps2"])
    if eps2 < 0:
        errors.append("eps2 must be non-negative")
    pts, errs = _exact_solution_points(cert, PA, PB)
    errors.extend(errs)
    cnt = 0
    for A in pts:
        d2 = _polyline_dist2(A, polylines)
        if d2 is None:
            errors.append("empty OCCT polyline")
        elif d2 > eps2:
            errors.append("coverage claim fails: dist^2 > eps^2")
        else:
            cnt += 1
    return errors, cnt


def _max_corner_dist2(enc, a, b):
    """3D 区间盒 8 顶点到线段 ab 的最大平方距离（精确有理数）"""
    best = None
    for sx in (enc[0].lo, enc[0].hi):
        for sy in (enc[1].lo, enc[1].hi):
            for sz in (enc[2].lo, enc[2].hi):
                d2 = _seg_dist2([sx, sy, sz], a, b)
                if best is None or d2 > best:
                    best = d2
    return best


def _corner_dist2_float(enc, a, b):
    best = 0.0
    for sx in (float(enc[0].lo), float(enc[0].hi)):
        for sy in (float(enc[1].lo), float(enc[1].hi)):
            for sz in (float(enc[2].lo), float(enc[2].hi)):
                dx, dy, dz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
                L2 = dx * dx + dy * dy + dz * dz
                if L2 == 0:
                    d2 = (sx - a[0]) ** 2 + (sy - a[1]) ** 2 + (sz - a[2]) ** 2
                else:
                    p = [sx, sy, sz]
                    dot = (p[0] - a[0]) * dx + (p[1] - a[1]) * dy + (p[2] - a[2]) * dz
                    if dot <= 0:
                        d2 = (p[0] - a[0]) ** 2 + (p[1] - a[1]) ** 2 + (p[2] - a[2]) ** 2
                    elif dot >= L2:
                        d2 = (p[0] - b[0]) ** 2 + (p[1] - b[1]) ** 2 + (p[2] - b[2]) ** 2
                    else:
                        d2 = (p[0] - a[0]) ** 2 + (p[1] - a[1]) ** 2 + (p[2] - a[2]) ** 2 - dot * dot / L2
                if d2 > best:
                    best = d2
    return best


def check_branch_coverage(cert, PA, PB):
    """整支严格覆盖：对每个输出叶，其 3D 区间包络到 OCCT 折线的最大距离 <= eps（精确有理数）"""
    errors = []
    polylines = [[[parse_rat(c) for c in q] for q in pl] for pl in cert["occt_polyline"]]
    segs = []
    for pl in polylines:
        for i in range(len(pl) - 1):
            segs.append((pl[i], pl[i + 1]))
    eps2 = parse_rat(cert["eps2"])
    leaves = []

    def walk(n):
        if n.get("kind") == "leaf":
            if n["status"] == "output":
                leaves.append(n)
            return
        for ch in n.get("children", []):
            walk(ch)

    walk(cert["tree"])
    cnt = 0
    fails = 0
    for n in leaves:
        box = [FI(parse_rat(n["box"][2 * i]), parse_rat(n["box"][2 * i + 1])) for i in range(4)]
        try:
            enc = PA.cart_iv(box[0], box[1])
        except ZeroDivisionError:
            errors.append("branch coverage: undefined enclosure")
            fails += 1
            continue
        # 浮点预筛：选使 max corner dist2 最小的线段（任一固定线段都给出上界，故 sound）
        best = None
        for a, b in segs:
            v = _corner_dist2_float(enc, a, b)
            if best is None or v < best[0]:
                best = (v, a, b)
        if best is None:
            errors.append("branch coverage: empty polyline")
            fails += 1
            continue
        exact = _max_corner_dist2(enc, best[1], best[2])
        if exact <= eps2:
            cnt += 1
        else:
            fails += 1
            errors.append("branch coverage fails in a leaf")
    return errors, cnt, fails


def _canonical_polyline_digest(cert):
    poly = cert.get("occt_polyline")
    canon = json.dumps(poly, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canon.encode()).hexdigest()


def check_claims(cert):
    """Optional claims block: bind the audited output polyline."""
    if "claims" not in cert:
        return []
    errors = []
    claims = cert["claims"]
    if not isinstance(claims, dict):
        return ["claims must be an object"]
    if claims.get("format", "polyline3d-v1") != "polyline3d-v1":
        errors.append("unsupported claims format")
    if "occt_polyline" not in cert:
        errors.append("claims present but no occt_polyline")
        return errors
    digest = claims.get("output_digest")
    if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
        errors.append("output_digest must be 64 lowercase hex characters")
        return errors
    actual = _canonical_polyline_digest(cert)
    if actual != digest:
        errors.append("output digest mismatch")
    return errors


def verify(cert, verbose=True):
    PA = Patch(cert["patches"][0])
    PB = Patch(cert["patches"][1])
    errors = []
    stats = {"nodes": 0, "excluded": 0, "output": 0, "cert_ok": 0, "cert_fail": 0, "p7_ok": 0, "k4_ok": 0, "k4_fail": 0}
    claims_errors = check_claims(cert)
    if claims_errors:
        errors.extend(claims_errors)
        return False, stats, errors
    if not PA.positive or not PB.positive:
        errors.append("patch weights must be positive")
        if verbose:
            print("FINAL: REJECT")
            for e in errors:
                print("  ERROR:", e)
        return False, stats, errors
    if "chunks" in cert:
        errors.append("chunked certificate: use --chunked")
        if verbose:
            print("FINAL: REJECT")
            for e in errors:
                print("  ERROR:", e)
        return False, stats, errors

    # ---- P7 证书：独立重做 NURBS->Bezier（精确有理数 Boehm），与声称块比对 ----
    if cert.get("nurbs"):
        for idx, nur in enumerate(cert["nurbs"]):
            try:
                blocks = p7_convert(nur)
            except Exception as e:
                errors.append("P7 conversion failed: %r" % (e,))
                continue
            claimed = [[parse_rat(x) for x in p] for p in cert["patches"][idx]["poles"]]
            claimed_nu = cert["patches"][idx]["nu"]
            claimed_nv = cert["patches"][idx]["nv"]
            if claimed_nu != nur["du"] + 1 or claimed_nv != nur["dv"] + 1:
                errors.append("P7 claimed block degrees mismatch")
                continue
            matched = any(all(claimed[k] == blk[k] for k in range(len(claimed))) for blk in blocks)
            if matched:
                stats["p7_ok"] += 1
            else:
                errors.append("P7 mismatch: claimed block not in exact conversion")

    def box_of(b):
        return [FI(parse_rat(b[2 * i]), parse_rat(b[2 * i + 1])) for i in range(4)]

    def rec(node, A, B):
        stats["nodes"] += 1
        box = box_of(node["box"])
        for c in range(4):
            if box[c].lo < 0 or box[c].hi > 1 or box[c].lo > box[c].hi:
                errors.append("box outside [0,1]^4 or invalid")
        if node["kind"] == "leaf":
            status = node.get("status")
            if status not in ("excluded", "output"):
                errors.append("unknown leaf status")
                return
            if status == "excluded":
                stats["excluded"] += 1
                la, ha = aabb(A)
                lb, hb = aabb(B)
                if not disjoint(la, ha, lb, hb):
                    errors.append("excluded leaf is NOT disjoint (exact AABB)")
            else:
                stats["output"] += 1
                if node.get("claimed_4d"):
                    if krawczyk4d_ok(PA, PB, box[0], box[1], box[2], box[3]):
                        stats["k4_ok"] += 1
                    else:
                        stats["k4_fail"] += 1
                        errors.append("claimed_4d fails exact 4D Krawczyk")
                if node.get("claimed_slice"):
                    if _slice_claim_ok(PA, PB, box[0], box[1], box[2], box[3]):
                        stats["cert_ok"] += 1
                    else:
                        stats["cert_fail"] += 1
                        errors.append("claimed_slice fails exact slice Krawczyk")
            return
        axis = node["axis"]
        ch = node["children"]
        if axis not in (0, 1, 2, 3) or len(ch) != 2:
            errors.append("invalid node")
            return
        mid = (box[axis].lo + box[axis].hi) / 2
        exp_lo = list(box); exp_lo[axis] = FI(box[axis].lo, mid)
        exp_hi = list(box); exp_hi[axis] = FI(mid, box[axis].hi)
        for k, exp in ((0, exp_lo), (1, exp_hi)):
            cb = box_of(ch[k]["box"])
            if any(cb[c].lo != exp[c].lo or cb[c].hi != exp[c].hi for c in range(4)):
                errors.append("children do not exactly partition parent")
        if axis == 0:
            L, R = split_u(A)
            rec(ch[0], L, B); rec(ch[1], R, B)
        elif axis == 1:
            L, R = split_v(A)
            rec(ch[0], L, B); rec(ch[1], R, B)
        elif axis == 2:
            L, R = split_u(B)
            rec(ch[0], A, L); rec(ch[1], A, R)
        else:
            L, R = split_v(B)
            rec(ch[0], A, L); rec(ch[1], A, R)

    if "tree" in cert:
        root_b = box_of(cert["tree"]["box"])
        if any(root_b[c].lo != 0 or root_b[c].hi != 1 for c in range(4)):
            errors.append("tree root is not [0,1]^4")
        rec(cert["tree"],
            {"nu": PA.nu, "nv": PA.nv, "P": PA.P},
            {"nu": PB.nu, "nv": PB.nv, "P": PB.P})
    for item in cert.get("certified", []):
        b = box_of(item["box"])
        for c in range(4):
            if b[c].lo < 0 or b[c].hi > 1 or b[c].lo > b[c].hi:
                errors.append("certified box outside [0,1]^4 or invalid")
        if krawczyk_ok(PA, PB, b[0], b[1], b[2], b[3]):
            stats["cert_ok"] += 1
        else:
            stats["cert_fail"] += 1
            errors.append("certified box fails exact Krawczyk")
    if "coincidence" in cert:
        errs, cnt = check_coincidence(cert, PA, PB)
        errors.extend(errs)
        stats["coinc_ok"] = stats.get("coinc_ok", 0) + cnt
    if "tangent_points" in cert:
        errs, cnt = check_tangency(cert, PA, PB)
        errors.extend(errs)
        stats["tan_ok"] = stats.get("tan_ok", 0) + cnt
    if "occt_polyline" in cert and "solution_points" in cert:
        if "delta2" in cert:
            errs, cnt = check_miss(cert, PA, PB)
            stats["miss_ok"] = stats.get("miss_ok", 0) + cnt
        else:
            errs, cnt = check_coverage(cert, PA, PB)
            stats["cover_ok"] = stats.get("cover_ok", 0) + cnt
        errors.extend(errs)
    if cert.get("kind") == "branch_coverage":
        errs, cnt, fails = check_branch_coverage(cert, PA, PB)
        stats["branch_ok"] = stats.get("branch_ok", 0) + cnt
        stats["branch_fail"] = stats.get("branch_fail", 0) + fails
        errors.extend(errs)
    ok = len(errors) == 0
    if verbose:
        print("nodes=%d excluded=%d output=%d cert_ok=%d cert_fail=%d p7_ok=%d k4_ok=%d k4_fail=%d coinc_ok=%d tan_ok=%d miss_ok=%d cover_ok=%d errors=%d" % (
            stats["nodes"], stats["excluded"], stats["output"],
            stats["cert_ok"], stats["cert_fail"], stats["p7_ok"],
            stats["k4_ok"], stats["k4_fail"],
            stats.get("coinc_ok", 0), stats.get("tan_ok", 0), stats.get("miss_ok", 0), stats.get("cover_ok", 0), len(errors)))
        for e in errors[:10]:
            print("  ERROR:", e)
    return ok, stats, errors


def _reconstruct(PA, PB, path):
    A = {"nu": PA.nu, "nv": PA.nv, "P": PA.P}
    B = {"nu": PB.nu, "nv": PB.nv, "P": PB.P}
    for axis, side in path:
        if axis == 0:
            L, R = split_u(A); A = L if side == 0 else R
        elif axis == 1:
            L, R = split_v(A); A = L if side == 0 else R
        elif axis == 2:
            L, R = split_u(B); B = L if side == 0 else R
        else:
            L, R = split_v(B); B = L if side == 0 else R
    return A, B


def _box_for_path(path):
    box = [FI(F(0), F(1)) for _ in range(4)]
    for axis, side in path:
        mid = (box[axis].lo + box[axis].hi) / 2
        box[axis] = FI(box[axis].lo, mid) if side == 0 else FI(mid, box[axis].hi)
    return box


def _chunk_verify(chunk, box, PA, PB):
    import hashlib as _h
    import json as _j
    stats = {"nodes": 0, "excluded": 0, "output": 0, "k4_ok": 0, "k4_fail": 0}
    errors = []
    canon = _j.dumps(chunk["tree"], sort_keys=True, separators=(",", ":"))
    if _h.sha256(canon.encode()).hexdigest() != chunk["sha256"]:
        errors.append("chunk sha256 mismatch")
    rb = [parse_rat(x) for x in chunk["root"]]
    for c in range(4):
        if rb[2 * c] != box[c].lo or rb[2 * c + 1] != box[c].hi:
            errors.append("chunk root box mismatch")
    A, B = _reconstruct(PA, PB, chunk["path"])

    def rec(node, A, B):
        stats["nodes"] += 1
        b = [FI(parse_rat(node["box"][2 * i]), parse_rat(node["box"][2 * i + 1])) for i in range(4)]
        for c in range(4):
            if b[c].lo < 0 or b[c].hi > 1 or b[c].lo > b[c].hi:
                errors.append("chunk node box outside [0,1]^4 or invalid")
        if node["kind"] == "leaf":
            status = node.get("status")
            if status not in ("excluded", "output"):
                errors.append("unknown leaf status")
                return
            if status == "excluded":
                stats["excluded"] += 1
                la, ha = aabb(A)
                lb, hb = aabb(B)
                if not disjoint(la, ha, lb, hb):
                    errors.append("excluded leaf is NOT disjoint (exact AABB)")
            else:
                stats["output"] += 1
                if node.get("claimed_4d"):
                    if krawczyk4d_ok(PA, PB, b[0], b[1], b[2], b[3]):
                        stats["k4_ok"] += 1
                    else:
                        stats["k4_fail"] += 1
                        errors.append("claimed_4d fails exact 4D Krawczyk")
                if node.get("claimed_slice"):
                    if _slice_claim_ok(PA, PB, b[0], b[1], b[2], b[3]):
                        stats["cert_ok"] = stats.get("cert_ok", 0) + 1
                    else:
                        stats["cert_fail"] = stats.get("cert_fail", 0) + 1
                        errors.append("claimed_slice fails exact slice Krawczyk")
            return
        axis = node["axis"]
        ch = node["children"]
        if axis not in (0, 1, 2, 3) or len(ch) != 2:
            errors.append("invalid chunk node")
            return
        mid = (b[axis].lo + b[axis].hi) / 2
        exp_lo = list(b); exp_lo[axis] = FI(b[axis].lo, mid)
        exp_hi = list(b); exp_hi[axis] = FI(mid, b[axis].hi)
        for k, exp in ((0, exp_lo), (1, exp_hi)):
            cb = [FI(parse_rat(ch[k]["box"][2 * i]), parse_rat(ch[k]["box"][2 * i + 1])) for i in range(4)]
            if any(cb[c].lo != exp[c].lo or cb[c].hi != exp[c].hi for c in range(4)):
                errors.append("children do not exactly partition parent")
        if axis == 0:
            L, R = split_u(A); rec(ch[0], L, B); rec(ch[1], R, B)
        elif axis == 1:
            L, R = split_v(A); rec(ch[0], L, B); rec(ch[1], R, B)
        elif axis == 2:
            L, R = split_u(B); rec(ch[0], A, L); rec(ch[1], A, R)
        else:
            L, R = split_v(B); rec(ch[0], A, L); rec(ch[1], A, R)

    rec(chunk["tree"], A, B)
    return stats, errors


def verify_chunked(cert, jobs=1, verbose=True):
    PA = Patch(cert["patches"][0])
    PB = Patch(cert["patches"][1])
    chunks = cert["chunks"]
    chunkmap = {tuple(tuple(p) for p in c["path"]): c for c in chunks}
    errors = []
    partition_errors = []
    if len(chunkmap) != len(chunks):
        partition_errors.append("duplicate chunk path")
    for c in chunks:
        for axis, side in c["path"]:
            if axis not in (0, 1, 2, 3) or side not in (0, 1):
                partition_errors.append("invalid chunk path")
    visited = set()

    def walk(prefix, box):
        if prefix in chunkmap:
            visited.add(prefix)
            c = chunkmap[prefix]
            rb = [parse_rat(x) for x in c["root"]]
            for k in range(4):
                if rb[2 * k] != box[k].lo or rb[2 * k + 1] != box[k].hi:
                    partition_errors.append("chunk root does not match frontier box")
            return
        cands = [p for p in chunkmap if p[:len(prefix)] == prefix]
        if not cands:
            partition_errors.append("frontier hole at %r" % (prefix,))
            return
        axis = cands[0][len(prefix)][0]
        mid = (box[axis].lo + box[axis].hi) / 2
        left = list(box); left[axis] = FI(box[axis].lo, mid)
        right = list(box); right[axis] = FI(mid, box[axis].hi)
        walk(prefix + ((axis, 0),), left)
        walk(prefix + ((axis, 1),), right)

    walk((), [FI(F(0), F(1)) for _ in range(4)])
    if len(visited) != len(chunks):
        partition_errors.append("unused chunks")

    results = []
    if jobs > 1:
        import multiprocessing as mp
        tasks = [(c, _box_for_path(c["path"]), PA, PB) for c in chunks]
        ctx = mp.get_context("fork")
        with ctx.Pool(jobs) as pool:
            results = pool.starmap(_chunk_verify, tasks)
    else:
        for c in chunks:
            results.append(_chunk_verify(c, _box_for_path(c["path"]), PA, PB))

    stats = {"nodes": 0, "excluded": 0, "output": 0, "k4_ok": 0, "k4_fail": 0}
    for st, errs in results:
        for k in stats:
            stats[k] += st[k]
        errors.extend(errs)
    errors.extend(partition_errors)
    errors.extend(check_claims(cert))

    # top-level fields that are not part of the chunk subtrees
    if cert.get("nurbs"):
        for idx, nur in enumerate(cert["nurbs"]):
            try:
                blocks = p7_convert(nur)
            except Exception as e:
                errors.append("P7 conversion failed: %r" % (e,))
                continue
            claimed = [[parse_rat(x) for x in p] for p in cert["patches"][idx]["poles"]]
            claimed_nu = cert["patches"][idx]["nu"]
            claimed_nv = cert["patches"][idx]["nv"]
            if claimed_nu != nur["du"] + 1 or claimed_nv != nur["dv"] + 1:
                errors.append("P7 claimed block degrees mismatch")
                continue
            matched = any(all(claimed[k] == blk[k] for k in range(len(claimed))) for blk in blocks)
            if matched:
                stats["p7_ok"] = stats.get("p7_ok", 0) + 1
            else:
                errors.append("P7 mismatch: claimed block not in exact conversion")
    for item in cert.get("certified", []):
        b = [FI(parse_rat(item["box"][2 * i]), parse_rat(item["box"][2 * i + 1])) for i in range(4)]
        if any(b[c].lo < 0 or b[c].hi > 1 or b[c].lo > b[c].hi for c in range(4)):
            errors.append("certified box outside [0,1]^4 or invalid")
        if krawczyk_ok(PA, PB, b[0], b[1], b[2], b[3]):
            stats["cert_ok"] = stats.get("cert_ok", 0) + 1
        else:
            stats["cert_fail"] = stats.get("cert_fail", 0) + 1
            errors.append("certified box fails exact Krawczyk")
    if "coincidence" in cert:
        errs, cnt = check_coincidence(cert, PA, PB)
        errors.extend(errs)
        stats["coinc_ok"] = stats.get("coinc_ok", 0) + cnt
    if "tangent_points" in cert:
        errs, cnt = check_tangency(cert, PA, PB)
        errors.extend(errs)
        stats["tan_ok"] = stats.get("tan_ok", 0) + cnt
    if "occt_polyline" in cert and "solution_points" in cert:
        if "delta2" in cert:
            errs, cnt = check_miss(cert, PA, PB)
            stats["miss_ok"] = stats.get("miss_ok", 0) + cnt
        else:
            errs, cnt = check_coverage(cert, PA, PB)
            stats["cover_ok"] = stats.get("cover_ok", 0) + cnt
        errors.extend(errs)
    if cert.get("kind") == "branch_coverage":
        errors.append("branch coverage certificate cannot be chunked")

    ok = len(errors) == 0
    if verbose:
        print("chunked: chunks=%d jobs=%d nodes=%d excluded=%d output=%d k4_ok=%d k4_fail=%d errors=%d" % (
            len(chunks), jobs, stats["nodes"], stats["excluded"], stats["output"],
            stats["k4_ok"], stats["k4_fail"], len(errors)))
        for e in errors[:10]:
            print("  ERROR:", e)
    return ok, stats, errors


def main():
    if len(sys.argv) < 2:
        print("usage: python3 ssi_verifier.py cert.json [--chunked] [--jobs N]")
        return 2
    cert = json.load(open(sys.argv[1]))
    if "chunks" in cert and "--chunked" not in sys.argv:
        print("note: certificate has chunks; using chunked verification")
    try:
        if "--chunked" in sys.argv or "chunks" in cert:
            jobs = 1
            if "--jobs" in sys.argv:
                jobs = int(sys.argv[sys.argv.index("--jobs") + 1])
            ok, stats, errors = verify_chunked(cert, jobs=jobs)
        else:
            ok, stats, errors = verify(cert)
    except Exception as e:
        print("FINAL: REJECT")
        print("  ERROR:", e)
        return 1
    print("FINAL:", "ACCEPT" if ok else "REJECT")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())















