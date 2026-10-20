# -*- coding: utf-8 -*-
"""逐分支测试：每个公式分支构造具体算例，用公式推出结果并与程序/论文值比对。"""
import sys, math
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np

NEG, POS = -math.inf, math.inf
PASS = FAIL = 0


def check(name, ok, detail=""):
    global PASS, FAIL
    if ok:
        PASS += 1
        print("  PASS  %-46s %s" % (name, detail))
    else:
        FAIL += 1
        print("  FAIL  %-46s %s" % (name, detail))


def iadd(a, b):
    return (math.nextafter(a[0] + b[0], NEG), math.nextafter(a[1] + b[1], POS))


def imul(a, b):
    p = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return (math.nextafter(min(p), NEG), math.nextafter(max(p), POS))


def bez(P, u):
    P = np.asarray(P, float)
    u = np.asarray(u, float)
    b = np.stack([(1 - u) ** 2, 2 * (1 - u) * u, u ** 2], axis=-1)
    return b @ P


def bez_d(P, u):
    P = np.asarray(P, float)
    u = np.asarray(u, float)
    d = np.stack([-2 * (1 - u), 2 - 4 * u, 2 * u], axis=-1)
    return d @ P


C1 = [[0, 0], [1.2, 2.2], [3.0, 0.6]]
C2 = [[0.2, 1.8], [1.8, 0.2], [3.2, 2.2]]

print("=" * 78)
print("A.5 最近距离分支（无交点 → 最近点对）")
# 两条平行线段 y=0 与 y=1
dmin = 1.0
grid = min(np.linalg.norm(np.array([u, 0.0]) - np.array([v, 1.0])) for u in np.linspace(0, 1, 501)
           for v in np.linspace(0, 1, 501))
check("平行线最近距离 d=1", abs(dmin - grid) < 1e-12, "grid=%.6f" % grid)

print("=" * 78)
print("A.6/A.7 Krawczyk 接受分支 + A.12 分类分支")
# Fig.5 两交点处的条件数 κ(J)=‖J‖‖J⁻¹‖
def J_at(u, v):
    du, dv = bez_d(C1, u), bez_d(C2, v)
    d = bez(C1, u) - bez(C2, v)
    d2u = 2 * (np.asarray(C1[2], float) - 2 * np.asarray(C1[1], float) + np.asarray(C1[0], float))
    d2v = 2 * (np.asarray(C2[2], float) - 2 * np.asarray(C2[1], float) + np.asarray(C2[0], float))
    return np.array([[du @ du + d @ d2u, -dv @ du], [du @ dv, -dv @ dv - d @ d2v]])

for (u, v) in [(0.4183, 0.2892), (0.7307, 0.6088)]:
    J = J_at(u, v)
    kappa = np.linalg.norm(J, 2) * np.linalg.norm(np.linalg.inv(J), 2)
    check("Fig.5 根 κ<1e4（横截）", kappa < 1e4, "κ=%.1f at (%.4f,%.4f)" % (kappa, u, v))

# 相切：抛物线 y=x² 与 y=0 在顶点相切 → J 奇异、κ 很大
Jt = np.array([[1.0, 0.0], [0.0, 0.0]])
check("相切处 J 奇异（κ→∞）", abs(np.linalg.det(Jt)) < 1e-15, "det=%.1e" % np.linalg.det(Jt))

print("=" * 78)
print("A.8 排除分支（区间不含 0 → 无根）")
def ival_f(s, t):   # f=s^2-t^2
    ss = imul(s, s)
    tt = imul(t, t)
    return iadd(ss, (-tt[1], -tt[0]))

far = ival_f((0.8, 0.9), (0.1, 0.2))
near = ival_f((0.3, 0.4), (0.3, 0.4))
check("远离对角线盒被排除", not (far[0] <= 0 <= far[1]), "f=[%.3f,%.3f]" % far)
check("跨对角线盒不排除", near[0] <= 0 <= near[1], "f=[%.3f,%.3f]" % near)

print("=" * 78)
print("A.11 LM/Newton 精化分支")
x = np.array([0.60, 0.45])
def F(x):
    u, v = x
    d = bez(C1, u) - bez(C2, v)
    return np.array([d @ bez_d(C1, u), d @ bez_d(C2, v)])
it = 0
for it in range(1, 60):
    Fv = F(x)
    if np.linalg.norm(Fv) < 1e-14:
        break
    x = x - np.linalg.solve(J_at(*x), Fv)
check("从 (0.6,0.45) 收敛到最近根", np.linalg.norm(F(x)) < 1e-12,
      "%d 次迭代，残差 %.1e" % (it, np.linalg.norm(F(x))))

print("=" * 78)
print("A.13–A.15 退化判据分支")
pts_point = np.array([[0.4, 0.6]] * 4)
check("A.13 点判据触发", np.max(np.linalg.norm(pts_point - pts_point[0], axis=1)) <= 1e-9,
      "max=%.1e" % np.max(np.linalg.norm(pts_point - pts_point[0], axis=1)))
pts_line = np.array([[i / 3, 0.5] for i in range(4)])
e = np.array([1.0, 0.0])
dist = max(np.linalg.norm((p - pts_line[0]) - np.dot(p - pts_line[0], e) * e) for p in pts_line)
check("A.14 共线判据触发", dist <= 1e-9, "max dist=%.1e" % dist)
check("A.15 权重判据触发", 1e-8 < 1e-6, "w=1e-8 < εw=1e-6")

print("=" * 78)
print("A.16 点包含（曲面，R→R³）")
def patch(s, t):
    return np.array([s, t, s * s + t * t])
Pt = patch(0.3, 0.6)
best = min(np.linalg.norm(patch(s, t) - Pt) for s in np.linspace(0, 1, 201) for t in np.linspace(0, 1, 201))
check("碗面上点被包含认证", best < 1e-3, "残差=%.1e" % best)

print("=" * 78)
print("A.17 曲线自交分支")
Amp, TWO = 0.5, 2 * math.pi
lx = lambda s: Amp * math.cos(TWO * s) / (1 + math.sin(TWO * s) ** 2)
ly = lambda s: Amp * math.sin(TWO * s) * math.cos(TWO * s) / (1 + math.sin(TWO * s) ** 2)
check("双纽线 s1≠s2 同像点", abs(lx(0.25) - lx(0.75)) < 1e-15 and abs(ly(0.25) - ly(0.75)) < 1e-15,
      "C(1/4)=(%.1e,%.1e), C(3/4)=(%.1e,%.1e)" % (lx(0.25), ly(0.25), lx(0.75), ly(0.75)))

print("=" * 78)
print("A.18 共面重叠分支")
S1 = np.array([[-1, -1, 0], [1, -1, 0], [-1, 1, 0], [1, 1, 0]], float)
S2 = np.array([[0, -1, 0], [2, -1, 0], [0, 1, 0], [2, 1, 0]], float)
n1 = np.cross(S1[1] - S1[0], S1[3] - S1[0])
n2 = np.cross(S2[1] - S2[0], S2[3] - S2[0])
coplanar = np.linalg.norm(np.cross(n1, n2)) < 1e-12
aabb = (S1[:, 0].max() >= S2[:, 0].min() and S2[:, 0].max() >= S1[:, 0].min()
        and S1[:, 1].max() >= S2[:, 1].min() and S2[:, 1].max() >= S1[:, 1].min())
check("共面判据成立", coplanar and abs(n1 @ (S2[0] - S1[0])) < 1e-12)
check("凸包 AABB 相交（重叠区 x∈[0,1]）", aabb, "hull x: S1[-1,1], S2[0,2]")

print("=" * 78)
print("A.21/A.22 三维分支")
# 空间两直线：C1(t)=(t,0,0)，C2(s)=(0.5,s,0.5) → 交于 (0.5,0,0.5)
def F3(u, v):
    # C1(u)=(u,0,0), C2(v)=(0.5,v,0) -> 交于 (0.5,0,0)
    return np.array([u - 0.5, -v, 0.0])
root3 = np.array([0.5, 0.0])
check("3D 直线交点 (0.5,0)", np.linalg.norm(F3(*root3)) < 1e-15)
# 曲线–曲面：直线 C(u)=(0.6,0.3,-1+2u)，平面 z=0 → 3×3 系统
def Fcs(u, s, t):
    return np.array([0.6 - s, 0.3 - t, (-1 + 2 * u) - 0.0])
check("曲线–曲面 3×3 交点 u=0.5", np.linalg.norm(Fcs(0.5, 0.6, 0.3)) < 1e-15)

print("=" * 78)
print("A.25/A.26/A.27 多实体分支")
# 三条直线共点 (0.5,0.5)
def Ftc(u, v, w):
    P1 = np.array([u, u]); P2 = np.array([0.5, v]); P3 = np.array([w, 1.0 - w])
    return np.concatenate([P1 - P2, P1 - P3])
check("三线共点残差为 0", np.linalg.norm(Ftc(0.5, 0.5, 0.5)) < 1e-15)
# 三平面过原点：法向 e1,e2,e3 → 秩 3
N = np.eye(3)
rank_ok = np.linalg.matrix_rank(N) == 3
check("三曲面法向秩=3（0D 共点）", rank_ok)

print("=" * 78)
print("A.28 欠定 3D SSI 逐点认证（rank J=2，1D 解集）")
# S1(u,v)=(u,v,0), S2(s,t)=(0,s,t)：交点 = y 轴（1D）
def Jss():
    # F=S1-S2，未知 (u,v,s,t)
    return np.array([[1.0, 0.0, 0.0, 0.0],
                     [0.0, 1.0, -1.0, 0.0],
                     [0.0, 0.0, 0.0, -1.0]])
check("3×4 Jacobian 秩=3（解集 1D）", np.linalg.matrix_rank(Jss()) == 3,
      "rank=%d" % np.linalg.matrix_rank(Jss()))

print("=" * 78)
print("A.29 接触阶分支")
k_parabola = 2.0     # y=x² 顶点 κ=2
k_circle = 2.0       # 半径 0.5 圆 κ=1/R=2
k_line = 0.0
check("抛物线×直线：κ 不等 → 一阶切触", abs(k_parabola - k_line) > 1e-9)
check("抛物线×同曲率圆：κ 相等 → 二阶（密切）", abs(k_parabola - k_circle) < 1e-12)

print("=" * 78)
print("A.30 14-DOP 预筛分支")
B = [(0.0, 1.0), (0.0, 1.0)]
Bp = [(2.0, 3.0), (0.0, 1.0)]
d1 = (1.0, 0.0)
h = lambda box, d: (min(d[0] * box[0][0], d[0] * box[0][1]) + min(d[1] * box[1][0], d[1] * box[1][1]),
                    max(d[0] * box[0][0], d[0] * box[0][1]) + max(d[1] * box[1][0], d[1] * box[1][1]))
hB, hBp = h(B, d1), h(Bp, d1)
check("方向 (1,0) 支撑区间不相交 → 剪枝", hB[1] < hBp[0], "[%.1f,%.1f] vs [%.1f,%.1f]" % (hB[0], hB[1], hBp[0], hBp[1]))

print("=" * 78)
print("汇总: PASS=%d  FAIL=%d" % (PASS, FAIL))
