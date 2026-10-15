# =============================================================================
# example_ctypes.py — Python 经 ctypes 调用 NURBS 算法库示例
# 运行：python example_ctypes.py   （自动加载 cpp/build/nurbs_shared.dll）
# =============================================================================
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nurbs_ctypes as nb

lib = nb.load()
nb.set_default_lib(lib)
print("[py] 已加载 DLL:", lib._name if hasattr(lib, "_name") else "nurbs_shared.dll")

# ---- 1) 构造曲线 ----
L1 = nb.Curve.from_points([(0, 0), (4, 4)], 1)           # y=x
L2 = nb.Curve.from_points([(0, 4), (4, 0)], 1)           # y=4-x  → 交 (2,2)
Q = nb.Curve.from_points([(0, 0), (1, 2), (2, 0)], 2)    # 二次样条
C = nb.Curve.from_points([(0, 0), (1, 3), (3, -1), (4, 1)], 3)

# ---- 2) 求值（与解析解核对） ----
print("\n=== 求值 ===")
print(f"  Q(0.0) = {Q.evaluate(0.0)}  期望 (0,0)")
print(f"  Q(0.5) = {Q.evaluate(0.5)}  期望 (1,1)")
print(f"  Q(1.0) = {Q.evaluate(1.0)}  期望 (2,0)")
print(f"  Q'(0.5) = {Q.derivative(0.5)}  期望 (2,0)")

# ---- 3) 求交 / 最近点 ----
print("\n=== unified_solve（求交 / 最近点）===")
for nm, a, b in [("L1×L2", L1, L2), ("Q×C", Q, C)]:
    r = nb.best_solve(a.unified_solve(b))
    print(f"  {nm:6s} 状态={nb.STATE_NAMES[r.state]}  u={r.u:.6f} v={r.v:.6f}  "
          f"P1=({r.point1_x:.4f},{r.point1_y:.4f}) d={r.distance:.3e} iter={r.iterations}")

# ---- 4) 全交点 ----
print("\n=== findAllIntersections ===")
for nm, a, b in [("L1×L2", L1, L2)]:
    rs = a.find_intersections(b)
    print(f"  {nm:6s} 根数={len(rs)}")
    for i, r in enumerate(rs):
        print(f"    root#{i} u={r.u:.9f} v={r.v:.9f} d={r.distance:.3e} P=({r.point1_x:.6f},{r.point1_y:.6f})")

# ---- 5) 错误处理（负权重曲线应报错） ----
print("\n=== 错误处理 ===")
try:
    nb.Curve.from_points([(0, 0), (1, 1)], 1, weights=[1, -1])
    print("  !! 未拦截非法曲线（不应发生）")
except ValueError as e:
    print(f"  非法曲线(负权重)被拦截: {e}")

for c in (L1, L2, Q, C):
    c.destroy()
print("\n[py] ctypes 示例完成")
