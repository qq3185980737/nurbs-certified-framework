# =============================================================================
# nurbs_ctypes.py — 通过 Python ctypes 加载 nurbs_shared.dll 的绑定模块
# -----------------------------------------------------------------------------
# 用法：
#   import nurbs_ctypes as nb
#   lib = nb.load()                 # 或 nb.load("路径/nurbs_shared.dll")
#   c1 = nb.Curve.from_points([(0,0),(2,2),(4,4)...], 1, weights=...)
#   res = c1.unified_solve(c2)      # 返回 NurbsSolveResult 列表
#   roots = c1.find_intersections(c2)
# 数值与 JS / C++ 核心逐位一致（见 cpp/tests/test_equivalence）。
# =============================================================================
import ctypes
import os

_NURBS = "nurbs_shared.dll"
STATE_INTERSECTION = 0
STATE_TANGENT = 1
STATE_SEPARATED = 2
STATE_COINCIDENT = 3
STATE_UNDETERMINED = 4
STATE_NAMES = {
    0: "INTERSECTION 相交", 1: "TANGENT 相切", 2: "SEPARATED 分离(最近点)",
    3: "COINCIDENT 重合", 4: "UNDETERMINED 无法判定",
}
OK = 0
ERR_INVALID = -1
ERR_RUNTIME = -2


class NurbsControlPoint(ctypes.Structure):
    _fields_ = [("x", ctypes.c_double), ("y", ctypes.c_double), ("w", ctypes.c_double)]


class NurbsSolveResult(ctypes.Structure):
    _fields_ = [
        ("state", ctypes.c_int), ("u", ctypes.c_double), ("v", ctypes.c_double),
        ("point1_x", ctypes.c_double), ("point1_y", ctypes.c_double),
        ("point2_x", ctypes.c_double), ("point2_y", ctypes.c_double),
        ("distance", ctypes.c_double), ("iterations", ctypes.c_int),
        ("converged", ctypes.c_int),
    ]


class NurbsRoot(ctypes.Structure):
    _fields_ = [
        ("u", ctypes.c_double), ("v", ctypes.c_double), ("distance", ctypes.c_double),
        ("point1_x", ctypes.c_double), ("point1_y", ctypes.c_double),
        ("point2_x", ctypes.c_double), ("point2_y", ctypes.c_double),
        ("iterations", ctypes.c_int),
    ]


def load(path=None):
    """加载 DLL 并配置函数签名；返回 NurbsLib（命名空间对象）。"""
    if path is None:
        here = os.path.dirname(os.path.abspath(__file__))
        cand = os.path.join(os.path.dirname(here), "build", _NURBS)
        if os.path.exists(cand):
            path = cand
    lib = ctypes.CDLL(path)

    lib.nurbs_curve_create.argtypes = [
        ctypes.POINTER(NurbsControlPoint), ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.nurbs_curve_create.restype = ctypes.c_int
    lib.nurbs_curve_destroy.argtypes = [ctypes.c_void_p]
    lib.nurbs_curve_destroy.restype = None

    for fn in ("nurbs_evaluate_point", "nurbs_evaluate_first_derivative",
               "nurbs_evaluate_second_derivative"):
        getattr(lib, fn).argtypes = [ctypes.c_void_p, ctypes.c_double,
                                     ctypes.POINTER(ctypes.c_double),
                                     ctypes.POINTER(ctypes.c_double)]
        getattr(lib, fn).restype = ctypes.c_int

    lib.nurbs_unified_solve.argtypes = [
        ctypes.c_void_p, ctypes.c_void_p,
        ctypes.POINTER(ctypes.POINTER(NurbsSolveResult)), ctypes.POINTER(ctypes.c_int),
    ]
    lib.nurbs_unified_solve.restype = ctypes.c_int
    lib.nurbs_solve_results_free.argtypes = [ctypes.POINTER(NurbsSolveResult)]
    lib.nurbs_solve_results_free.restype = None

    lib.nurbs_find_all_intersections.argtypes = [
        ctypes.c_void_p, ctypes.c_void_p,
        ctypes.POINTER(ctypes.POINTER(NurbsRoot)), ctypes.POINTER(ctypes.c_int),
    ]
    lib.nurbs_find_all_intersections.restype = ctypes.c_int
    lib.nurbs_roots_free.argtypes = [ctypes.POINTER(NurbsRoot)]
    lib.nurbs_roots_free.restype = None
    return lib


def _clamped_knots(n, p):
    k = [0.0] * (p + 1)
    mid = n - p - 1
    for i in range(1, mid + 1):
        k.append(i / (mid + 1))
    k += [1.0] * (p + 1)
    return k


class Curve:
    """不透明曲线句柄的高层封装。"""

    def __init__(self, lib, handle):
        self.lib = lib
        self.handle = handle
        self._closed = False

    @classmethod
    def from_points(cls, points, degree, weights=None):
        """points: [(x,y), ...]；weights 可选（默认 1）。构建 clamped NURBS 曲线。"""
        n = len(points)
        if weights is None:
            weights = [1.0] * n
        arr = (NurbsControlPoint * n)()
        for i, (pt, w) in enumerate(zip(points, weights)):
            arr[i].x, arr[i].y, arr[i].w = float(pt[0]), float(pt[1]), float(w)
        knots = _clamped_knots(n, degree)
        karr = (ctypes.c_double * len(knots))(*knots)
        lib = cls._default_lib
        h = ctypes.c_void_p()
        rc = lib.nurbs_curve_create(arr, n, karr, len(knots), degree, ctypes.byref(h))
        if rc != OK:
            raise ValueError(f"曲线创建失败，错误码={rc}（非法输入/权重/节点）")
        return cls(lib, h)

    def _ensure_alive(self):
        if self._closed:
            raise RuntimeError("曲线句柄已释放")

    def evaluate(self, u):
        self._ensure_alive()
        x, y = ctypes.c_double(), ctypes.c_double()
        rc = self.lib.nurbs_evaluate_point(self.handle, u, ctypes.byref(x), ctypes.byref(y))
        if rc != OK:
            raise ValueError(f"求值失败 rc={rc}")
        return x.value, y.value

    def derivative(self, u):
        self._ensure_alive()
        x, y = ctypes.c_double(), ctypes.c_double()
        rc = self.lib.nurbs_evaluate_first_derivative(self.handle, u, ctypes.byref(x), ctypes.byref(y))
        if rc != OK:
            raise ValueError(f"一阶导失败 rc={rc}")
        return x.value, y.value

    def unified_solve(self, other):
        self._ensure_alive()
        res = ctypes.POINTER(NurbsSolveResult)()
        n = ctypes.c_int()
        rc = self.lib.nurbs_unified_solve(self.handle, other.handle,
                                          ctypes.byref(res), ctypes.byref(n))
        if rc != OK:
            raise ValueError(f"求解失败 rc={rc}")
        out = [res[i] for i in range(n.value)]
        self.lib.nurbs_solve_results_free(res)
        return out

    def find_intersections(self, other):
        self._ensure_alive()
        roots = ctypes.POINTER(NurbsRoot)()
        n = ctypes.c_int()
        rc = self.lib.nurbs_find_all_intersections(self.handle, other.handle,
                                                    ctypes.byref(roots), ctypes.byref(n))
        if rc != OK:
            raise ValueError(f"全交点失败 rc={rc}")
        out = [roots[i] for i in range(n.value)]
        self.lib.nurbs_roots_free(roots)
        return out

    def destroy(self):
        if not self._closed:
            self.lib.nurbs_curve_destroy(self.handle)
            self._closed = True

    def __del__(self):
        self.destroy()


def best_solve(results):
    """unified_solve 结果里取距离最小者。"""
    if not results:
        return None
    return min(results, key=lambda r: r.distance)


# 默认库（懒加载）：供 Curve.from_points 使用
_Curve_default_lib_holder = {}


def set_default_lib(lib):
    Curve._default_lib = lib


Curve._default_lib = None
