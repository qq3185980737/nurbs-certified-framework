# -*- coding: utf-8 -*-
"""P0-2：装配干涉场景（支架 + 横杆）→ Bezier 数据 + OCCT 边-面交点，导出纯文本。"""
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCC.Core.BRepCheck import BRepCheck_Analyzer
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE
from OCC.Core.TopoDS import topods
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface, BRepAdaptor_Curve
from OCC.Core.GeomConvert import GeomConvert_BSplineSurfaceToBezierSurface, GeomConvert_BSplineCurveToBezierCurve
from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.gp import gp_Pnt, gp_Ax2, gp_Dir
from OCC.Core.BRep import BRep_Tool

OUT = "/mnt/c/Users/15632/Documents/deepseek-harness/default-workspace/cad_review12/step_case"

# ---- 零件 A：底板 + 凸台 + 圆角 + 两孔；零件 B：贯穿凸台的横杆 ----
plate = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(80, 60, 10)).Shape()
boss = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(40, 30, 10), gp_Dir(0, 0, 1)), 12.0, 18.0).Shape()
A = BRepAlgoAPI_Fuse(plate, boss).Shape()
mk = BRepFilletAPI_MakeFillet(A)
exp = TopExp_Explorer(A, TopAbs_EDGE)
c = 0
while exp.More() and c < 6:
    try:
        mk.Add(2.0, topods.Edge(exp.Current())); c += 1
    except Exception:
        pass
    exp.Next()
try:
    A = mk.Shape()
except Exception:
    pass
for (cx, cy) in ((15, 15), (65, 45)):
    A = BRepAlgoAPI_Cut(A, BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(cx, cy, -5), gp_Dir(0, 0, 1)), 4.0, 20.0).Shape()).Shape()
B = BRepPrimAPI_MakeBox(gp_Pnt(-5, 26, 14), gp_Pnt(90, 34, 24)).Shape()
print("A valid:", BRepCheck_Analyzer(A).IsValid(), "| B valid:", BRepCheck_Analyzer(B).IsValid(), flush=True)

NA = BRepBuilderAPI_NurbsConvert(A).Shape()
NB = BRepBuilderAPI_NurbsConvert(B).Shape()

shapes = [("A", NA), ("B", NB)]
faces, edges = [], []
for tag, sh in shapes:
    exp = TopExp_Explorer(sh, TopAbs_FACE)
    while exp.More():
        faces.append((tag, topods.Face(exp.Current()))); exp.Next()
    exp = TopExp_Explorer(sh, TopAbs_EDGE)
    while exp.More():
        edges.append((tag, topods.Edge(exp.Current()))); exp.Next()
print("faces:", len(faces), "edges:", len(edges), flush=True)


def conv_surf(face):
    ad = BRepAdaptor_Surface(face)
    try:
        cv = GeomConvert_BSplineSurfaceToBezierSurface(ad.BSpline())
        nu, nv = cv.NbUPatches(), cv.NbVPatches()
    except Exception:
        return [], (0, 1, 0, 1)
    out = []
    for k in range(1, nu + 1):
        for l in range(1, nv + 1):
            try:
                out.append(cv.Patch(k, l))
            except Exception:
                pass
    return out, (ad.FirstUParameter(), ad.LastUParameter(), ad.FirstVParameter(), ad.LastVParameter())


def conv_curve(edge):
    ad = BRepAdaptor_Curve(edge)
    try:
        cv = GeomConvert_BSplineCurveToBezierCurve(ad.BSpline())
        n = cv.NbArcs()
    except Exception:
        return [], (0, 1)
    return [cv.Arc(k) for k in range(1, n + 1)], (ad.FirstParameter(), ad.LastParameter())


def sp_rows(bz):
    rat = bz.IsURational() or bz.IsVRational()
    rows = []
    for i in range(1, bz.NbUPoles() + 1):
        for j in range(1, bz.NbVPoles() + 1):
            p = bz.Pole(i, j); w = bz.Weight(i, j) if rat else 1.0
            rows.append((p.X(), p.Y(), p.Z(), w))
    return bz.UDegree(), bz.VDegree(), bz.NbUPoles(), bz.NbVPoles(), rows


def cu_rows(bz):
    try:
        rat = bz.IsRational()
    except Exception:
        rat = False
    rows = []
    for i in range(1, bz.NbPoles() + 1):
        p = bz.Pole(i); w = bz.Weight(i) if rat else 1.0
        rows.append((p.X(), p.Y(), p.Z(), w))
    return bz.Degree(), bz.NbPoles(), rows


fdat, edat = [], []
for fi, (tag, f) in enumerate(faces):
    ad = BRepAdaptor_Surface(f)
    patches, rng = conv_surf(f)
    fdat.append({"id": fi, "tag": tag, "ad": ad, "surf": ad.Surface().Surface(),
                 "patches": patches, "range": rng})
for ei, (tag, e) in enumerate(edges):
    ad = BRepAdaptor_Curve(e)
    segs, rng = conv_curve(e)
    edat.append({"id": ei, "tag": tag, "ad": ad, "curve": ad.Curve(),
                 "segs": segs, "range": rng})


def bb_curve(curve, rng):
    xs, ys, zs = [], [], []
    for i in range(21):
        t = rng[0] + (rng[1] - rng[0]) * i / 20.0
        p = curve.Value(t); xs.append(p.X()); ys.append(p.Y()); zs.append(p.Z())
    return min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)


def bb_surf(surf, rng):
    xs, ys, zs = [], [], []
    for i in range(9):
        for j in range(9):
            u = rng[0] + (rng[1] - rng[0]) * i / 8.0
            v = rng[2] + (rng[3] - rng[2]) * j / 8.0
            p = surf.Value(u, v); xs.append(p.X()); ys.append(p.Y()); zs.append(p.Z())
    return min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)


pairs = []
for e in edat:
    if e["curve"] is None:
        continue
    be = bb_curve(e["curve"], e["range"])
    for f in fdat:
        if e["tag"] == f["tag"]:
            continue          # 只取跨零件对（装配干涉）
        bf = bb_surf(f["surf"], f["range"])
        if not (be[1] >= bf[0] - 1e-9 and bf[1] >= be[0] - 1e-9 and
                be[3] >= bf[2] - 1e-9 and bf[3] >= be[2] - 1e-9 and
                be[5] >= bf[4] - 1e-9 and bf[5] >= be[4] - 1e-9):
            continue
        try:
            res = GeomAPI_IntCS(e["curve"], f["surf"])
        except Exception:
            continue
        if not res.IsDone():
            continue
        pts = []
        for k in range(1, res.NbPoints() + 1):
            try:
                p = res.Point(k)
                seq = list(res.Parameters(k))
                uu, vv, ww = seq[0], seq[1], (seq[2] if len(seq) > 2 else 0.0)
            except Exception:
                continue
            r, fr = e["range"], f["range"]
            if not (r[0] - 1e-7 <= uu <= r[1] + 1e-7):
                continue
            if not (fr[0] - 1e-7 <= vv <= fr[1] + 1e-7 and fr[2] - 1e-7 <= ww <= fr[3] + 1e-7):
                continue
            pts.append((p.X(), p.Y(), p.Z()))
        if pts:
            pairs.append({"edge": e["id"], "face": f["id"], "pts": pts})

print("跨零件边-面对（OCCT 有交点）:", len(pairs), "| 总交点:", sum(len(p["pts"]) for p in pairs), flush=True)

with open(OUT + "/assembly_data.txt", "w", encoding="utf-8") as fh:
    fh.write("NSURF %d\n" % len(fdat))
    for f in fdat:
        for pi, bz in enumerate(f["patches"]):
            du, dv, nu, nv, rows = sp_rows(bz)
            fh.write("SURF %d %d %d %d %d %d\n" % (f["id"], pi, du, dv, nu, nv))
            for (x, y, z, w) in rows:
                fh.write("%.17g %.17g %.17g %.17g\n" % (x, y, z, w))
            fh.write("ENDSURF\n")
    fh.write("NCURVE %d\n" % len(edat))
    for e in edat:
        for si, bz in enumerate(e["segs"]):
            deg, n, rows = cu_rows(bz)
            fh.write("CURVE %d %d %d %d\n" % (e["id"], si, deg, n))
            for (x, y, z, w) in rows:
                fh.write("%.17g %.17g %.17g %.17g\n" % (x, y, z, w))
            fh.write("ENDCURVE\n")
    fh.write("NPAIR %d\n" % len(pairs))
    for p in pairs:
        fh.write("PAIR %d %d %d\n" % (p["edge"], p["face"], len(p["pts"])))
        for (x, y, z) in p["pts"]:
            fh.write("%.17g %.17g %.17g\n" % (x, y, z))
        fh.write("ENDPAIR\n")
print("已导出", OUT + "/assembly_data.txt", flush=True)
