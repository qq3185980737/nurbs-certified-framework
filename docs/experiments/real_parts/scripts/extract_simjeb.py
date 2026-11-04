#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""extract_simjeb.py — SimJEB STEP -> 有理 Bezier patch 文本（只导 SURF，不导 CURVE）。

重建依据：docs/experiments/step_case/extract_bracket.py 的
    STEPControl_Reader -> [可选 BRepBuilderAPI_NurbsConvert] -> 逐 face
    BRepAdaptor_Surface.BSpline() -> GeomConvert_BSplineSurfaceToBezierSurface
思路；SimJEB 提取只保留曲面（本基准只需要 patch-pair SSI 认证）。

输出文本格式（与 docs/experiments/real_parts/data/*.txt 完全一致）：
    NSURF <全部 patch 总数>
    SURF <fid> <pi> <du> <dv> <nu> <nv>
    nu*nv 行 "x y z w"   （i 为外层、j 为内层；x,y,z 为笛卡尔控制点，w 为权）
    ENDSURF

face 顺序 = TopExp_Explorer(shape, TopAbs_FACE) 顺序，fid 从 0 连续编号；
pi = (k-1)*NbVPatches + (l-1)，k/l 为 GeomConvert 的 U/V patch 序号。

用法（WSL）：
  /root/miniconda3/envs/occ/bin/python extract_simjeb.py \
      --stp ../../stp/228.stp --out /tmp/228.txt [--meta /tmp/228_meta.csv] \
      [--no-nurbs-convert] [--max-faces N]
"""
import argparse
import sys
import time

from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.TopoDS import topods
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.GeomConvert import GeomConvert_BSplineSurfaceToBezierSurface

try:
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
    _HAS_NURBS_CONVERT = True
except Exception:
    _HAS_NURBS_CONVERT = False

# 可选 metadata：曲面类型 / 参数域 / 面积
_SURF_TYPE_NAMES = None
try:
    from OCC.Core.GeomAbs import (GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone,
                                  GeomAbs_Sphere, GeomAbs_Torus, GeomAbs_BezierSurface,
                                  GeomAbs_BSplineSurface, GeomAbs_SurfaceOfRevolution,
                                  GeomAbs_SurfaceOfExtrusion, GeomAbs_OffsetSurface,
                                  GeomAbs_OtherSurface)
    _SURF_TYPE_NAMES = {
        GeomAbs_Plane: "Plane", GeomAbs_Cylinder: "Cylinder", GeomAbs_Cone: "Cone",
        GeomAbs_Sphere: "Sphere", GeomAbs_Torus: "Torus",
        GeomAbs_BezierSurface: "BezierSurface", GeomAbs_BSplineSurface: "BSplineSurface",
        GeomAbs_SurfaceOfRevolution: "SurfaceOfRevolution",
        GeomAbs_SurfaceOfExtrusion: "SurfaceOfExtrusion",
        GeomAbs_OffsetSurface: "OffsetSurface", GeomAbs_OtherSurface: "OtherSurface",
    }
except Exception:
    pass

try:
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRepGProp import brepgprop
    _HAS_GPROP = True
except Exception:
    _HAS_GPROP = False


def surf_type_name(ad):
    if _SURF_TYPE_NAMES is None:
        return "?"
    try:
        return _SURF_TYPE_NAMES.get(ad.GetType(), "Unknown")
    except Exception:
        return "?"


def face_area(face):
    if not _HAS_GPROP:
        return float("nan")
    try:
        props = GProp_GProps()
        brepgprop.SurfaceProperties(face, props)
        return props.Mass()
    except Exception:
        return float("nan")


def read_step(path):
    r = STEPControl_Reader()
    status = r.ReadFile(path)
    if status != IFSelect_RetDone:
        raise RuntimeError("STEPControl_Reader.ReadFile failed on %s (status=%s)" % (path, status))
    r.TransferRoots()
    sh = r.OneShape()
    if sh.IsNull():
        raise RuntimeError("STEP read produced null shape: %s" % path)
    return sh


def surf_rows(bz):
    """返回 (du, dv, nu, nv, rows) ；rows 按 i 外层 j 内层，元素 (x,y,z,w)。"""
    rat = bz.IsURational() or bz.IsVRational()
    nu, nv = bz.NbUPoles(), bz.NbVPoles()
    rows = []
    for i in range(1, nu + 1):
        for j in range(1, nv + 1):
            p = bz.Pole(i, j)
            w = bz.Weight(i, j) if rat else 1.0
            rows.append((p.X(), p.Y(), p.Z(), w))
    return bz.UDegree(), bz.VDegree(), nu, nv, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--meta", default=None)
    ap.add_argument("--no-nurbs-convert", action="store_true",
                    help="不先做 BRepBuilderAPI_NurbsConvert（直接对原 face 取 BSpline）")
    ap.add_argument("--max-faces", type=int, default=0)
    args = ap.parse_args()

    t0 = time.time()
    sh = read_step(args.stp)
    t_read = time.time() - t0
    used_convert = False
    if not args.no_nurbs_convert and _HAS_NURBS_CONVERT:
        try:
            sh = BRepBuilderAPI_NurbsConvert(sh).Shape()
            used_convert = True
        except Exception as e:
            print("WARN NurbsConvert failed: %r" % (e,), flush=True)

    exp = TopExp_Explorer(sh, TopAbs_FACE)
    faces = []
    while exp.More():
        faces.append(topods.Face(exp.Current()))
        exp.Next()
        if args.max_faces and len(faces) >= args.max_faces:
            break

    meta_rows = []
    patch_blocks = []   # (fid, pi, du, dv, nu, nv, rows)
    n_conv_fail = 0
    t_conv0 = time.time()
    for fid, face in enumerate(faces):
        ad = BRepAdaptor_Surface(face)
        try:
            cv = GeomConvert_BSplineSurfaceToBezierSurface(ad.BSpline())
            nu_k, nv_l = cv.NbUPatches(), cv.NbVPatches()
        except Exception:
            n_conv_fail += 1
            meta_rows.append((fid, surf_type_name(ad), 0, 0, 0, 0, 0.0,
                              ad.FirstUParameter(), ad.LastUParameter(),
                              ad.FirstVParameter(), ad.LastVParameter(),
                              face_area(face)))
            continue
        pi = 0
        for k in range(1, nu_k + 1):
            for l in range(1, nv_l + 1):
                try:
                    bz = cv.Patch(k, l)
                except Exception:
                    continue
                du, dv, nu, nv, rows = surf_rows(bz)
                patch_blocks.append((fid, pi, du, dv, nu, nv, rows))
                pi += 1
        meta_rows.append((fid, surf_type_name(ad), nu_k, nv_l, pi,
                          sum(b[4] * b[5] for b in patch_blocks[-pi:]) if pi else 0,
                          0.0,
                          ad.FirstUParameter(), ad.LastUParameter(),
                          ad.FirstVParameter(), ad.LastVParameter(),
                          face_area(face)))
    t_conv = time.time() - t_conv0

    n_patch = len(patch_blocks)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("NSURF %d\n" % n_patch)
        for (fid, pi, du, dv, nu, nv, rows) in patch_blocks:
            fh.write("SURF %d %d %d %d %d %d\n" % (fid, pi, du, dv, nu, nv))
            for (x, y, z, w) in rows:
                fh.write("%.17g %.17g %.17g %.17g\n" % (x, y, z, w))
            fh.write("ENDSURF\n")

    if args.meta:
        with open(args.meta, "w", encoding="utf-8") as fh:
            fh.write("fid,type,nupat,nvpat,npoles,area,u0,u1,v0,v1\n")
            for (fid, typ, nup, nvp, npole, _dummy, _a, u0, u1, v0, v1, area) in meta_rows:
                fh.write("%d,%s,%d,%d,%d,%.12g,%.12g,%.12g,%.12g,%.12g\n"
                         % (fid, typ, nup, nvp, npole, area, u0, u1, v0, v1))

    n_faces = len(faces)
    print("EXTRACT stp=%s faces=%d patches=%d conv_fail=%d nurbs_convert=%s "
          "read=%.2fs convert=%.2fs total=%.2fs"
          % (args.stp, n_faces, n_patch, n_conv_fail, used_convert,
             t_read, t_conv, time.time() - t0), flush=True)
    print("WROTE %s%s" % (args.out, (" meta=" + args.meta) if args.meta else ""), flush=True)


if __name__ == "__main__":
    sys.exit(main())
