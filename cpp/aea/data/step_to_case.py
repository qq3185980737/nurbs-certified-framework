# -*- coding: utf-8 -*-
"""Convert an ABC STEP model into bench_cert_step data format (SURF patches)."""
import sys
from pathlib import Path
from OCP.STEPControl import STEPControl_Reader
from OCP.IFSelect import IFSelect_RetDone
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE
from OCP.TopoDS import TopoDS
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCP.GeomConvert import GeomConvert_BSplineSurfaceToBezierSurface

step = Path(sys.argv[1])
out = Path(sys.argv[2])
reader = STEPControl_Reader()
if reader.ReadFile(str(step)) != IFSelect_RetDone:
    raise SystemExit('STEP read failed')
reader.TransferRoots()
shape = reader.OneShape()
nshape = BRepBuilderAPI_NurbsConvert(shape).Shape()

faces = []
ex = TopExp_Explorer(nshape, TopAbs_FACE)
while ex.More():
    faces.append(TopoDS.Face(ex.Current()))
    ex.Next()

# B-rep topological adjacency: faces sharing an edge
face_edges = []
ex = TopExp_Explorer(nshape, TopAbs_EDGE)
# collect edges per face by exploring each face
from OCP.TopAbs import TopAbs_EDGE
for face in faces:
    es = []
    ee = TopExp_Explorer(face, TopAbs_EDGE)
    while ee.More():
        es.append(TopoDS.Edge(ee.Current()))
        ee.Next()
    face_edges.append(es)

adjacent = []
for i in range(len(faces)):
    for j in range(i + 1, len(faces)):
        if any(a.IsSame(b) for a in face_edges[i] for b in face_edges[j]):
            adjacent.append((i, j))

lines = ['NSURF %d' % len(faces)]
for fi, face in enumerate(faces):
    ad = BRepAdaptor_Surface(face)
    try:
        bs = ad.BSpline()
        conv = GeomConvert_BSplineSurfaceToBezierSurface(bs)
    except Exception as e:
        print('face %d skip: %s' % (fi, e), file=sys.stderr)
        continue
    for i in range(1, conv.NbUPatches() + 1):
        for j in range(1, conv.NbVPatches() + 1):
            bz = conv.Patch(i, j)
            rat = bz.IsURational() or bz.IsVRational()
            lines.append('SURF %d %d %d %d %d %d' % (
                fi, (i - 1) * conv.NbVPatches() + (j - 1),
                bz.UDegree(), bz.VDegree(), bz.NbUPoles(), bz.NbVPoles()))
            for a in range(1, bz.NbUPoles() + 1):
                for b in range(1, bz.NbVPoles() + 1):
                    p = bz.Pole(a, b)
                    w = bz.Weight(a, b) if rat else 1.0
                    lines.append('%.17g %.17g %.17g %.17g' % (p.X(), p.Y(), p.Z(), w))
            lines.append('ENDSURF')
lines.append('NADJ %d' % len(adjacent))
for i, j in adjacent:
    lines.append('%d %d' % (i, j))
out.write_text('\n'.join(lines) + '\n', encoding='utf-8')
print('faces=%d lines=%d -> %s' % (len(faces), len(lines), out))
