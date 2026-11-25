#!/usr/bin/env python3
"""Face-level OCCT section comparison for a directory of STEP files.

For every pair of non-adjacent B-rep faces, build a BRepAlgoAPI_Section and count
intersecting pairs. This is a face-level check; patch pairs are not one-to-one with
face pairs.
"""
import argparse
from pathlib import Path

from OCP.STEPControl import STEPControl_Reader
from OCP.IFSelect import IFSelect_RetDone
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE
from OCP.TopoDS import TopoDS
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section


def process(path: Path):
    r = STEPControl_Reader()
    if r.ReadFile(str(path)) != IFSelect_RetDone:
        return None, 0, 0
    r.TransferRoots()
    shape = r.OneShape()
    faces = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        faces.append(TopoDS.Face(ex.Current()))
        ex.Next()
    edges = []
    for f in faces:
        es = []
        ee = TopExp_Explorer(f, TopAbs_EDGE)
        while ee.More():
            es.append(TopoDS.Edge(ee.Current()))
            ee.Next()
        edges.append(es)

    def adjacent(i, j):
        return any(a.IsSame(b) for a in edges[i] for b in edges[j])

    pairs = inter = 0
    for i in range(len(faces)):
        for j in range(i + 1, len(faces)):
            if adjacent(i, j):
                continue
            pairs += 1
            try:
                sec = BRepAlgoAPI_Section(faces[i], faces[j], False)
                sec.Build()
                if not sec.IsDone():
                    continue
                ne = 0
                ee = TopExp_Explorer(sec.Shape(), TopAbs_EDGE)
                while ee.More():
                    ne += 1
                    ee.Next()
                if ne > 0:
                    inter += 1
            except Exception:
                pass
    return len(faces), pairs, inter


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, default=Path("data/raw/fusion"))
    ap.add_argument("--pattern", default="*.stp")
    ap.add_argument("--out", type=Path, default=Path("results/occt_pairs.tsv"))
    args = ap.parse_args()
    files = sorted(args.input_dir.glob(args.pattern))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    rows = ["file\tfaces\tnonadj_pairs\tocct_intersecting"]
    total = [0, 0, 0]
    for f in files:
        nf, pairs, inter = process(f)
        if nf is None:
            continue
        rows.append(f"{f.name}\t{nf}\t{pairs}\t{inter}")
        total[0] += nf
        total[1] += pairs
        total[2] += inter
        print(f"{f.name}: faces={nf} pairs={pairs} occt_inter={inter}")
    rows.append(f"TOTAL\t{total[0]}\t{total[1]}\t{total[2]}")
    args.out.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print("TOTAL", total, "->", args.out)


if __name__ == "__main__":
    main()
