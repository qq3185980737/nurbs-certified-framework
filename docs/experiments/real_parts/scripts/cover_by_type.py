#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""cover_by_type.py — 把覆盖检验的"点落在 cover 盒像内"比例按 pair 类型拆分。

背景：cover_check 对"OCCT 判定有交线"的抽样对直接调用 cover()（绕过 classify_pair
的 type 判定）。type=4(共面) 的对在 classify_pair 里走提前返回，cover() 对其会触顶；
type=3(undecided) 的对本身就没有完整 certified 覆盖。因此必须按 type 拆开看：
  - type=1（1D certified）：cover 必须覆盖真交点 —— 这才是覆盖检验的检验对象。
  - type=3/4：cover 触顶/未认证，miss 属预期，需单独报告。
输出 results/cover_by_type.csv
"""
import csv
import glob
import os

RP = "/mnt/d/qqdownload/nurbs-certified-framework/docs/experiments/real_parts"
RES = os.path.join(RP, "results")


def main():
    rows = []
    for cf in sorted(glob.glob(os.path.join(RES, "cover", "*_cover.csv"))):
        part = os.path.basename(cf).replace("_cover.csv", "")
        pairs = {}
        pf = os.path.join(RES, "probe", part + "_pairs.csv")
        if os.path.exists(pf):
            for r in csv.DictReader(open(pf, encoding="utf-8")):
                pairs[(r["fidA"], r["piA"], r["fidB"], r["piB"])] = r
        agg = {}
        for r in csv.DictReader(open(cf, encoding="utf-8")):
            key = (r["fidA"], r["piA"], r["fidB"], r["piB"])
            p = pairs.get(key)
            t = int(p["type"]) if p else -1
            a = agg.setdefault(t, {"pairs": 0, "pts": 0, "both": 0, "pairs_all": 0, "cover_boxes": 0})
            a["pairs"] += 1
            a["pts"] += int(r["npoints"])
            a["both"] += int(r["inside_both"])
            a["cover_boxes"] += int(r["ncover"])
            if int(r["inside_both"]) == int(r["npoints"]):
                a["pairs_all"] += 1
        for t, a in sorted(agg.items()):
            rows.append({"part": part, "type": t, "pairs": a["pairs"], "points": a["pts"],
                         "points_inside_both": a["both"], "pairs_all_inside": a["pairs_all"],
                         "point_cover_rate": round(a["both"] / max(1, a["pts"]), 4),
                         "avg_cover_boxes": round(a["cover_boxes"] / max(1, a["pairs"]), 1)})
    out = os.path.join(RES, "cover_by_type.csv")
    cols = ["part", "type", "pairs", "points", "points_inside_both", "pairs_all_inside",
            "point_cover_rate", "avg_cover_boxes"]
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    # 全局按 type 汇总
    tot = {}
    for r in rows:
        a = tot.setdefault(r["type"], {"pairs": 0, "points": 0, "both": 0, "pairs_all": 0})
        a["pairs"] += r["pairs"]; a["points"] += r["points"]
        a["both"] += r["points_inside_both"]; a["pairs_all"] += r["pairs_all_inside"]
    print("type  pairs  points  inside_both  rate    pairs_all_inside")
    for t, a in sorted(tot.items()):
        print("%4d  %5d  %6d  %7d  %6.2f%%  %d/%d"
              % (t, a["pairs"], a["points"], a["both"],
                 100.0 * a["both"] / max(1, a["points"]), a["pairs_all"], a["pairs"]))
    print("WROTE", out)


if __name__ == "__main__":
    main()
