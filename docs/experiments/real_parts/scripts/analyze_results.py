#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""analyze_results.py — 汇总真实零件基准的所有 CSV，生成可写入论文的表格。

输入（全部由本目录脚本产生）：
  results/batch/{A,H,B,C}/results.csv       有界批量运行
  results/probe/<part>_probe.log            probe_pairs 的 PROBE 行
  results/probe/<part>_pairs.csv            逐对 classify_pair 结果
  results/probe/<part>_occt.csv             GeomAPI_IntSS 结果
  results/cover/<part>_cover.csv            cover_check 覆盖检验
  logs/<part>.log                           历史运行日志
输出：
  results/batch/all_results.csv
  results/repro_dualbinary.csv
  results/probe_summary.csv
  results/occt_confusion.csv
  results/occt_conflicts.csv
  results/cover_summary.csv
  results/stats.json
"""
import csv
import glob
import json
import os
import re
import sys

ROOT = "/mnt/d/qqdownload/nurbs-certified-framework"
RP = os.path.join(ROOT, "docs/experiments/real_parts")
RES = os.path.join(RP, "results")
HEAD_COMMIT = "f60a193e8d1fb5cf8b566e19b7729b00e7c1f699"

RESULT_RE = re.compile(
    r"RESULT faces=(\d+) facepairs=(\d+) \| patchpairs=(\d+) bboxskip=(\d+) try=(\d+) \| "
    r"disjoint=(\d+) \| 1D=(\d+) \(comps=(\d+)\) \| 2D=(\d+) \| coplanar=(\d+) \| undecided=(\d+)")
BOXES_RE = re.compile(r"certified_boxes=(\d+) undecided_boxes=(\d+) time=([\d.]+) ms")
PARSED_RE = re.compile(r"parsed (\d+) faces \((\d+) patches\), tau=([\d.]+)")
PAIRSTAT_RE = re.compile(
    r"PAIRSTAT tau=[\d.]+ cap=(\d+) n=(\d+) med_ms=([\d.]+) p95_ms=([\d.]+) p99_ms=([\d.]+) max_ms=([\d.]+) "
    r"cap_hit=(\d+) cover_gt2500=(\d+) pair_t0=(\d+) pair_t1=(\d+) pair_t2=(\d+) pair_t3=(\d+) pair_t4=(\d+) "
    r"und_pair_ratio=([\d.]+) t0_ms=([\d.]+) t1_ms=([\d.]+) t2_ms=([\d.]+) t3_ms=([\d.]+) t4_ms=([\d.]+)")
PROBE_RE = re.compile(
    r"PROBE data=(\S+) faces=(\d+) patches=(\d+) facepairs=(\d+) overlap_facepairs=(\d+) "
    r"patchpairs=(\d+) overlap_patchpairs=(\d+) cand=(\d+) pos=(\d+) neg=(\d+) tau=([\d.]+) cap=(\d+)")


def parse_log_text(text):
    d = {}
    m = PARSED_RE.search(text)
    if m:
        d["faces"] = int(m.group(1)); d["patches"] = int(m.group(2))
    m = RESULT_RE.search(text)
    if m:
        keys = ["faces", "facepairs", "patchpairs", "bboxskip", "try", "disjoint",
                "d1", "comps", "d2", "coplanar", "undecided"]
        for k, v in zip(keys, m.groups()):
            d[k] = int(v)
    m = BOXES_RE.search(text)
    if m:
        d["certified_boxes"] = int(m.group(1)); d["undecided_boxes"] = int(m.group(2))
        d["method_ms"] = float(m.group(3))
    m = PAIRSTAT_RE.search(text)
    if m:
        keys = ["node_cap", "pair_n", "pair_med_ms", "pair_p95_ms", "pair_p99_ms",
                "pair_max_ms", "cap_hit", "cover_gt2500", "pair_t0", "pair_t1",
                "pair_t2", "pair_t3", "pair_t4", "und_pair_ratio", "t0_ms", "t1_ms",
                "t2_ms", "t3_ms", "t4_ms"]
        for k, v in zip(keys, m.groups()):
            d[k] = float(v) if "." in v else int(v)
    # 定时输出会有多条 PAIRSTAT：取最后一条（超时运行 = 最近累计统计）
    allm = list(PAIRSTAT_RE.finditer(text))
    if len(allm) > 1:
        m = allm[-1]
        keys = ["node_cap", "pair_n", "pair_med_ms", "pair_p95_ms", "pair_p99_ms",
                "pair_max_ms", "cap_hit", "cover_gt2500", "pair_t0", "pair_t1",
                "pair_t2", "pair_t3", "pair_t4", "und_pair_ratio", "t0_ms", "t1_ms",
                "t2_ms", "t3_ms", "t4_ms"]
        for k, v in zip(keys, m.groups()):
            d[k] = float(v) if "." in v else int(v)
    return d


def read_batch(group):
    path = os.path.join(RES, "batch", group, "results.csv")
    if not os.path.exists(path):
        return []
    return list(csv.DictReader(open(path, encoding="utf-8")))


def fnum(x, nd=2):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return round(v, nd)


# ---------------------------------------------------------------- 1. merge batch
def merge_batch():
    allrows = []
    for g in ("A", "H", "B", "C"):
        allrows.extend(read_batch(g))
    out = os.path.join(RES, "batch", "all_results.csv")
    if allrows:
        keys = list(allrows[0].keys())
        with open(out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(allrows)
    return allrows


# ------------------------------------------------- 2. HEAD vs WORK vs HIST 双跑
def repro_dual():
    versions = {"HEAD": ("H", "git HEAD " + HEAD_COMMIT),
                "WORK": ("A", "git HEAD + 未提交修改 (+25/-3)")}
    rows = []
    for part in ("228", "387", "107", "411"):
        hist = parse_log_text(open(os.path.join(RP, "logs", part + ".log"), encoding="utf-8",
                                   errors="replace").read()) if os.path.exists(
            os.path.join(RP, "logs", part + ".log")) else {}
        if hist:
            r = {"part": part, "version": "HIST", "commit": "历史 logs/ 记录（基线二进制）", "rc": "", "timeout": ""}
            r.update({k: hist.get(k, "") for k in
                      ["faces", "facepairs", "patchpairs", "bboxskip", "try", "disjoint", "d1",
                       "d2", "coplanar", "undecided", "certified_boxes", "undecided_boxes"]})
            r["method_s"] = round(hist.get("method_ms", 0) / 1000.0, 2) if hist.get("method_ms") else ""
            rows.append(r)
        for ver, (grp, commit) in versions.items():
            for r0 in read_batch(grp):
                if r0["part"] != part:
                    continue
                d = parse_log_text(open(r0["log"], encoding="utf-8", errors="replace").read()) \
                    if r0.get("log") and os.path.exists(r0["log"]) else {}
                r = {"part": part, "version": ver, "commit": commit,
                     "rc": r0.get("rc", ""), "timeout": r0.get("timeout", "")}
                r.update({k: d.get(k, "") for k in
                          ["faces", "facepairs", "patchpairs", "bboxskip", "try", "disjoint", "d1",
                           "d2", "coplanar", "undecided", "certified_boxes", "undecided_boxes"]})
                r["method_s"] = round(d["method_ms"] / 1000.0, 2) if d.get("method_ms") else ""
                rows.append(r)
    cols = ["part", "version", "commit", "rc", "timeout", "faces", "facepairs", "patchpairs",
            "bboxskip", "try", "disjoint", "d1", "d2", "coplanar", "undecided",
            "certified_boxes", "undecided_boxes", "method_s"]
    out = os.path.join(RES, "repro_dualbinary.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return rows


# --------------------------------------------------------- 3. probe + OCCT 汇总
def probe_analysis():
    parts = sorted(set(os.path.basename(p).replace("_pairs.csv", "")
                       for p in glob.glob(os.path.join(RES, "probe", "*_pairs.csv"))))
    rows = []
    occt_rows = []
    conflicts = []
    for part in parts:
        log = os.path.join(RES, "probe", part + "_probe.log")
        ptxt = open(log, encoding="utf-8", errors="replace").read() if os.path.exists(log) else ""
        m = PROBE_RE.search(ptxt)
        pinfo = {}
        if m:
            keys = ["data", "faces", "patches", "facepairs", "overlap_facepairs", "patchpairs",
                    "overlap_patchpairs", "cand", "pos", "neg", "tau", "cap"]
            pinfo = dict(zip(keys, m.groups()))
        pairs = list(csv.DictReader(open(os.path.join(RES, "probe", part + "_pairs.csv"), encoding="utf-8")))
        occ = {}
        opath = os.path.join(RES, "probe", part + "_occt.csv")
        if os.path.exists(opath):
            for r in csv.DictReader(open(opath, encoding="utf-8")):
                occ[(r["fidA"], r["piA"], r["fidB"], r["piB"])] = r
        pos = [r for r in pairs if r["src"] == "pos"]
        neg = [r for r in pairs if r["src"] == "neg"]

        def tm(rs):
            v = sorted(float(r["ms"]) for r in rs if r.get("ms"))
            if not v:
                return (0.0, 0.0, 0.0)
            n = len(v)
            med = v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])
            return (round(med, 2), round(sum(v) / n, 2), round(v[-1], 2))

        pos_med, pos_mean, pos_max = tm(pos)
        neg_med, neg_mean, neg_max = tm(neg)
        types_pos = {k: 0 for k in range(5)}
        for r in pos:
            types_pos[int(r["type"])] += 1

        def occ_stat(rs):
            has = no = err = 0
            for r in rs:
                o = occ.get((r["fidA"], r["piA"], r["fidB"], r["piB"]))
                if o is None:
                    continue
                ln = int(o["occt_lines"])
                if ln < 0:
                    err += 1
                elif ln >= 1:
                    has += 1
                else:
                    no += 1
            return has, no, err

        ph, pn, pe = occ_stat(pos)
        nh, nn, ne = occ_stat(neg)
        # 需求 3(b)：bbox 相交但认证 disjoint
        b = [r for r in pos if int(r["type"]) == 0]
        bh, bn, be = occ_stat(b)
        rows.append({
            "part": part, "faces": pinfo.get("faces", ""), "patches": pinfo.get("patches", ""),
            "facepairs": pinfo.get("facepairs", ""), "overlap_facepairs": pinfo.get("overlap_facepairs", ""),
            "patchpairs_total": pinfo.get("patchpairs", ""), "patchpairs_bbox_overlap": pinfo.get("overlap_patchpairs", ""),
            "n_pos": len(pos), "n_neg": len(neg),
            "pos_type0": types_pos[0], "pos_type1": types_pos[1], "pos_type2": types_pos[2],
            "pos_type3": types_pos[3], "pos_type4": types_pos[4],
            "pos_med_ms": pos_med, "pos_mean_ms": pos_mean, "pos_max_ms": pos_max,
            "neg_med_ms": neg_med, "neg_mean_ms": neg_mean, "neg_max_ms": neg_max,
            "occt_pos_has_line": ph, "occt_pos_no_line": pn, "occt_pos_err": pe,
            "occt_neg_has_line": nh, "occt_neg_no_line": nn, "occt_neg_err": ne,
            "b_bbox_overlap_cert_disjoint": len(b),
            "b_occt_no_line": bn, "b_occt_has_line": bh, "b_occt_err": be,
        })
        # 逐类混淆（用于 occt_confusion.csv）
        for cls, rs in (("pos_bbox_overlap", pos), ("neg_bbox_disjoint", neg)):
            for t in range(5):
                sub = [r for r in rs if int(r["type"]) == t]
                if not sub:
                    continue
                h, n0, e0 = occ_stat(sub)
                occt_rows.append({"part": part, "sample_class": cls, "method_type": t,
                                  "n_pairs": len(sub), "occt_has_line": h,
                                  "occt_no_line": n0, "occt_error": e0,
                                  "agree_rate": round((h if t in (1, 2, 4) else n0) / max(1, h + n0), 4)})
        # 冲突清单
        for r in pos:
            o = occ.get((r["fidA"], r["piA"], r["fidB"], r["piB"]))
            if not o:
                continue
            t = int(r["type"]); ln = int(o["occt_lines"])
            tag = ""
            if t == 0 and ln >= 1:
                tag = "SEVERE_method_disjoint_but_occt_line"
            elif t == 1 and ln == 0:
                tag = "certified_1D_but_occt_no_line"
            elif t == 3 and ln >= 1:
                tag = "undecided_but_occt_has_line"
            elif t == 4 and ln == 0:
                tag = "coplanar_but_occt_no_line"
            if tag:
                conflicts.append({"part": part, "fidA": r["fidA"], "piA": r["piA"],
                                  "fidB": r["fidB"], "piB": r["piB"], "type": t,
                                  "surface_bbox_overlap": r["bbox_overlap"], "cover": r["cover"],
                                  "certified": r["certified"], "undecided": r["undecided"],
                                  "coincidence": r["coincidence"], "ms": r["ms"],
                                  "occt_lines": ln, "tag": tag})
    outcsv = os.path.join(RES, "probe_summary.csv")
    if rows:
        with open(outcsv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    with open(os.path.join(RES, "occt_confusion.csv"), "w", newline="", encoding="utf-8") as fh:
        cols = ["part", "sample_class", "method_type", "n_pairs", "occt_has_line",
                "occt_no_line", "occt_error", "agree_rate"]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(occt_rows)
    with open(os.path.join(RES, "occt_conflicts.csv"), "w", newline="", encoding="utf-8") as fh:
        cols = ["part", "fidA", "piA", "fidB", "piB", "type", "surface_bbox_overlap", "cover",
                "certified", "undecided", "coincidence", "ms", "occt_lines", "tag"]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(conflicts)
    return rows, occt_rows, conflicts


# ------------------------------------------------------------ 4. cover check 汇总
def cover_summary():
    files = sorted(glob.glob(os.path.join(RES, "cover", "*_cover.csv")))
    rows = []
    for f in files:
        part = os.path.basename(f).replace("_cover.csv", "")
        tot_pairs = tot_pts = tot_both = pairs_all = pairs_miss = 0
        for r in csv.DictReader(open(f, encoding="utf-8")):
            tot_pairs += 1
            tot_pts += int(r["npoints"])
            tot_both += int(r["inside_both"])
            if int(r["inside_both"]) == int(r["npoints"]):
                pairs_all += 1
            else:
                pairs_miss += 1
        rows.append({"part": part, "pairs_checked": tot_pairs, "points": tot_pts,
                     "points_inside_both": tot_both, "pairs_all_inside": pairs_all,
                     "pairs_with_miss": pairs_miss,
                     "point_cover_rate": round(tot_both / max(1, tot_pts), 4)})
    if rows:
        with open(os.path.join(RES, "cover_summary.csv"), "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    return rows


def main():
    merged = merge_batch()
    repro = repro_dual()
    probe_rows, occt_rows, conflicts = probe_analysis()
    cover = cover_summary()
    stats = {"batch_rows": len(merged), "repro_rows": len(repro),
             "probe_parts": len(probe_rows), "conflicts": len(conflicts),
             "cover_parts": len(cover)}
    with open(os.path.join(RES, "stats.json"), "w", encoding="utf-8") as fh:
        json.dump(stats, fh, ensure_ascii=False, indent=1)
    print(json.dumps(stats, ensure_ascii=False))
    # 控制台摘要
    print("\n== repro dual binary ==")
    for r in repro:
        print("  %-4s %-5s rc=%-3s timeout=%-2s 1D=%-4s 2D=%-3s cop=%-3s und=%-4s certB=%-7s undB=%-8s t=%ss"
              % (r["part"], r["version"], r["rc"], r["timeout"], r["d1"], r["d2"],
                 r["coplanar"], r["undecided"], r["certified_boxes"], r["undecided_boxes"], r["method_s"]))
    print("\n== probe / occt ==")
    for r in probe_rows:
        print("  %-4s pos=%d type%s med=%.1fms max=%.1fms | bbox-ov-cert-disjoint=%d occt_no_line=%d has=%d"
              % (r["part"], r["n_pos"], [r["pos_type%d" % k] for k in range(5)],
                 r["pos_med_ms"], r["pos_max_ms"], r["b_bbox_overlap_cert_disjoint"],
                 r["b_occt_no_line"], r["b_occt_has_line"]))
    print("\n== conflicts: %d ==" % len(conflicts))
    from collections import Counter
    for tag, c in Counter(x["tag"] for x in conflicts).most_common():
        print("  %-40s %d" % (tag, c))


if __name__ == "__main__":
    sys.exit(main())
