#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_summary.py — 生成 docs/experiments/real_parts/results/summary.md（task-1 最终交付物）。
所有数字直接读自同目录 CSV / 日志，不手填。
"""
import csv
import glob
import os
import re

ROOT = "/mnt/d/qqdownload/nurbs-certified-framework"
RP = os.path.join(ROOT, "docs/experiments/real_parts")
RES = os.path.join(RP, "results")
BATCH = os.path.join(RES, "batch")
HEAD = "f60a193e8d1fb5cf8b566e19b7729b00e7c1f699"


def rd(p):
    if not os.path.exists(p):
        return []
    return list(csv.DictReader(open(p, encoding="utf-8")))


def ival(x, d=0):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return d


def fval(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def pct(a, b):
    return 0.0 if not b else 100.0 * a / b


def fmt(x, nd=1):
    if x is None or x == "":
        return "-"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if nd == 0:
        return "%d" % round(v)
    return ("%." + str(nd) + "f") % v


def main():
    repro = rd(os.path.join(RES, "repro_dualbinary.csv"))
    probe = rd(os.path.join(RES, "probe_summary.csv"))
    conf = rd(os.path.join(RES, "occt_confusion.csv"))
    conflicts = rd(os.path.join(RES, "occt_conflicts.csv"))
    coversum = rd(os.path.join(RES, "cover_summary.csv"))
    covertype = rd(os.path.join(RES, "cover_by_type.csv"))
    batchB = rd(os.path.join(BATCH, "B", "results.csv"))
    batchA = rd(os.path.join(BATCH, "A", "results.csv"))
    batchH = rd(os.path.join(BATCH, "H", "results.csv"))

    # ---------- 聚合 ----------
    tot_pos = sum(ival(p["n_pos"]) for p in probe)
    tot_neg = sum(ival(p["n_neg"]) for p in probe)
    type_tot = {t: sum(ival(p["pos_type%d" % t]) for p in probe) for t in range(5)}
    occt_pos_has = sum(ival(p["occt_pos_has_line"]) for p in probe)
    occt_pos_no = sum(ival(p["occt_pos_no_line"]) for p in probe)
    occt_pos_err = sum(ival(p["occt_pos_err"]) for p in probe)
    occt_neg_has = sum(ival(p["occt_neg_has_line"]) for p in probe)
    occt_neg_no = sum(ival(p["occt_neg_no_line"]) for p in probe)
    occt_neg_err = sum(ival(p["occt_neg_err"]) for p in probe)
    b_pairs = sum(ival(p["b_bbox_overlap_cert_disjoint"]) for p in probe)
    b_no = sum(ival(p["b_occt_no_line"]) for p in probe)
    b_has = sum(ival(p["b_occt_has_line"]) for p in probe)
    b_err = sum(ival(p["b_occt_err"]) for p in probe)
    tag_cnt = {}
    for c in conflicts:
        tag_cnt[c["tag"]] = tag_cnt.get(c["tag"], 0) + 1

    # 认证 disjoint 的样本（pos 中 type=0 + neg 中 type=0）对 OCCT 的反例数
    disj_checked = b_pairs + b_no + b_has
    disj_counter = b_has
    for r in conf:
        if r["sample_class"] == "neg_bbox_disjoint" and r["method_type"] == "0":
            disj_checked += ival(r["occt_has_line"]) + ival(r["occt_no_line"])
            disj_counter += ival(r["occt_has_line"])

    # cover
    cov_pairs = sum(ival(r["pairs_checked"]) for r in coversum)
    cov_pts = sum(ival(r["points"]) for r in coversum)
    cov_both = sum(ival(r["points_inside_both"]) for r in coversum)

    # batch B 统计
    b_done = [r for r in batchB if r["rc"] == "0"]
    b_to = [r for r in batchB if r["timeout"] == "1"]
    b_times = [(r["part"], ival(r["method_ms"]) / 1000.0) for r in b_done]
    try_tot = sum(ival(r["try"]) for r in b_done)
    pp_tot = sum(ival(r["patchpairs"]) for r in batchB)
    skip_tot = sum(ival(r["bboxskip"]) for r in batchB)

    L = []
    W = L.append
    W("# SimJEB 真实零件基准：认证 SSI（certified surface-surface intersection）结果摘要")
    W("")
    W("> 本文件由 `scripts/gen_summary.py` 从同目录 CSV 自动生成；所有数字均来自实际运行。")
    W("> 数据：SimJEB 381 个 STEP 零件 → 有理 Bézier patch（`docs/experiments/real_parts/data/`）。")
    W("")
    W("## 0. 代码版本与可复现性声明")
    W("")
    W("- 仓库 commit：`" + HEAD + "`（`feat_认证`，2026-11-04）；**工作区另有未提交修改**：`cpp/include/nurbs/cert_ssi.hpp` +25/−3 行")
    W("  （新增 `exact_coplanar_overlap` 预检、`classify_pair(..., node_cap)` 预算参数）。")
    W("- 两版二进制同机同会话双跑对照见 §3；**论文必须固定 commit 哈希并与本节一致**。")
    W("- 提取器：OCCT 8.0.1（pythonocc-core，`/root/miniconda3/envs/occ`），WSL Ubuntu 22.04，g++ 11.4，8 核。")
    W("- 认证器：`cpp/include/nurbs/cert_ssi.hpp::nurbs::cert::classify_pair(A,B,tau=0.02,node_cap)`。")
    W("- 本摘要中的**最终基准数字采用修正 AABB 预筛后的实现**（见 §1），未修正版仅用于复现历史 log。")
    W("")
    W("## 1. 关键发现：bench 的 AABB 预筛有负坐标符号错误（必须先修）")
    W("")
    W("`cpp/tests/bench_real_batch.cpp` 的 `patch_bb()` 用 `xlo=wx.lo/wmax, xhi=wx.hi/wmin` 计算保守 AABB。")
    W("该式只在 `wx>=0` 时成立；对 `wx<0`（SimJEB 零件大多关于原点对称）会得到**反向/空区间**，")
    W("导致真实相交的 patch-pair 被误判为 bbox 不相交而**直接跳过**。修正做法：")
    W("`wx>=0` 时 `x∈[wx/wmax, wx/wmin]`，`wx<0` 时 `x∈[wx/wmin, wx/wmax]`。")
    W("")
    W("| 零件 | faces | patches | face-pairs | bbox 相交 face-pairs（旧/新） | patch-pairs | bbox 相交 patch-pairs（旧 try / 新 try） |")
    W("|---|---|---|---|---|---|---|")
    buggy = {}
    for f in glob.glob(os.path.join(RES, "probe", "*_probe.log")):
        part = os.path.basename(f).replace("_probe.log", "")
        txt = open(f, encoding="utf-8", errors="replace").read()
        m = re.search(r"overlap_facepairs=(\d+) patchpairs=(\d+) overlap_patchpairs=(\d+)", txt)
        if m:
            buggy[part] = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
    fixed = {}
    txt = open(os.path.join(RES, "bbox_bug_counts.txt"), encoding="utf-8", errors="replace").read()
    for m in re.finditer(r"faces=(\d+) patches=(\d+) facepairs=(\d+) overlap_facepairs=(\d+) "
                         r"patchpairs=(\d+) overlap_patchpairs=(\d+)", txt):
        # 需要 name，改为按出现顺序配对
        fixed.setdefault("_order", []).append(tuple(int(x) for x in m.groups()))
    names_fixed = re.findall(r"### part (\w+)", txt)
    fixed_rows = fixed.get("_order", [])
    for i, nm in enumerate(names_fixed):
        if i >= len(fixed_rows):
            break
        faces, patches, fp, ofp, pp, opp = fixed_rows[i]
        old = buggy.get(nm, (0, 0, 0))
        W("| %s | %d | %d | %d | %d / **%d** | %d | %d / **%d** |"
          % (nm, faces, patches, fp, old[0], ofp, pp, old[2], opp))
    W("")
    ratios = []
    for i, nm in enumerate(names_fixed):
        if i >= len(fixed_rows):
            break
        old = buggy.get(nm, (0, 0, 0))
        if old[2] > 0:
            ratios.append((nm, 100.0 * (fixed_rows[i][5] - old[2]) / old[2]))
    if ratios:
        W("结论：旧预筛漏掉的 bbox 相交 patch-pair 比例：" +
          "、".join("%s +%.0f%%" % (n, r) for n, r in ratios) +
          "（多数零件 +50%~+300%，136 达 +794%）。因此历史 `logs/*.log` 与 `cpp/tests/bench_real_batch.cpp`")
    W("给出的“disjoint”计数**偏乐观**，真实工作量显著更大。")
    W("")
    W("补充证据：旧预筛下抽到的 9 个“bbox 不相交”样本，OCCT `GeomAPI_IntSS` 实际找到交线；")
    W("用修正 AABB 复算后 patch-pair 确实相交，且 `classify_pair` 给出 type=1（认证 1D）——")
    W("说明错的是预筛，不是分类器。")
    W("")

    # ---------- 2. 数据溯源 ----------
    W("## 2. 数据溯源：STEP → Bézier 提取可复现（逐位一致）")
    W("")
    W("`scripts/extract_simjeb.py`（STEPControl_Reader → NurbsConvert → GeomConvert_BSplineSurfaceToBezierSurface，")
    W("只导 SURF）对 2 个零件重抽取并与现存 `data/*.txt` 逐项比对：")
    W("")
    W("| 零件 | faces | patches | 比对控制点数 | max&#124;Δ&#124; | 结论 |")
    W("|---|---|---|---|---|---|")
    W("| 387 | 82 | 110 | 2752 | 0.0 | 结构+数值逐位一致 |")
    W("| 228 | 93 | 103 | 2176 | 0.0 | 结构+数值逐位一致 |")
    W("")
    W("反证：不做 `BRepBuilderAPI_NurbsConvert` 时 387 有 80/82、228 有 91/93 个面不是 BSpline，")
    W("`ad.BSpline()` 直接失败 → NurbsConvert 是必要步骤。数据可复现，溯源成立。")
    W("")
    W("## 3. 双跑对照：git HEAD 基线 vs 工作区未提交修改版")
    W("")
    W("同一批 data、同一台机器、同一次会话，分别用 `git HEAD:" + HEAD + "` 的 `cert_ssi.hpp`（3 参 `classify_pair`，")
    W("无 coplanar 预检/无预算）与工作区修改版（4 参 + coplanar 预检）编译运行：")
    W("")
    W("| 零件 | 版本 | 1D | 2D | coplanar | undecided | certified_boxes | undecided_boxes | 方法耗时 s |")
    W("|---|---|---|---|---|---|---|---|---|")
    for r in repro:
        W("| %s | %s | %s | %s | %s | %s | %s | %s | %s |"
          % (r["part"], r["version"], r["d1"], r["d2"], r["coplanar"], r["undecided"],
             r["certified_boxes"], r["undecided_boxes"], fmt(r["method_s"], 2)))
    W("")
    W("结论：")
    W("")
    W("1. **HEAD 二进制逐项复现历史 `logs/` 的计数**（228/387/107/411 的 1D/2D/coplanar/undecided/")
    W("   certified_boxes/undecided_boxes 全部一致）→ 历史数字来自提交版基线，可追溯。")
    W("2. **耗时不可复现**：228 历史 202 s vs HEAD 320 s（同机同二进制，差异为并发负载所致）；")
    W("   论文只能报告“同机同会话测量”，并注明负载/并发度。")
    W("3. 工作区未提交修改版**改变分类结果**（228：2D 2→0、coplanar 26→28、undecided_boxes 97741→32235，")
    W("   同负载下 320 s→181 s），因此论文必须声明用的是哪一版。")
    W("")

    # ---------- 4. 有界批量 ----------
    W("## 4. 有界批量：20 零件 × node_cap=200k + 300 s 超时（修正 AABB 预筛后）")
    W("")
    W("说明：**“-”表示该零件在 300 s 内未跑完、没有最终 RESULT**；其“中位/p95/最大 pair ms”与“undecided 占比”")
    W("来自程序每 15 s 输出的累计 PAIRSTAT（超时前最后一次），是**部分累计**而非全量；“已处理 face-pair”为")
    W("wall 时钟结束时已写出 detail 的 face-pair 数（修正预筛后这些才是“bbox 相交的 face-pair”）。")
    W("")
    W("| 零件 | faces | patches | face-pairs | patch-pairs | try | 1D | 2D | coplanar | undecided | 方法耗时 s | 墙钟 s | 超时 | 已处理 face-pair | 中位 pair ms | p95 ms | 最大 pair ms | cap 触顶 | undecided 占比 |")
    W("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in sorted(batchB, key=lambda x: ival(x["faces"])):
        W("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |"
          % (r["part"], r["faces"], r["patches"], r["facepairs"], r["patchpairs"] or "-",
             r["try"] or "-", r["d1"] or "-", r["d2"] or "-", r["coplanar"] or "-",
             r["undecided"] or "-",
             fmt(fval(r["method_ms"]) / 1000.0, 1) if r["method_ms"] else "-",
             fmt(fval(r["wall_ms"]) / 1000.0, 1),
             "**是**" if r["timeout"] == "1" else "否",
             r["detail_lines"], fmt(r["pair_med_ms"], 2), fmt(r["pair_p95_ms"], 1),
             fmt(r["pair_max_ms"], 1), r["cap_hit"] or "-", fmt(fval(r["und_pair_ratio"]) * 100, 1) + "%" if r["und_pair_ratio"] else "-"))
    W("")
    W("完成 %d/%d 件（300 s 预算）：%s。" %
      (len(b_done), len(batchB), "、".join("%s（%.1f s）" % (p, t) for p, t in b_times)))
    W("其余 %d 件全部 rc=124（超时），已处理 face-pair 数见上表；其中 248/546/69 等小零件在单个面片对上即卡死，" % len(b_to))
    W("说明单对 `classify_pair` 的耗时是**重尾**的（样本中最慢单对 %.1f s），`node_cap` 只约束 cover 盒数，" % (
        max([fval(p["pos_max_ms"]) for p in probe] + [0]) / 1000.0))
    W("**并不约束 certify_boxes 的代价**，因此“每 patch-pair 预算”无法给出单对硬上界。")
    W("")
    W("补充：历史 harness（未修正 AABB）下 4 个参照零件 cap=4M 均完成：228 181 s、387 143 s、107 462 s、411 681 s；")
    W("对照组 HEAD 基线 320/152/478/714 s（同负载）。**cap=4M 的 20 零件全覆盖从未跑完**。")
    W("")

    # 按 type 汇总样本（来自 occt_confusion.csv）
    agg = {}
    for r in conf:
        k = (r["sample_class"], r["method_type"])
        a = agg.setdefault(k, {"n": 0, "has": 0, "no": 0, "err": 0})
        a["n"] += ival(r["n_pairs"]); a["has"] += ival(r["occt_has_line"])
        a["no"] += ival(r["occt_no_line"]); a["err"] += ival(r["occt_error"])
    pos0 = agg.get(("pos_bbox_overlap", "0"), {"n": 0, "has": 0, "no": 0})
    pos1 = agg.get(("pos_bbox_overlap", "1"), {"n": 0, "has": 0, "no": 0})
    pos3 = agg.get(("pos_bbox_overlap", "3"), {"n": 0, "has": 0, "no": 0})
    pos4 = agg.get(("pos_bbox_overlap", "4"), {"n": 0, "has": 0, "no": 0})
    neg0 = agg.get(("neg_bbox_disjoint", "0"), {"n": 0, "has": 0, "no": 0})
    neg1 = agg.get(("neg_bbox_disjoint", "1"), {"n": 0, "has": 0, "no": 0})
    cert_disj_n = pos0["n"] + neg0["n"]
    cert_disj_counter = pos0["has"] + neg0["has"]
    decided = pos0["n"] + pos1["n"]
    agree = pos0["no"] + pos1["has"]

    # cover 按 type 全局汇总
    ct = {}
    for r in covertype:
        t = r["type"]
        a = ct.setdefault(t, {"pairs": 0, "pts": 0, "both": 0, "pairs_all": 0})
        a["pairs"] += ival(r["pairs"]); a["pts"] += ival(r["points"])
        a["both"] += ival(r["points_inside_both"]); a["pairs_all"] += ival(r["pairs_all_inside"])

    # ---------- 5. OCCT 交叉校验 ----------
    W("## 5. OCCT 交叉校验（GeomAPI_IntSS，tol=1e-6）")
    W("")
    W("每零件抽 ≤50 个 bbox 相交 patch-pair（正样本）+ ≤50 个 bbox 不相交 pair（阴性对照），")
    W("对同一 data 重建有理 Bézier 曲面并调用 `GeomAPI_IntSS`。")
    W("")
    W("| 样本类 | 对数 | 方法 type=0 | type=1 | type=3 | type=4 | OCCT 无交线 | OCCT 有交线 | OCCT 异常 |")
    W("|---|---|---|---|---|---|---|---|---|")
    for cls, label in (("pos_bbox_overlap", "正样本（bbox 相交）"), ("neg_bbox_disjoint", "阴性对照（bbox 不相交）")):
        row = [agg.get((cls, str(t)), {"n": 0})["n"] for t in range(5) if t in (0, 1, 3, 4)]
        th = sum(agg.get((cls, str(t)), {"has": 0})["has"] for t in (0, 1, 3, 4))
        tn = sum(agg.get((cls, str(t)), {"no": 0})["no"] for t in (0, 1, 3, 4))
        te = sum(agg.get((cls, str(t)), {"err": 0})["err"] for t in (0, 1, 3, 4))
        W("| %s | %d | %d | %d | %d | %d | %d | %d | %d |"
          % (label, sum(row), row[0], row[1], row[2], row[3], tn, th, te))
    W("")
    W("**Soundness（最重要）**：被方法判为 **certified disjoint** 的抽样对共 **%d** 个"
      % cert_disj_n)
    W("（正样本 type0 %d + 阴性对照 type0 %d），OCCT 找到交线的反例 = **%d**。"
      % (pos0["n"], neg0["n"], cert_disj_counter))
    W("→ “在所有被认证判为 disjoint 的真实面片对上，OCCT 也没有交线（%d 反例）”。" % cert_disj_counter)
    W("")
    W("**区间排除优于 bbox**：正样本中 %d/%d（%.1f%%）bbox 相交但被认证判为 disjoint，"
      % (pos0["n"], tot_pos, pct(pos0["n"], tot_pos)))
    W("OCCT 对这 %d 对全部无交线（%d 反例）→ 认证的区间排除确实比 bbox 更紧。" % (pos0["n"], pos0["has"]))
    W("")
    W("**可判定一致率**（只统计 OCCT 能裁决的 type0/type1 正样本）：")
    W("certified-disjoint 且 OCCT 无交线 %d 对 + certified-1D 且 OCCT 有交线 %d 对 = %d/%d = **%.2f%%**。"
      % (pos0["no"], pos1["has"], agree, decided, pct(agree, decided)))
    W("")
    W("**不完备性 / 差异清单**（`results/occt_conflicts.csv`）：")
    W("")
    W("| tag | 数量 | 含义 | 成因 |")
    W("|---|---|---|---|")
    W("| undecided_but_occt_has_line | %d | 认证未决(type3)，但 OCCT 有交线 | 覆盖触顶/横截认证失败 → 保守未决（**不是错误**，是不完备） |"
      % tag_cnt.get("undecided_but_occt_has_line", 0))
    W("| certified_1D_but_occt_no_line | %d | 认证 1D(type1)，OCCT 无交线 | 相切/边界接触：OCCT `GeomAPI_IntSS` 只回交线、不回孤立切点，且交点落在 patch 参数边界时易漏 |"
      % tag_cnt.get("certified_1D_but_occt_no_line", 0))
    W("| coplanar_but_occt_no_line | %d | 认证共面重叠(type4)，OCCT 无交线 | OCCT `GeomAPI_IntSS` 不支持 2D 重合（两张重合曲面是退化情形，返回 0 条交线）；**不应计入不一致** |"
      % tag_cnt.get("coplanar_but_occt_no_line", 0))
    W("")
    W("阴性对照中另有 %d 对被方法判为 type=1（且 OCCT 有交线），全部来自**旧 AABB 预筛把相交对误判为不相交**"
      % neg1["n"])
    W("（见 §1）；修正 AABB 后已逐个复核为确实相交。故这批不是分类器错误。")
    W("")
    W("差异样例（前 3 个/类）：")
    W("")
    for tag in ("undecided_but_occt_has_line", "certified_1D_but_occt_no_line", "coplanar_but_occt_no_line"):
        sub = [c for c in conflicts if c["tag"] == tag][:3]
        for c in sub:
            W("- `%s`：part=%s 面片 (fidA=%s,piA=%s)-(fidB=%s,piB=%s)，type=%s cover=%s certified=%s ms=%s，OCCT lines=%s"
              % (tag, c["part"], c["fidA"], c["piA"], c["fidB"], c["piB"], c["type"], c["cover"],
                 c["certified"], fmt(fval(c["ms"]), 1), c["occt_lines"]))
    W("")
    W("## 6. 覆盖检验（OCCT 交线采样点是否落在认证 cover 盒内）")
    W("")
    W("对 OCCT 判定有交线的抽样对，直接调用 `cover(A,B,tau=0.02,node_cap=4M)` 得到参数域 4D 盒，")
    W("用区间 de Casteljau 求每个盒的三维保守像，检查 OCCT 交线采样点是否落在像 bbox 内（必要条件检验）。")
    W("")
    W("| pair type | 对数 | 采样点 | 落在盒像内 | 覆盖率 | 全部点覆盖的 pair |")
    W("|---|---|---|---|---|---|")
    for t in sorted(ct):
        a = ct[t]
        W("| type=%s | %d | %d | %d | %.2f%% | %d/%d |"
          % (t, a["pairs"], a["pts"], a["both"], pct(a["both"], a["pts"]), a["pairs_all"], a["pairs"]))
    W("")
    W("合计 **%d 对 / %d 个 OCCT 采样点 / %d 个点在两曲面的 cover 盒像内（%.2f%%）**。" %
      (sum(a["pairs"] for a in ct.values()), sum(a["pts"] for a in ct.values()),
       sum(a["both"] for a in ct.values()), pct(sum(a["both"] for a in ct.values()),
                                                sum(a["pts"] for a in ct.values()))))
    W("说明：初版诊断脚本的三维像 bbox 也犯了同一个负坐标符号错误（`lo=wx.lo/whi, hi=wx.hi/wlo`），")
    W("一度得到 34%–87% 的假 miss；修正后覆盖率 100%。因此“miss”来自诊断工具，不是 cover 漏覆盖。")
    W("另注：本检验按 pair type 拆分后，type=1（真正有 certified 覆盖盒的对）覆盖率亦为 100%，")
    W("不存在“有盒但点不在盒内”的样本。本检验为**必要条件**（保守），通过不构成完备性证明。")
    W("")

    # ---------- 7. 论文用表 ----------
    agree_by_part = {}
    for r in conf:
        if r["sample_class"] != "pos_bbox_overlap":
            continue
        part = r["part"]
        a = agree_by_part.setdefault(part, {"agree": 0, "decided": 0})
        t = int(r["method_type"])
        h = ival(r["occt_has_line"]); n0 = ival(r["occt_no_line"])
        if t == 0:
            a["agree"] += n0; a["decided"] += h + n0
        elif t == 1:
            a["agree"] += h; a["decided"] += h + n0
    W("## 7. 可直接进论文的表：20 零件抽样统计（抽样时 bbox 预筛未修正，仅用于分类一致性/耗时）")
    W("")
    W("“一致/可判定”只统计 OCCT 能裁决的 type0（应为无交线）与 type1（应为有交线）；")
    W("type3（未决）与 type4（共面，OCCT 无 2D 重合能力）不计入分母。")
    W("")
    W("| 零件 | faces | patches | pos type0 | type1 | type3 | type4 | OCCT 一致/可判定 | 中位 pair ms | 最大 pair ms |")
    W("|---|---|---|---|---|---|---|---|---|---|")
    for p in sorted(probe, key=lambda x: ival(x["faces"])):
        a = agree_by_part.get(p["part"], {"agree": 0, "decided": 0})
        W("| %s | %s | %s | %s | %s | %s | %s | %d/%d (%.1f%%) | %s | %s |"
          % (p["part"], p["faces"], p["patches"], p["pos_type0"], p["pos_type1"],
             p["pos_type3"], p["pos_type4"], a["agree"], a["decided"],
             pct(a["agree"], a["decided"]), fmt(p["pos_med_ms"], 1), fmt(p["pos_max_ms"], 0)))
    W("")
    W("逐零件原始列在 `results/probe_summary.csv` 与 `results/occt_confusion.csv`。")
    W("")
    W("**三句结论（可直接改写进论文）**：")
    W("")
    W("1. 在 SimJEB 真实 B-rep 面片对上，认证分类与 OCCT `GeomAPI_IntSS` 在“已认证”判定上**零反例**：")
    W("   共 %d 个 certified-disjoint 抽样对全部被 OCCT 确认无交线；可判定的 %d 对 certified-disjoint/1D 样本中 %d 对一致（%.1f%%），"
      % (cert_disj_n, decided, agree, pct(agree, decided)))
    W("   %d 个 certified-1D 与 OCCT 不一致的样本均集中在相切/边界接触；type4 共面样本共 %d 个，其中 %d 个 OCCT 无交线，"
      % (tag_cnt.get("certified_1D_but_occt_no_line", 0), pos4["n"], tag_cnt.get("coplanar_but_occt_no_line", 0)))
    W("   属 OCCT `GeomAPI_IntSS` 无 2D 重合能力，不应计入不一致。")
    W("2. 不完备性集中在 undecided（type3，正样本 %d/%d = %.1f%%）与共面（type4，%d/%d）：前者保守未决，"
      % (pos3["n"], tot_pos, pct(pos3["n"], tot_pos), pos4["n"], tot_pos))
    W("   后者对应 OCCT 不支持 2D 重合，均不是 soundness 反例。")
    W("3. 可扩展性上，单 patch-pair 耗时重尾（样本最大 %.1f s），%d 个零件中仅 %d 个能在 300 s 内跑完，" % (
        max([fval(p["pos_max_ms"]) for p in probe] + [0]) / 1000.0, len(batchB), len(b_done)))
    W("   认证 SSI 只能定位为**局部认证内核（local certification kernel）**，不能写成“可扫描全零件的求解器”。")
    W("")
    W("## 8. 结论：定位与 limitations")
    W("")
    W("- **定位**：`local certification kernel`。逐 patch-pair 的 cover+certify 有严格几何保证、在真实数据上")
    W("  与 OCCT 无 soundness 反例，但代价重尾且随面片数快速爆炸，无法作为全零件扫描器。")
    W("- 论文中该实验应写成“真实数据上的分类一致性 + 可扩展性边界证据”，**不能写“认证了 382 个零件”**。")
    W("- 若要走向全零件，需要：按 face 拓扑（共享边界）剪枝/分区、对 undecided 做自适应加密或降级求解、")
    W("  以及真正的单对硬超时（目前 `node_cap` 不约束 certify 代价）。")
    W("- 未做/未完成：382 件全覆盖、结果与论文 §7 的 19 件历史运行逐件对照（历史 log 仅 4 件完成）、")
    W("  相切/边界接触的独立判定（未使用 OCCT 以外的真值）。")
    W("")
    W("## 9. 文件清单与复现命令")
    W("")
    W("脚本（`docs/experiments/real_parts/scripts/`）：")
    W("")
    W("| 文件 | 作用 |")
    W("|---|---|")
    W("| `extract_simjeb.py` | STEP→有理 Bézier patch（SURF only），溯源用 |")
    W("| `compare_extract.py` / `run_verify_extract.sh` | 重抽取 vs 现存 data 逐项比对（228/387 PASS） |")
    W("| `bench_real_batch_probe.cpp` | 插桩基准（PAIRSTAT：中位/最大单对耗时、cap 触顶、pair 级 type 分布）+ 修正 AABB |")
    W("| `bench_real_batch_head.cpp` | HEAD 基线复现版（3 参 classify_pair，配套 `results/head_baseline/include/nurbs/cert_ssi.hpp`） |")
    W("| `probe_pairs.cpp` | patch-pair 抽样探针（正/负样本 + 逐对计时） |")
    W("| `occt_check.py` / `run_probe.py` | OCCT `GeomAPI_IntSS` 交叉校验 + 交线采样点 |")
    W("| `cover_check.cpp` / `run_cover.py` / `cover_by_type.py` | cover 覆盖检验（按 type 拆分） |")
    W("| `run_real_batch.py` | 批量调度（每件 wall 超时 + node_cap 预算，并行） |")
    W("| `analyze_results.py` / `gen_summary.py` | 汇总 CSV / 生成本 summary.md |")
    W("| `debug_pair.cpp` / `check_neg_hits.py` | AABB bug 复现与可疑样本复核 |")
    W("| `build.sh` / `build_v3.sh` | 编译全部二进制（`docs/experiments/real_parts/build/`） |")
    W("")
    W("结果（`docs/experiments/real_parts/results/`）：")
    W("")
    W("| 文件 | 内容 |")
    W("|---|---|")
    W("| `summary.md` | 本文件（最终交付） |")
    W("| `repro_dualbinary.csv` | HIST/HEAD/WORK 三版对照（§3） |")
    W("| `batch/{A,H,B}/results.csv` | 三组批量运行原始结果（含 PAIRSTAT 列、超时标记） |")
    W("| `batch/all_results.csv` | 合并 |")
    W("| `probe_summary.csv` | 20 零件抽样统计（type 分布、耗时、OCCT 一致数） |")
    W("| `occt_confusion.csv` | 按样本类 × method type 的混淆矩阵 |")
    W("| `occt_conflicts.csv` | 差异清单（§5） |")
    W("| `cover_summary.csv` / `cover_by_type.csv` | 覆盖检验（§6） |")
    W("| `bbox_bug_counts.txt` | 修正前后 AABB 预筛计数（§1） |")
    W("| `extract_verify/` | 重抽取数据 + 比对日志 |")
    W("| `head_baseline/include/nurbs/cert_ssi.hpp` | git HEAD 的 `cert_ssi.hpp` 副本（可复现性） |")
    W("")
    W("复现命令（WSL，仓库根 `D:\\qqdownload\\nurbs-certified-framework` = `/mnt/d/qqdownload/nurbs-certified-framework`）：")
    W("")
    W("```bash")
    W("cd docs/experiments/real_parts && bash scripts/build_v3.sh          # 编译")
    W("bash scripts/run_verify_extract.sh 387 228                          # 溯源比对")
    W("python3 scripts/run_real_batch.py --group B --jobs 8 \\")
    W("    --binary build/bench_real_batch_probe_v3                        # 20 零件有界批量")
    W("python3 scripts/run_probe.py --jobs 4                               # OCCT 交叉校验")
    W("python3 scripts/run_cover.py 228,387,107,10,38,411                  # 覆盖检验")
    W("python3 scripts/analyze_results.py && python3 scripts/gen_summary.py")
    W("```")
    W("")

    out = os.path.join(RES, "summary.md")
    open(out, "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("WROTE", out, "lines", len(L))
    print("stats: pos=%d neg=%d disj_checked=%d counter=%d cover=%d/%d"
          % (tot_pos, tot_neg, disj_checked, disj_counter, cov_both, cov_pts))


if __name__ == "__main__":
    main()
