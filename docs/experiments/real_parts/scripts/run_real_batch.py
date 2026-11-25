#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_real_batch.py — 真实零件有界批量认证（每 patch-pair 预算 + 每件墙钟超时）。

调用 cpp/build/bench_real_batch_linux，按 (part, tau, node_cap, timeout) 逐件运行，
解析 RESULT 行写入 CSV；超时件记录已完成 face-pair 进度（detail 行数 / 全部 face-pair 数）。

用法（WSL）：
  python3 run_real_batch.py --group A|B|C --jobs 7
输出：
  results/batch/<group>/<part>_cap<cap>.log
  results/batch/<group>/results.csv
"""
import argparse
import csv
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

ROOT = "/mnt/d/qqdownload/nurbs-certified-framework"
RP = os.path.join(ROOT, "docs/experiments/real_parts")
BIN = os.path.join(ROOT, "cpp/build/bench_real_batch_linux")
SHM = "/dev/shm"

# 20 件，覆盖 manifest faces 分档 56/93/120/282/499/1000+，含历史已完成的 228/387/411/107
PARTS_B = ["248", "115", "69", "387", "228", "411", "506", "70", "538", "570",
           "113", "205", "546", "107", "534", "200", "614", "38", "10", "136"]

# 预算曲线：同一零件 x 不同 node_cap
PARTS_C = ["228", "107", "411"]

GROUPS = {
    # A: 默认预算（node_cap=4M，即几乎无预算限制），工作区修改版，与历史 logs 对照
    "A": [("228", 0.02, 4000000, 420), ("387", 0.02, 4000000, 420),
          ("107", 0.02, 4000000, 520), ("411", 0.02, 4000000, 700)],
    # H: 同上零件，git HEAD 基线二进制（3 参 classify_pair、无 coplanar 预检/无预算）
    "H": [("228", 0.02, 4000000, 420), ("387", 0.02, 4000000, 420),
          ("107", 0.02, 4000000, 560), ("411", 0.02, 4000000, 760)],
    # B: 有界预算 node_cap=200k，20 件全覆盖（插桩二进制，行缓冲 + 周期 PAIRSTAT）
    "B": [(p, 0.02, 200000, 300) for p in PARTS_B],
    # C: 预算曲线 node_cap = 50k / 500k（插桩二进制）
    "C": [(p, 0.02, cap, 300) for p in PARTS_C for cap in (50000, 500000)],
}

RESULT_RE = re.compile(
    r"RESULT faces=(\d+) facepairs=(\d+) \| patchpairs=(\d+) bboxskip=(\d+) try=(\d+) \| "
    r"disjoint=(\d+) \| 1D=(\d+) \(comps=(\d+)\) \| 2D=(\d+) \| coplanar=(\d+) \| undecided=(\d+)")
BOXES_RE = re.compile(r"certified_boxes=(\d+) undecided_boxes=(\d+) time=([\d.]+) ms")
PARSED_RE = re.compile(r"parsed (\d+) faces \((\d+) patches\), tau=([\d.]+)")

FIELDS = ["group", "part", "tau", "node_cap", "timeout_s", "rc", "timeout",
          "faces", "patches", "facepairs", "patchpairs", "bboxskip", "try",
          "disjoint", "d1", "comps", "d2", "coplanar", "undecided",
          "certified_boxes", "undecided_boxes", "method_ms", "wall_ms",
          "detail_lines", "progress_pct",
          "pair_n", "pair_med_ms", "pair_p95_ms", "pair_p99_ms", "pair_max_ms",
          "cap_hit", "cover_gt2500", "pair_t0", "pair_t1", "pair_t2", "pair_t3", "pair_t4",
          "und_pair_ratio", "t0_ms", "t1_ms", "t2_ms", "t3_ms", "t4_ms", "log"]

PAIRSTAT_RE = re.compile(
    r"PAIRSTAT tau=[\d.]+ cap=(\d+) n=(\d+) med_ms=([\d.]+) p95_ms=([\d.]+) p99_ms=([\d.]+) max_ms=([\d.]+) "
    r"cap_hit=(\d+) cover_gt2500=(\d+) pair_t0=(\d+) pair_t1=(\d+) pair_t2=(\d+) pair_t3=(\d+) pair_t4=(\d+) "
    r"und_pair_ratio=([\d.]+) t0_ms=([\d.]+) t1_ms=([\d.]+) t2_ms=([\d.]+) t3_ms=([\d.]+) t4_ms=([\d.]+)")


def run_one(job):
    group, name, tau, cap, to, binary = job
    logdir = os.path.join(RP, "results/batch", group)
    os.makedirs(logdir, exist_ok=True)
    tag = "%s_cap%d" % (name, cap)
    log = os.path.join(logdir, tag + ".log")
    det = os.path.join(SHM, "rp_%s_%s.detail" % (group, tag))
    if os.path.exists(det):
        os.remove(det)
    cmd = ["timeout", str(to), binary, os.path.join(RP, "data", name + ".txt"),
           str(tau), str(cap), det]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=to + 60)
        out = p.stdout + p.stderr
        rc = p.returncode
    except subprocess.TimeoutExpired as e:
        out = (e.stdout or "") + (e.stderr or "")
        rc = 124
    wall_ms = (time.time() - t0) * 1000.0
    with open(log, "w", encoding="utf-8") as fh:
        fh.write(out)
    row = {k: "" for k in FIELDS}
    row.update({"group": group, "part": name, "tau": tau, "node_cap": cap,
                "timeout_s": to, "rc": rc, "log": log})
    m = PARSED_RE.search(out)
    nf = 0
    if m:
        nf = int(m.group(1))
        row["faces"] = int(m.group(1)); row["patches"] = int(m.group(2))
    m = RESULT_RE.search(out)
    if m:
        g = [int(x) for x in m.groups()]
        (row["faces"], row["facepairs"], row["patchpairs"], row["bboxskip"], row["try"],
         row["disjoint"], row["d1"], row["comps"], row["d2"], row["coplanar"],
         row["undecided"]) = g
    m = BOXES_RE.search(out)
    if m:
        row["certified_boxes"] = int(m.group(1))
        row["undecided_boxes"] = int(m.group(2))
        row["method_ms"] = float(m.group(3))
    m = PAIRSTAT_RE.search(out)
    if m:
        (row["node_cap"], row["pair_n"], row["pair_med_ms"], row["pair_p95_ms"],
         row["pair_p99_ms"], row["pair_max_ms"], row["cap_hit"], row["cover_gt2500"],
         row["pair_t0"], row["pair_t1"], row["pair_t2"], row["pair_t3"], row["pair_t4"],
         row["und_pair_ratio"], row["t0_ms"], row["t1_ms"], row["t2_ms"], row["t3_ms"],
         row["t4_ms"]) = [float(x) if "." in x else int(x) for x in m.groups()]
    row["wall_ms"] = round(wall_ms, 1)
    nlines = 0
    if os.path.exists(det):
        with open(det, "r", encoding="utf-8", errors="replace") as fh:
            nlines = sum(1 for _ in fh)
    row["detail_lines"] = nlines
    row["timeout"] = 1 if rc == 124 else 0
    if not row["facepairs"] and nf:
        row["facepairs"] = nf * (nf - 1) // 2
    if row["facepairs"] and nf:
        row["progress_pct"] = round(100.0 * nlines / row["facepairs"], 2)
    print("done %-3s %-4s cap=%-8d rc=%d wall=%7.1fms detail=%d"
          % (group, name, cap, rc, wall_ms, nlines), flush=True)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True, choices=sorted(GROUPS))
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--binary", default=BIN)
    args = ap.parse_args()
    jobs = [(args.group, p, tau, cap, to, args.binary) for (p, tau, cap, to) in GROUPS[args.group]]
    os.makedirs(os.path.join(RP, "results/batch", args.group), exist_ok=True)
    rows = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for r in ex.map(run_one, jobs):
            rows.append(r)
    outcsv = os.path.join(RP, "results/batch", args.group, "results.csv")
    with open(outcsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("WROTE", outcsv, flush=True)


if __name__ == "__main__":
    sys.exit(main())
