#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Assemble the cross-baseline comparison CSVs from the raw result files:
  results/sisl_ssi.csv            (C++ SISL harness)
  results/occt_ssi_raw.txt        (pythonocc OCCT harness)
  results/certified_cover_probe.txt (interval-cover certified experiment)
  results/matrix_baseline_raw.txt (minimal matrix-representation reproduction)
No numbers are invented here: every value is parsed from those files.
"""
import csv
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")


def read_lines(name):
    with open(os.path.join(RES, name), encoding="utf-8", errors="replace") as f:
        return [ln.rstrip("\n") for ln in f]


def parse_sisl():
    """fixture -> {epsge: dict}"""
    out = {}
    for ln in read_lines("sisl_ssi.csv"):
        if not ln or ln.startswith("fixture"):
            continue
        # CSV with quoted truth fields
        parts = next(csv.reader([ln]))
        if len(parts) < 15:
            continue
        fx, epsge = parts[0], parts[1]
        out.setdefault(fx, {})[epsge] = {
            "stat": int(parts[2]),
            "points": int(parts[3]),
            "curves": int(parts[4]),
            "guide": int(parts[5]),
            "traced": int(parts[6]),
            "trace_err": int(parts[8]),
            "res_point": float(parts[9]),
            "res_guide": float(parts[10]),
            "res_trace": float(parts[11]),
            "t_topo": float(parts[12]),
            "t_trace": float(parts[13]),
            "truth": parts[14],
        }
    return out


def parse_occt():
    out = {}
    for ln in read_lines("occt_ssi_raw.txt"):
        if not ln or ln.startswith("fixture") or ln.startswith("OCCT_VERSION"):
            continue
        parts = next(csv.reader([ln]))
        if len(parts) < 10:
            continue
        fx, method, tol = parts[0], parts[1], parts[2]
        out.setdefault(fx, {})[tol] = {
            "method": method,
            "done": int(parts[3]),
            "curves": int(parts[4]),
            "points": int(parts[5]),
            "samples": int(parts[6]),
            "residual": float(parts[7]),
            "time_ms": float(parts[8]),
            "truth": parts[9],
        }
    return out


def parse_certified():
    out = {}
    for ln in read_lines("certified_cover_probe.txt"):
        if not ln or ln.startswith("#") or ln.startswith("fixture"):
            continue
        parts = next(csv.reader([ln]))
        if len(parts) < 6:
            continue
        out[parts[0]] = {
            "excl": int(parts[1]),
            "cover": int(parts[2]),
            "subdiv": int(parts[3]),
            "time_ms": float(parts[4]),
            "truth": parts[5],
        }
    return out


def parse_matrix():
    """Return (curve_surface_rows, certified_rows, ss_rows)."""
    cs, cert, ss = [], [], []
    sec = None
    for ln in read_lines("matrix_baseline_raw.txt"):
        if ln.startswith("==="):
            sec = ln
            continue
        if not ln or ln.startswith("#") or ln.startswith("case,") or \
                ln.startswith("fixture,"):
            continue
        parts = next(csv.reader([ln]))
        if sec and "(1)" in sec:
            cs.append(parts)
        elif sec and "(1b)" in sec:
            cert.append(parts)
        elif sec and "(2)" in sec:
            ss.append(parts)
    return cs, cert, ss


def main():
    sisl = parse_sisl()
    occt = parse_occt()
    cert = parse_certified()
    mcs, mcert, mss = parse_matrix()

    # ---------------- results/ssi_baselines_compare.csv ----------------
    notes = {
        "A": "transverse curve; all three methods agree",
        "B": "no intersection; all agree",
        "C": "corner tangency (isolated point); OCCT reports nothing (its "
             "GeomAPI_IntSS public result has curves only), SISL reports 1 point",
        "D": "degenerate tangential curve; OCCT reports 0 curves (misses the "
             "whole curve), SISL fragments it into many micro-curves + points "
             "with s1310 errors, certified cover covers the diagonal",
        "E": "bicubic x bicubic interior curve; SISL and OCCT agree",
        "F": "rational cylinder x line; surface/curve (GeomAPI_IntCS)",
    }
    truth = {
        "A": "1 curve: diagonal s=t",
        "B": "0 intersections",
        "C": "1 point: (0,0) corner tangency",
        "D": "1 tangential curve: diagonal s=t",
        "E": "1 interior curve: u^3+v^3=0.1",
        "F": "1 point: (0.866025,0.5,0.5)",
    }
    rows = []
    for fx in ["A", "B", "C", "D", "E", "F"]:
        s = sisl[fx]["1.0e-06"]
        o = occt[fx]["1e-06"]
        c = cert.get(fx)
        rows.append({
            "fixture": fx,
            "truth": truth[fx],
            "SISL_epsge": "1e-6",
            "SISL_points": s["points"],
            "SISL_curves": s["curves"],
            "SISL_trace_errors": s["trace_err"],
            "SISL_residual": f"{max(s['res_point'], s['res_guide'], s['res_trace']):.3e}",
            "SISL_time_ms": f"{s['t_topo'] + s['t_trace']:.3f}",
            "OCCT_tol": "1e-6",
            "OCCT_curves": o["curves"],
            "OCCT_points": o["points"],
            "OCCT_residual": f"{o['residual']:.3e}",
            "OCCT_time_ms": f"{o['time_ms']:.3f}",
            "certified_param_tol": "1e-3" if c else "n/a",
            "certified_cover_boxes": c["cover"] if c else "n/a",
            "certified_excluded_boxes": c["excl"] if c else "n/a",
            "certified_time_ms": f"{c['time_ms']:.2f}" if c else "n/a",
            "note": notes[fx],
        })
    out = os.path.join(RES, "ssi_baselines_compare.csv")
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote", out)

    # ---------------- results/matrix_curve_surface.csv ----------------
    out2 = os.path.join(RES, "matrix_curve_surface.csv")
    with open(out2, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["case", "n_real_roots_in_[0,1]", "n_distinct_roots",
                    "n_eigenvalues_total", "n_roots_truth", "max_residual",
                    "time_ms", "roots"])
        for r in mcs:
            w.writerow(r)
    print("wrote", out2)

    # ---------------- results/matrix_certified_1d.csv ----------------
    out3 = os.path.join(RES, "matrix_certified_1d.csv")
    with open(out3, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["case", "cover_boxes", "excluded_boxes",
                    "certified_root_clusters", "time_ms"])
        for r in mcert:
            w.writerow(r)
    print("wrote", out3)

    # ---------------- results/matrix_surface_surface.csv ----------------
    out4 = os.path.join(RES, "matrix_surface_surface.csv")
    with open(out4, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["fixture", "n_t_samples", "real_roots_total",
                    "curve_points", "complex_pts", "max_residual", "time_ms",
                    "roots_per_t_min", "roots_per_t_max"])
        for r in mss:
            w.writerow(r)
    print("wrote", out4)


if __name__ == "__main__":
    main()
