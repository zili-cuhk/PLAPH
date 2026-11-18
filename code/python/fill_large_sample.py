"""Rows of Supplementary Tables tab:pl-limits and tab:baseline-sieve from their result files.

Usage: python3 fill_large_sample.py RESULTS_SIMULATION_DIR [SUPP_TEX]
Prints the table rows; with SUPP_TEX, also compares them with the rows in that file and reports
any difference.
"""
import csv, os, re, sys


def num(x, d):
    s = f"{float(x):.{d}f}"
    return s.replace("-", "$-$") if s.startswith("-") else s


def pl_rows(path):
    out = []
    for r in csv.DictReader(open(path)):
        lab = r["design"].split()[0]
        ratio = (float(r["ratio_Z1"]) + float(r["ratio_Z3"]) + float(r["ratio_Z6"])) / 3
        z = [max(float(r[f"z_{m}_n{n}"]) for m in ("Z2", "Z4", "Z7")) for n in (500, 700)]
        out.append(" & ".join([lab, f"{float(r['tau']):.1f}", f"{float(r['censoring_rate']):.2f}",
                               num(r["limit_Z2"], 3), num(r["limit_Z4"], 3), num(r["limit_Z7"], 3),
                               f"{ratio:.3f}", f"{z[0]:.1f}", f"{z[1]:.1f}"]) + r" \\")
    return out


def sieve_rows(path):
    return [" & ".join([r["design"].split()[0]] + [f"{float(r[k]):.3f}" for k in
            ("ratio_k10", "ratio_k20", "ratio_k40", "ratio_partial_likelihood")]) + r" \\"
            for r in csv.DictReader(open(path))]


def body(tex, label):
    a = tex.index(r"\label{" + label + "}")
    m = tex.index(r"\midrule", a) + len(r"\midrule"); b = tex.index(r"\bottomrule", a)
    return [L.strip() for L in tex[m:b].strip().split("\n") if L.strip()]


if __name__ == "__main__":
    res = sys.argv[1]
    tabs = {"tab:pl-limits": pl_rows(os.path.join(res, "pl_limits", "pl_limits.csv")),
            "tab:baseline-sieve": sieve_rows(os.path.join(res, "baseline_sieve", "baseline_sieve_limits.csv"))}
    for lab, rows in tabs.items():
        print(f"% {lab}"); print("\n".join(rows))
    if len(sys.argv) > 2:
        tex = open(sys.argv[2]).read(); bad = 0
        for lab, rows in tabs.items():
            have = body(tex, lab)
            for i, (x, y) in enumerate(zip(rows, have)):
                if re.sub(r"\s+", "", x) != re.sub(r"\s+", "", y):
                    bad += 1; print(f"DIFF {lab} row {i+1}:\n  file: {y}\n  csv:  {x}")
            if len(rows) != len(have):
                bad += 1; print(f"DIFF {lab}: {len(have)} rows in file, {len(rows)} from csv")
        print(f"{bad} differences")
