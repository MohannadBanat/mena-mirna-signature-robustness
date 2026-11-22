#!/usr/bin/env python3
"""Plan v5, R4: library-quality results by group (descriptive). Per-metric AUCs (cases vs controls, DeLong) and
direction, from the main-analysis QC table; classifier AUCs are copied from the existing output, not recomputed."""
import numpy as np, pandas as pd
W = "/home/claude/rev"; OUT = f"{W}/revision_2026-11-21/out"
q = pd.read_csv(f"{W}/main/results/qc_full.tsv", sep="\t", index_col=0).sort_index()
y = np.array([s.startswith("BC") for s in q.index])
def delong(y, s):
    y = np.asarray(y).astype(bool); s = np.asarray(s, float); pos, neg = s[y], s[~y]; m, n = len(pos), len(neg)
    V10 = np.array([((p > neg).sum() + 0.5 * (p == neg).sum()) / n for p in pos]); V01 = np.array([((pos > q_).sum() + 0.5 * (pos == q_).sum()) / m for q_ in neg])
    auc = V10.mean(); se = np.sqrt(V10.var(ddof=1) / m + V01.var(ddof=1) / n); return auc, max(0, auc - 1.96 * se), min(1, auc + 1.96 * se)
metrics = {
    "Raw reads": q.reads,
    "Fraction without adapter": q.no_adapter / q.reads,
    "Fraction with insert outside 16-30 nt": q.bad_len / q.reads,
    "Retained inserts": q.kept,
    "Alignment rate": q.reads_aligned / q.kept,
    "Duplication rate": q.dup_rate,
    "UMI molecules": q.umi_total,
    "miRNAs detected": q.mirnas_detected,
}
rows = []
for k, v in metrics.items():
    a, lo, hi = delong(y, v.values)
    rows.append([k, "higher in cases" if a >= 0.5 else "lower in cases", round(max(a, 1 - a), 3),
                 round(lo if a >= 0.5 else 1 - hi, 3), round(hi if a >= 0.5 else 1 - lo, 3)])
t = pd.DataFrame(rows, columns=["metric", "direction", "AUC (either direction)", "lo", "hi"])
t.to_csv(f"{OUT}/r4_metric_aucs.tsv", sep="\t", index=False); print(t.to_string(index=False))
