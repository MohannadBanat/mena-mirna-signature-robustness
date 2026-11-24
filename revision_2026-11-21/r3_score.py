#!/usr/bin/env python3
"""Plan v5, R3 (with amendment v5a, A2): score the alignment-sensitivity matrices.
Normalisation (a) CPM; scores (i) mean z-score and (ii) mean log2 CPM; panels Up8, Down11, LASSO5, High3."""
import json, sys, warnings
import numpy as np, pandas as pd
from scipy.stats import spearmanr
warnings.filterwarnings("ignore")
W = "/home/claude/rev"; D = sys.argv[1]; OUT = f"{W}/revision_2026-11-21/out"
sig = pd.read_csv(f"{W}/main/results/locked_signature_v1.csv").set_index("miRNA ID")
UP = sig.index[sig.direction == "up"].tolist(); DOWN = sig.index[sig.direction == "down"].tolist()
LASSO = sig.index[sig.LASSO5].tolist(); HIGH3 = ["hsa-let-7g-5p", "hsa-miR-24-3p", "hsa-miR-142-3p"]
PANELS = [("Up8", UP, []), ("Down11", [], DOWN), ("LASSO5", LASSO, []), ("High3", HIGH3, [])]
def delong(y, s):
    y = np.asarray(y).astype(bool); s = np.asarray(s, float); pos, neg = s[y], s[~y]; m, n = len(pos), len(neg)
    V10 = np.array([((p > neg).sum() + 0.5 * (p == neg).sum()) / n for p in pos]); V01 = np.array([((pos > q).sum() + 0.5 * (pos == q).sum()) / m for q in neg])
    auc = V10.mean(); se = np.sqrt(V10.var(ddof=1) / m + V01.var(ddof=1) / n); return auc, max(0, auc - 1.96 * se), min(1, auc + 1.96 * se)
rows = []
for X in "ABC":
    for unit in ["umi", "reads"]:
        M = pd.read_csv(f"{D}/r3_{X}_{unit}.tsv", sep="\t", index_col=0); S = sorted(M.columns); M = M[S].astype(float)
        assert len(S) == 200, f"{X} {unit}: {len(S)} libraries"
        y = np.array([s.startswith("BC") for s in S])
        lc = np.log2(M / M.sum() * 1e6 + 1); z = lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, np.nan), axis=0)
        for pname, up, dn in PANELS:
            u = [m for m in up if m in z.index]; d = [m for m in dn if m in z.index]
            nd = int(z.loc[u + d].notna().all(1).sum())
            for sname, s in [("z-score", (z.loc[u].mean() if u else 0) - (z.loc[d].mean() if d else 0)),
                             ("mean log", (lc.loc[u].mean() if u else 0) - (lc.loc[d].mean() if d else 0))]:
                a, lo, hi = delong(y, s.values); rows.append([X, unit, pname, sname, nd, round(a, 3), round(lo, 3), round(hi, 3)])
        rows.append([X, unit, "_total_molecules_or_reads", "", int(M.values.sum()), np.nan, np.nan, np.nan])
g = pd.DataFrame(rows, columns=["setting", "unit", "panel", "score", "n_scored_or_total", "AUC", "lo", "hi"])
g.to_csv(f"{OUT}/r3_alignment_grid.tsv", sep="\t", index=False)
# reproduction checks for setting A
A = g[(g.setting == "A") & (g.unit == "umi") & (g.panel == "Down11") & (g.score == "z-score")].AUC.iloc[0]
Au = pd.read_csv(f"{D}/r3_A_umi.tsv", sep="\t", index_col=0); Ru = pd.read_csv(f"{W}/rerun/counts_matrix_umi.tsv", sep="\t", index_col=0)
common = Au.index.intersection(Ru.index); rhos = [spearmanr(Au.loc[common, s], Ru.loc[common, s]).statistic for s in sorted(Au.columns)]
chk = {"A_down11_umi_z": A, "main_manuscript": 0.660, "within_0.02_of_main": bool(abs(A - 0.660) <= 0.02), "rerun_2026-10-24": 0.633,
       "A_vs_rerun_spearman_median": round(float(np.median(rhos)), 4), "A_vs_rerun_spearman_min": round(float(np.min(rhos)), 4), "miRNAs_common": len(common)}
dn = g[(g.setting.isin(["B", "C"])) & (g.panel == "Down11")].AUC
chk["down11_BC_min"], chk["down11_BC_max"] = float(dn.min()), float(dn.max()); chk["A_down11_umi_z"] = float(chk["A_down11_umi_z"])
chk["verdict"] = "pipeline-sensitive" if dn.max() >= 0.75 else ("robust to alignment" if dn.max() < 0.70 else "partial")
json.dump(chk, open(f"{OUT}/r3_decision.json", "w"), indent=1); print(json.dumps(chk, indent=1))
print(g[g.panel != "_total_molecules_or_reads"].pivot_table(index=["setting", "unit", "score"], columns="panel", values="AUC").round(3).to_string())
