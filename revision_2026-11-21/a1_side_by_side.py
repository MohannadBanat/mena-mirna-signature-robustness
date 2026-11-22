#!/usr/bin/env python3
"""Amendment v5a, A1: CPM vs DESeq2 median-of-ratios, side by side (descriptive)."""
import json, warnings
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
W = "/home/claude/rev"; OUT = f"{W}/revision_2026-11-21/out"; R = f"{W}/repo/data"
M = pd.read_csv(f"{W}/main/results/counts_matrix_umi_mirbase22.tsv", sep="\t", index_col=0); S = sorted(M.columns); M = M[S].astype(float)
y = np.array([s.startswith("BC") for s in S]).astype(int)
sig = pd.read_csv(f"{W}/main/results/locked_signature_v1.csv").set_index("miRNA ID"); L5 = sig.index[sig.LASSO5].tolist()
pub = {r.split(",")[0]: float(r.split(",")[1]) for r in open(f"{R}/scores.txt").read().strip().split(";")}
qc = {r.split(",")[0]: [float(x) for x in r.split(",")[1:]] for r in open(f"{R}/qc_full.txt").read().strip().split(";")}
g = lambda j: np.array([qc[s][j] for s in S]); reads, kept, noad, badl, ali, umi, det = [g(j) for j in range(7)]
NA = noad / reads * 100
Q = np.c_[np.log(reads), noad / reads, badl / reads, np.log(kept), ali / kept, 1 - umi / ali, np.log(umi), det]

def mor(M):
    pos = M[(M > 0).all(1)]; lg = np.log(pos); return M / np.exp(lg.sub(lg.mean(1), axis=0).median())
def score(lc):
    z = lc.sub(lc.mean(1), axis=0).div(lc.std(1), axis=0); return z.loc[L5].mean().values
def boot(yy, ss, n=2000):
    rng = np.random.default_rng(0); a = []
    for _ in range(n):
        i = rng.integers(0, len(yy), len(yy))
        if yy[i].min() == yy[i].max(): continue
        a.append(roc_auc_score(yy[i], ss[i]))
    return [round(roc_auc_score(yy, ss), 3), *np.round(np.percentile(a, [2.5, 97.5]), 3)]
def cv(X, seed):
    m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000, class_weight="balanced"))
    return roc_auc_score(y, cross_val_predict(m, X, y, cv=StratifiedKFold(10, shuffle=True, random_state=seed), method="predict_proba")[:, 1])

res = {}
for name, lc in [("CPM", np.log2(M / M.sum() * 1e6 + 1)), ("median-of-ratios", np.log2(mor(M) + 1))]:
    L = score(lc); r = {}
    if name == "CPM":
        d = max(abs(L[i] - pub[s]) for i, s in enumerate(S)); r["reproduction_max_abs_diff"] = round(float(d), 4)
        if d > 0.001: raise SystemExit(f"CPM scores do not reproduce published scores (max diff {d})")
    r["AUC_all"] = boot(y, L)
    lo, hi = NA[y == 0].min(), NA[y == 0].max(); m = (NA >= lo) & (NA <= hi)
    r["matched"] = {"cases": int((m & (y == 1)).sum()), "controls": int((m & (y == 0)).sum()), "AUC": boot(y[m], L[m])}
    for nm, X in [("LASSO5 alone", L[:, None]), ("quality + LASSO5", np.c_[Q, L])]:
        a = np.array([cv(X, s) for s in range(20)]); r[nm] = [round(a.mean(), 3), round(a.min(), 3), round(a.max(), 3)]
    cl = np.array([s in [f"BC{k}" for k in range(144, 153)] for s in S]); other = (y == 1) & ~cl
    cm = L[cl].mean(); rng = np.random.default_rng(0); cases = np.where(y == 1)[0]
    draws = np.array([L[rng.choice(cases, 9, replace=False)].mean() for _ in range(100000)])
    keep = ~cl; rho = spearmanr(L[y == 1], umi[y == 1])
    r["cluster"] = {"mean_BC144_152": round(cm, 2), "all_below_zero": bool((L[cl] < 0).all()), "mean_other_cases": round(L[other].mean(), 2),
                    "mean_controls": round(L[y == 0].mean(), 2), "pct_cases_below0": round((L[y == 1] < 0).mean() * 100), "pct_controls_below0": round((L[y == 0] < 0).mean() * 100),
                    "draws_le_cluster_mean_of_100000": int((draws <= cm).sum()), "AUC_excl_cluster": boot(y[keep], L[keep]),
                    "spearman_score_umi_cases": [round(rho.statistic, 2), float(f"{rho.pvalue:.2g}")]}
    res[name] = r
    pd.Series(L, index=S).to_csv(f"{OUT}/a1_lasso5_scores_{name}.tsv", sep="\t", header=["LASSO5"])
json.dump(res, open(f"{OUT}/a1_side_by_side.json", "w"), indent=1); print(json.dumps(res, indent=1))
