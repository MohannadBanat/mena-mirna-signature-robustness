#!/usr/bin/env python3
"""Plan v5, R2: normalisation sensitivity of the panel scores (152 cases vs 48 controls)."""
import json, warnings
import numpy as np, pandas as pd
from scipy.stats import rankdata
from rnanorm import TMM
warnings.filterwarnings("ignore")
W = "/home/claude/rev"; OUT = f"{W}/revision_2026-11-21/out"
sig = pd.read_csv(f"{W}/main/results/locked_signature_v1.csv").set_index("miRNA ID")
UP = sig.index[sig.direction == "up"].tolist(); DOWN = sig.index[sig.direction == "down"].tolist()
LASSO = sig.index[sig.LASSO5].tolist(); HIGH3 = ["hsa-let-7g-5p", "hsa-miR-24-3p", "hsa-miR-142-3p"]
PANELS = [("Up8", UP, []), ("Down11", [], DOWN), ("LASSO5", LASSO, []), ("High3", HIGH3, [])]

def delong(y, s):
    y = np.asarray(y).astype(bool); s = np.asarray(s, float); pos, neg = s[y], s[~y]; m, n = len(pos), len(neg)
    V10 = np.array([((p > neg).sum() + 0.5 * (p == neg).sum()) / n for p in pos]); V01 = np.array([((pos > q).sum() + 0.5 * (pos == q).sum()) / m for q in neg])
    auc = V10.mean(); se = np.sqrt(V10.var(ddof=1) / m + V01.var(ddof=1) / n); return auc, max(0, auc - 1.96 * se), min(1, auc + 1.96 * se)

def cpm(M): return M / M.sum() * 1e6
def mor(M):
    pos = M[(M > 0).all(1)]; lg = np.log(pos); gm = lg.mean(1)
    sf = np.exp(lg.sub(gm, axis=0).median()); return M / sf
def tmm(M):
    f = TMM().fit(M.T.values)                        # rnanorm: samples x genes
    nf = pd.Series(f.get_norm_factors(M.T.values), index=M.columns)
    return M / (M.sum() * nf) * 1e6
def uq(M):
    E = M[(M > 0).any(axis=1)]; q = E.quantile(0.75)
    return M / q * q.mean()
def quantile_log(M):
    L = np.log2(cpm(M) + 1); ref = np.sort(L.values, axis=0).mean(1); n = len(ref)
    out = np.empty_like(L.values)
    for j in range(L.shape[1]):
        r = rankdata(L.values[:, j], method="average") - 1
        out[:, j] = np.interp(r, np.arange(n), ref)
    return pd.DataFrame(out, index=L.index, columns=L.columns)
NORMS = [("CPM", lambda M: np.log2(cpm(M) + 1)), ("median-of-ratios", lambda M: np.log2(mor(M) + 1)),
         ("TMM", lambda M: np.log2(tmm(M) + 1)), ("upper-quartile", lambda M: np.log2(uq(M) + 1)), ("quantile", quantile_log)]

def run(tag, files):
    rows = []
    for unit, path in files:
        M = pd.read_csv(path, sep="\t", index_col=0); S = sorted(M.columns); M = M[S].astype(float)
        y = np.array([s.startswith("BC") for s in S])
        for nname, f in NORMS:
            lc = f(M); z = lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, np.nan), axis=0)
            for pname, up, dn in PANELS:
                u = [m for m in up if m in z.index]; d = [m for m in dn if m in z.index]
                nd = int(z.loc[u + d].notna().all(1).sum())
                sz = (z.loc[u].mean() if u else 0) - (z.loc[d].mean() if d else 0)
                sm = (lc.loc[u].mean() if u else 0) - (lc.loc[d].mean() if d else 0)
                for sname, s in [("z-score", sz), ("mean log", sm)]:
                    a, lo, hi = delong(y, s.values)
                    rows.append([tag, unit, nname, pname, sname, nd, round(a, 3), round(lo, 3), round(hi, 3)])
    return pd.DataFrame(rows, columns=["matrices", "unit", "normalisation", "panel", "score", "n_scored", "AUC", "lo", "hi"])

prim = run("main", [("UMI", f"{W}/main/results/counts_matrix_umi_mirbase22.tsv"), ("reads", f"{W}/main/results/counts_matrix_reads_mirbase22.tsv")])
sec = run("reproduction", [("UMI", f"{W}/rerun/counts_matrix_umi.tsv"), ("reads", f"{W}/rerun/counts_matrix_reads.tsv")])
allr = pd.concat([prim, sec]); allr.to_csv(f"{OUT}/r2_normalisation_grid.tsv", sep="\t", index=False)

chk = prim[(prim.unit == "UMI") & (prim.normalisation == "CPM") & (prim.score == "z-score")].set_index("panel").AUC
print("Reproduction check vs manuscript (main UMI, CPM, z-score):", chk.to_dict(), "expected Up8 0.925, Down11 0.660, LASSO5 0.922, High3 0.890")
d = prim[prim.panel == "Down11"]
mx = d.AUC.max(); mn = d.AUC.min()
verdict = "pipeline-sensitive" if mx >= 0.75 else ("robust to normalisation" if mx < 0.70 else "partial")
res = {"down_primary_min": mn, "down_primary_max": mx, "verdict": verdict,
       "down_reproduction_min": sec[sec.panel == "Down11"].AUC.min(), "down_reproduction_max": sec[sec.panel == "Down11"].AUC.max()}
print(json.dumps(res)); json.dump(res, open(f"{OUT}/r2_decision.json", "w"))
print(allr.pivot_table(index=["matrices", "unit", "normalisation", "score"], columns="panel", values="AUC").round(3).to_string())
