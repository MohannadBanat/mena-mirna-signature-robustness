#!/usr/bin/env python3
"""Plan v5, R1: permutation test of the nested cross-validated L1-logistic model (plan v4, B1).
Usage: r1_permutation.py <stage> where stage = reproduce | observe | permute | matched
"""
import sys, json, time, warnings
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")

W = "/home/claude/rev"
OUT = f"{W}/revision_2026-11-21/out"
umi = pd.read_csv(f"{W}/rerun/counts_matrix_umi.tsv", sep="\t", index_col=0)
samples = sorted(umi.columns); umi = umi[samples]
y = np.array([s.startswith("BC") for s in samples])
keep = ((umi >= 10).mean(1) >= 0.2); feats = umi.index[keep]
lcpm = np.log2(umi / umi.sum() * 1e6 + 1)
X = lcpm.loc[feats].T.values
Cs = np.logspace(-2, 1, 12)

def nested_auc(X, y, seed):
    p = np.zeros(len(y))
    for tr, te in StratifiedKFold(10, shuffle=True, random_state=seed).split(X, y):
        pipe = make_pipeline(StandardScaler(), LogisticRegression(penalty="l1", solver="liblinear", class_weight="balanced", max_iter=5000))
        gs = GridSearchCV(pipe, {"logisticregression__C": Cs}, cv=StratifiedKFold(5, shuffle=True, random_state=100 + seed), scoring="roc_auc", n_jobs=-1)
        gs.fit(X[tr], y[tr]); p[te] = gs.predict_proba(X[te])[:, 1]
    return roc_auc_score(y, p)

stage = sys.argv[1]
import os; os.makedirs(OUT, exist_ok=True)
if stage == "reproduce":
    a = np.array([nested_auc(X, y, s) for s in range(10)])
    r = {"features": int(X.shape[1]), "auc_mean": round(a.mean(), 4), "auc_min": round(a.min(), 4), "auc_max": round(a.max(), 4)}
    print(json.dumps(r)); json.dump(r, open(f"{OUT}/r1_reproduce.json", "w"))
elif stage == "observe":
    r = {"observed_auc_seed0": nested_auc(X, y, 0)}
    print(json.dumps(r)); json.dump(r, open(f"{OUT}/r1_observed.json", "w"))
elif stage == "permute":
    rng = np.random.default_rng(2026); null = []
    t0 = time.time()
    for i in range(200):
        yp = rng.permutation(y); null.append(nested_auc(X, yp, 0))
        if (i + 1) % 20 == 0:
            print(f"{i+1}/200 done, {time.time()-t0:.0f}s, running mean {np.mean(null):.3f}", flush=True)
            np.savetxt(f"{OUT}/r1_null_partial.txt", null)
    np.savetxt(f"{OUT}/r1_null.txt", null)
    obs = json.load(open(f"{OUT}/r1_observed.json"))["observed_auc_seed0"]
    null = np.array(null)
    r = {"observed": round(obs, 4), "n_perm": 200, "null_mean": round(null.mean(), 4), "null_sd": round(null.std(ddof=1), 4),
         "null_p95": round(np.percentile(null, 95), 4), "null_max": round(null.max(), 4),
         "n_null_ge_obs": int((null >= obs).sum()), "p": round((1 + (null >= obs).sum()) / 201, 5)}
    print(json.dumps(r)); json.dump(r, open(f"{OUT}/r1_permutation.json", "w"))
elif stage == "matched":
    q = pd.read_csv(f"{W}/main/results/qc_full.tsv", sep="\t", index_col=0).reindex(samples)
    na = q.no_adapter / q.reads * 100
    lo, hi = na[~y].min(), na[~y].max(); m = ((na >= lo) & (na <= hi)).values
    r = {"range_pct": [round(lo, 2), round(hi, 2)], "cases": int((m & y).sum()), "controls": int((m & ~y).sum()),
         "nested_auc_seed0": round(nested_auc(X[m], y[m], 0), 4)}
    print(json.dumps(r)); json.dump(r, open(f"{OUT}/r1_matched.json", "w"))
