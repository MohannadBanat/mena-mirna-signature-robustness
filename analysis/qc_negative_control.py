"""Exploratory technical-only negative control (added 24 October 2026, after the main results were known).
Classifiers on per-library quality metrics only (no miRNA expression) vs the locked LASSO5 score."""
import numpy as np, warnings
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
D = "/home/claude/ms/data/"
rd = lambda f: {r.split(",")[0]: [float(x) for x in r.split(",")[1:]] for r in open(D + f).read().strip().split(";")}
qc, sc = rd("qc_full.txt"), rd("scores.txt")
S = sorted(qc); y = np.array([s.startswith("BC") for s in S]).astype(int)
g = lambda j: np.array([qc[s][j] for s in S])
reads, kept, noad, badl, ali, umi, det = [g(j) for j in range(7)]
L = np.array([sc[s][0] for s in S]); NA = noad / reads * 100
Q = np.c_[np.log(reads), noad / reads, badl / reads, np.log(kept), ali / kept, 1 - umi / ali, np.log(umi), det]
def boot(yy, ss, n=2000):
    rng = np.random.default_rng(0); a = []
    for _ in range(n):
        i = rng.integers(0, len(yy), len(yy))
        if yy[i].min() == yy[i].max(): continue
        a.append(roc_auc_score(yy[i], ss[i]))
    return roc_auc_score(yy, ss), *np.percentile(a, [2.5, 97.5])
def cv(X, yy, seed):
    m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000, class_weight="balanced"))
    return roc_auc_score(yy, cross_val_predict(m, X, yy, cv=StratifiedKFold(10, shuffle=True, random_state=seed), method="predict_proba")[:, 1])
rng = np.random.default_rng(1)
for nm, X in [("library quality (8 metrics)", Q), ("quality ratios only (4)", Q[:, [1, 2, 4, 5]]), ("LASSO5 score", L[:, None]), ("library quality + LASSO5", np.c_[Q, L])]:
    a = np.array([cv(X, y, s) for s in range(20)])
    line = f"{nm}: CV AUC mean {a.mean():.3f}, range {a.min():.3f}-{a.max():.3f} (20 x 10-fold)"
    if "LASSO5" not in nm:
        perm = np.array([cv(X, rng.permutation(y), s) for s in range(200)])
        line += f"; permuted labels: mean {perm.mean():.3f}, max {perm.max():.3f}, n>=observed {(perm >= a.mean()).sum()}/200"
    print(line)
print("fraction without adapter, AUC (95% CI):", np.round(boot(y, NA), 3))
print("LASSO5, AUC (95% CI):", np.round(boot(y, L), 3))
for grp, nm in [(1, "cases"), (0, "controls")]:
    r = spearmanr(L[y == grp], NA[y == grp]); print(f"Spearman LASSO5 vs no-adapter, {nm}: rho {r.statistic:.2f}, p {r.pvalue:.3f}")
lo, hi = NA[y == 0].min(), NA[y == 0].max(); m = (NA >= lo) & (NA <= hi)
print(f"control range {lo:.2f}-{hi:.2f}%: cases {(m & (y == 1)).sum()}, controls {(m & (y == 0)).sum()}, LASSO5 AUC", np.round(boot(y[m], L[m]), 3))
