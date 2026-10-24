import gzip, numpy as np, pandas as pd
from sklearn.metrics import roc_curve
W = "/root/mena"; ann = pd.read_csv(f"{W}/results/gse73002_annotation.tsv", sep="\t", index_col=0)
rows, hdr = [], None
with gzip.open(f"{W}/GSE73002_series_matrix.txt.gz", "rt", errors="ignore") as fh:
    on = False
    for l in fh:
        if l.startswith("!series_matrix_table_begin"): on = True; continue
        if l.startswith("!series_matrix_table_end"): break
        if on: rows.append(l.rstrip("\n").split("\t"))
ex = pd.DataFrame([r[1:] for r in rows[1:]], index=[r[0].strip('"') for r in rows[1:]],
                  columns=[x.strip('"') for x in rows[0][1:]]).apply(pd.to_numeric, errors="coerce")
mp = {}
for l in open(f"{W}/ref/mature.fa"):
    if l.startswith(">hsa-"): p = l[1:].split(); mp[p[1]] = p[0]
ex.index = [mp.get(i, i) for i in ex.index]
c = ann["class"].reindex(ex.columns)
L = ["hsa-let-7g-5p", "hsa-miR-138-5p", "hsa-miR-142-3p", "hsa-miR-24-3p", "hsa-miR-32-3p"]
ctl = c == "control"
z = ex.loc[L].sub(ex.loc[L, ctl].mean(1), axis=0).div(ex.loc[L, ctl].std(1), axis=0)
s = z.mean()
out = []
for comp, case, neg in [("BC vs control", c == "breast_cancer", ctl), ("BC vs benign", c == "breast_cancer", c == "benign_breast"),
                        ("Benign vs control", c == "benign_breast", ctl)]:
    m = case | neg
    f, t, _ = roc_curve(case[m].astype(int), s[m])
    idx = np.unique(np.linspace(0, len(f) - 1, min(len(f), 150)).astype(int))
    for i in idx: out.append(f"{comp}\t{f[i]:.4f}\t{t[i]:.4f}")
for k in ["breast_cancer", "benign_breast", "control"]:
    v = s[c == k].values
    out.append(f"SCORE_{k}\t" + ",".join(f"{x:.3f}" for x in np.quantile(v, [0.05, 0.25, 0.5, 0.75, 0.95])))
open(f"{W}/results/gse_roc_points.tsv", "w").write("\n".join(out) + "\n")
print("roc points", len(out))
