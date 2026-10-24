#!/usr/bin/env python3
"""GSE73002 benign breast disease specificity test of the locked MENA signature."""
import os, re, sys, gzip, urllib.request, warnings
import numpy as np, pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
warnings.filterwarnings("ignore")
W = sys.argv[1] if len(sys.argv) > 1 else "/root/mena"
os.makedirs(f"{W}/results", exist_ok=True)
R = lambda n: f"{W}/results/{n}"
log = open(R("gse73002_summary.txt"), "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); log.write(s + "\n"); log.flush()

f = f"{W}/GSE73002_series_matrix.txt.gz"
if not os.path.exists(f):
    urllib.request.urlretrieve("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73002/matrix/GSE73002_series_matrix.txt.gz", f)
hdr, chars, rows = {}, [], []
with gzip.open(f, "rt", errors="ignore") as fh:
    intable = False
    for line in fh:
        if line.startswith("!series_matrix_table_begin"): intable = True; continue
        if line.startswith("!series_matrix_table_end"): break
        if intable: rows.append(line.rstrip("\n").split("\t")); continue
        if line.startswith("!Sample_"):
            k, *v = line.rstrip("\n").split("\t"); v = [x.strip('"') for x in v]
            if k == "!Sample_characteristics_ch1": chars.append(v)
            else: hdr.setdefault(k, v)
ids = hdr["!Sample_geo_accession"]
ann = pd.DataFrame({"gsm": ids, "title": hdr.get("!Sample_title", [""] * len(ids))})
for i, c in enumerate(chars):
    key = next((x.split(":")[0].strip() for x in c if ":" in x), f"char{i}")
    ann[key] = [x.split(":", 1)[1].strip() if ":" in x else x for x in c]
ann = ann.set_index("gsm")
say("Characteristic fields:", list(ann.columns))
for c in ann.columns:
    if c == "title": continue
    vc = ann[c].value_counts()
    if len(vc) <= 40: say(f"\n[{c}]\n" + vc.to_string())
    else: say(f"\n[{c}] {len(vc)} distinct values, e.g. {list(vc.index[:8])}")

ex = pd.DataFrame([r[1:] for r in rows[1:]], index=[r[0].strip('"') for r in rows[1:]],
                  columns=[x.strip('"') for x in rows[0][1:]]).apply(pd.to_numeric, errors="coerce")
say(f"\nExpression: {ex.shape[0]} probes x {ex.shape[1]} samples; value range {np.nanmin(ex.values):.2f}..{np.nanmax(ex.values):.2f}")
say("Example probe IDs:", list(ex.index[:5]))
# map MIMAT -> name via miRBase mature.fa if needed
if ex.index.str.startswith("MIMAT").mean() > 0.5:
    mp = {}
    for line in open(f"{W}/ref/mature.fa"):
        if line.startswith(">hsa-"):
            p = line[1:].split(); mp[p[1]] = p[0]
    ex.index = [mp.get(i, i) for i in ex.index]
if np.nanmax(ex.values) > 40: ex = np.log2(ex.clip(lower=1))

# ---- classify ----
txt = ann.drop(columns=["title"]).astype(str).agg(" | ".join, axis=1).str.lower() + " | " + ann.title.str.lower()
def cls(t):
    if "benign" in t and "breast" in t: return "benign_breast"
    if "breast cancer" in t and "benign" not in t: return "breast_cancer"
    if re.search(r"non-cancer|non cancer|healthy|no cancer|control", t) and "cancer patient" not in t: return "control"
    return "other"
ann["class"] = txt.map(cls)
say("\nClass counts:\n" + ann["class"].value_counts().to_string())
say("\nRaw text examples per class:")
for k in ann["class"].unique(): say(k, "->", txt[ann["class"] == k].iloc[0][:200])
stage_col = next((c for c in ann.columns if re.search(r"stage|^t$|pt|tnm", c, re.I)), None)
say("Stage column:", stage_col)
if stage_col: say(ann.loc[ann["class"] == "breast_cancer", stage_col].value_counts().to_string())
ann.to_csv(R("gse73002_annotation.tsv"), sep="\t")

# ---- scores ----
LASSO = ["hsa-let-7g-5p", "hsa-miR-138-5p", "hsa-miR-142-3p", "hsa-miR-24-3p", "hsa-miR-32-3p"]
HIGH3 = ["hsa-let-7g-5p", "hsa-miR-24-3p", "hsa-miR-142-3p"]
SHIM = ["hsa-miR-1246", "hsa-miR-1307-3p", "hsa-miR-4634", "hsa-miR-6861-5p", "hsa-miR-6875-5p"]
for nm, L in [("LASSO5", LASSO), ("SHIM5", SHIM)]:
    say(f"{nm} present on array: {[m for m in L if m in ex.index]}  missing: {[m for m in L if m not in ex.index]}")
ex = ex[ann.index.intersection(ex.columns)]
c = ann["class"].reindex(ex.columns)

def delong(y, s):
    y = np.asarray(y).astype(bool); s = np.asarray(s, float)
    pos, neg = s[y], s[~y]; m, n = len(pos), len(neg)
    V10 = np.array([((p > neg).sum() + 0.5 * (p == neg).sum()) / n for p in pos])
    V01 = np.array([((pos > q).sum() + 0.5 * (pos == q).sum()) / m for q in neg])
    auc = V10.mean(); se = np.sqrt(V10.var(ddof=1) / m + V01.var(ddof=1) / n)
    return auc, max(0, auc - 1.96 * se), min(1, auc + 1.96 * se)

def spec_at_sens(s_case, s_neg, sens=0.90, B=2000, rng=np.random.default_rng(7)):
    thr = np.quantile(s_case, 1 - sens); sp = (s_neg < thr).mean()
    bs = []
    for _ in range(B):
        a = rng.choice(s_case, len(s_case)); b = rng.choice(s_neg, len(s_neg))
        t = np.quantile(a, 1 - sens); bs.append((b < t).mean())
    return sp, np.quantile(bs, .025), np.quantile(bs, .975), thr

zx = ex.sub(ex.loc[:, c == "control"].mean(1), axis=0).div(ex.loc[:, c == "control"].std(1).replace(0, np.nan), axis=0)
scores = {}
for nm, L in [("LASSO5_direction", LASSO), ("HIGH3_direction", HIGH3)]:
    L = [m for m in L if m in zx.index]; scores[nm] = zx.loc[L].mean()
train = c.isin(["breast_cancer", "control"])  # benign held out
for nm, L in [("LASSO5_refit", LASSO), ("HIGH3_refit", HIGH3), ("SHIM5_refit_posctrl", SHIM)]:
    L = [m for m in L if m in ex.index]
    if not L: continue
    X = ex.loc[L].T.fillna(ex.loc[L].T.mean())
    lr = LogisticRegression(max_iter=2000).fit(X[train], (c[train] == "breast_cancer").astype(int))
    scores[nm] = pd.Series(lr.decision_function(X), index=X.index)

early = None
if stage_col:
    st = ann[stage_col].reindex(ex.columns).astype(str).str.lower()
    early = (c == "breast_cancer") & st.str.contains(r"\b(0|1|t0|t1|tis|is|i|ia|ib)\b")
CLIN_PREV, SCREEN_PREV = 0.10, 0.005
res = []
for nm, s in scores.items():
    for comp, case, neg in [("BC vs benign", c == "breast_cancer", c == "benign_breast"),
                            ("Early BC vs benign", early if early is not None else None, c == "benign_breast"),
                            ("Benign vs control", c == "benign_breast", c == "control"),
                            ("BC vs control", c == "breast_cancer", c == "control")]:
        if case is None or case.sum() < 5 or neg.sum() < 5: continue
        m = case | neg
        auc, lo, hi = delong(case[m], s[m])
        sp, slo, shi, thr = spec_at_sens(s[case].values, s[neg].values)
        prev = CLIN_PREV if "benign" in comp.split(" vs ")[1] else SCREEN_PREV
        ppv = 0.9 * prev / (0.9 * prev + (1 - sp) * (1 - prev)) if comp != "Benign vs control" else np.nan
        res.append([nm, comp, int(case.sum()), int(neg.sum()), auc, lo, hi, sp, slo, shi, prev, ppv])
rt = pd.DataFrame(res, columns=["score", "comparison", "n_case", "n_neg", "AUC", "AUC_lo", "AUC_hi",
                                "spec_at_90sens", "spec_lo", "spec_hi", "prevalence_assumed", "PPV_at_90sens"])
rt.to_csv(R("gse73002_results.tsv"), sep="\t", index=False)
say("\n== GSE73002 results =="); say(rt.round(3).to_string(index=False))
say("DONE")
