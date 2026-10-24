#!/usr/bin/env python3
"""Plan v4 + amendment v4a: (A) read counts vs UMI molecules; (B) nested cross-validated LASSO; run after analysis.py."""
import os, sys, time, warnings
import numpy as np, pandas as pd
from collections import Counter
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
W = sys.argv[1] if len(sys.argv) > 1 else "/root/mena"
SEEDS = int(os.environ.get("SEEDS", 10))
R = lambda n: f"{W}/results/{n}"
os.makedirs(f"{W}/results", exist_ok=True)
log = open(R("extra_summary.txt"), "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); log.write(s + "\n"); log.flush()

umi = pd.read_csv(f"{W}/counts_matrix_umi.tsv", sep="\t", index_col=0)
rds = pd.read_csv(f"{W}/counts_matrix_reads.tsv", sep="\t", index_col=0)
sig = pd.read_csv(f"{W}/locked_signature_v1.csv").set_index("miRNA ID")
pub = pd.read_csv(f"{W}/pub_s1.tsv", sep="\t").set_index("mirna")
if os.path.exists(f"{W}/ref/mature22.fa"):
    def fa_map(p):
        m = {}
        for l in open(p):
            if l.startswith(">hsa-"): x = l[1:].split(); m[x[0]] = x[1]
        return m
    new_acc = fa_map(f"{W}/ref/mature.fa"); old = {v: k for k, v in fa_map(f"{W}/ref/mature22.fa").items()}
    ren = {n: old[a] for n, a in new_acc.items() if a in old and old[a] != n}
    umi = umi.rename(index=ren).groupby(level=0).sum(); rds = rds.rename(index=ren).groupby(level=0).sum()
samples = sorted(umi.columns); umi = umi[samples]; rds = rds.reindex(index=umi.index, columns=samples).fillna(0)
grp = pd.Series(["BC" if s.startswith("BC") else "HC" for s in samples], index=samples); y = (grp == "BC").values
say(f"Rerun matrix: {umi.shape[0]} miRNAs x {umi.shape[1]} libraries; UMI molecules {int(umi.values.sum()):,}; aligned reads {int(rds.values.sum()):,}")
LASSO = sig.index[sig.LASSO5].tolist(); UP = sig.index[sig.direction == "up"].tolist(); DOWN = sig.index[sig.direction == "down"].tolist()
HIGH3 = ["hsa-let-7g-5p", "hsa-miR-24-3p", "hsa-miR-142-3p"]; LOW2 = ["hsa-miR-138-5p", "hsa-miR-32-3p"]

def delong(yy, s):
    yy = np.asarray(yy).astype(bool); s = np.asarray(s, float); pos, neg = s[yy], s[~yy]; m, n = len(pos), len(neg)
    V10 = np.array([((p > neg).sum() + 0.5 * (p == neg).sum()) / n for p in pos]); V01 = np.array([((pos > q).sum() + 0.5 * (pos == q).sum()) / m for q in neg])
    auc = V10.mean(); se = np.sqrt(V10.var(ddof=1) / m + V01.var(ddof=1) / n); return auc, max(0, auc - 1.96 * se), min(1, auc + 1.96 * se)
def lcpm(M): return np.log2(M / M.sum() * 1e6 + 1)
def zsc(lc): return lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, np.nan), axis=0)

# ---------- A1/A2/A5 panel scores, reads vs UMI ----------
say("\n== A1/A2/A5: panel AUCs, 152 vs 48, identical scoring on UMI molecules vs aligned reads ==")
rows = []
for unit, M in [("UMI", umi), ("reads", rds)]:
    lc = lcpm(M); z = zsc(lc)
    for name, up, dn in [("Up8", UP, []), ("Down11", [], DOWN), ("LASSO5", LASSO, []), ("High3", HIGH3, []), ("Low-count pair", LOW2, [])]:
        u = [m for m in up if m in z.index]; d = [m for m in dn if m in z.index]
        sz = (z.loc[u].mean() if u else 0) - (z.loc[d].mean() if d else 0)            # z-scored (this study)
        sm = (lc.loc[u].mean() if u else 0) - (lc.loc[d].mean() if d else 0)          # original style: mean log2 CPM
        nd = int(z.loc[u + d].notna().all(1).sum())
        a, lo, hi = delong(y, sz.values); a2, lo2, hi2 = delong(y, sm.values)
        rows.append([unit, name, nd, round(a, 3), round(lo, 3), round(hi, 3), round(a2, 3), round(lo2, 3), round(hi2, 3)])
pa = pd.DataFrame(rows, columns=["unit", "panel", "n_scored", "AUC_z", "lo", "hi", "AUC_meanlog", "lo2", "hi2"])
pa.to_csv(R("extra_panel_auc_reads_vs_umi.tsv"), sep="\t", index=False); say(pa.to_string(index=False))
dr = pa[(pa.unit == "reads") & (pa.panel == "Down11")].AUC_z.iloc[0]; du = pa[(pa.unit == "UMI") & (pa.panel == "Down11")].AUC_z.iloc[0]
verdict = "SUPPORTED" if dr >= 0.75 else ("NOT SUPPORTED" if dr < 0.70 else "PARTIAL")
say(f"A1 decision (plan v4): down-panel AUC reads {dr:.3f} vs UMI {du:.3f} -> hypothesis {verdict}")

# ---------- A3 DE with reads vs UMI ----------
say("\n== A3: female-only DESeq2, reads vs UMI ==")
de = {u: pd.read_csv(R(f"de_female_only_{u}.tsv"), sep="\t", index_col=0) for u in ["umi", "reads"]}
for u, d in de.items():
    t = [m for m in sig.index if m in d.index]
    rep = [m for m in t if d.loc[m, "padj"] < 0.05 and np.sign(d.loc[m, "log2FoldChange"]) == (1 if sig.loc[m, "direction"] == "up" else -1)]
    sh = d.join(pub[["pub_lfc", "pub_baseMean"]], how="inner", rsuffix="_p").dropna(subset=["log2FoldChange", "pub_lfc"])
    r_all = stats.spearmanr(sh.log2FoldChange, sh.pub_lfc); lowc = sh[sh.pub_baseMean < 10]
    r_low = stats.spearmanr(lowc.log2FoldChange, lowc.pub_lfc) if len(lowc) > 5 else (np.nan, np.nan)
    pubsig = pub[(pub.pub_padj < 0.05) & (pub.pub_lfc.abs() >= 1)]
    tested = pubsig.index[pubsig.index.isin(d.index)]
    repl = sum((d.loc[m, "padj"] < 0.05) and np.sign(d.loc[m, "log2FoldChange"]) == np.sign(pubsig.loc[m, "pub_lfc"]) for m in tested)
    say(f"{u}: miRNAs tested {len(d)}; panel tested {len(t)}/19, replicated {len(rep)} {sorted(rep)}; published DE tested {len(tested)}, replicated {repl}; "
        f"rho vs published log2FC all {r_all[0]:.3f} (n={len(sh)}), published baseMean<10 {r_low[0]:.3f} (n={len(lowc)})")

# ---------- A4 reads per molecule ----------
say("\n== A4: reads per molecule (PCR copies) ==")
ratio = (rds / umi.replace(0, np.nan))
mu = umi.mean(1)
for lab, sel in [("mean <10 molecules", mu < 10), ("10-100", (mu >= 10) & (mu < 100)), (">=100", mu >= 100)]:
    rr = ratio.loc[sel]; med = rr.median(1); cv = rr.std(1) / rr.mean(1)
    say(f"{lab}: n={int(sel.sum())}, median reads/molecule {np.nanmedian(med):.2f}, median CV across libraries {np.nanmedian(cv):.2f}")
lib_bc = (rds.sum() / umi.sum())[grp == "BC"]; lib_hc = (rds.sum() / umi.sum())[grp == "HC"]
say(f"library reads/molecule: cases median {lib_bc.median():.2f}, controls {lib_hc.median():.2f}, MWU p={stats.mannwhitneyu(lib_bc, lib_hc).pvalue:.3g}")
for m in [x for x in sig.index if x in umi.index and mu.get(x, 0) < 20]:
    a, b = ratio.loc[m, grp == "BC"].dropna(), ratio.loc[m, grp == "HC"].dropna()
    p = stats.mannwhitneyu(a, b).pvalue if len(a) > 2 and len(b) > 2 else np.nan
    say(f"  {m}: reads/molecule cases {a.median() if len(a) else np.nan:.2f} (n={len(a)}), controls {b.median() if len(b) else np.nan:.2f} (n={len(b)}), p={p:.3g}; "
        f"mean reads cases {rds.loc[m, grp=='BC'].mean():.1f} vs controls {rds.loc[m, grp=='HC'].mean():.1f}; mean molecules {umi.loc[m, grp=='BC'].mean():.1f} vs {umi.loc[m, grp=='HC'].mean():.1f}")

# ---------- B nested CV LASSO ----------
say("\n== B: nested cross-validated L1-logistic re-derivation (outer 10-fold x %d seeds, inner 5-fold) ==" % SEEDS)
keep = ((umi >= 10).mean(1) >= 0.2); feats = umi.index[keep]
Xu = lcpm(umi).loc[feats].T.values; Xr = lcpm(rds).loc[feats].T.values
qf = f"{W}/results/qc_full.tsv"
if os.path.exists(qf):
    q = pd.read_csv(qf, sep="\t", index_col=0).reindex(samples)
    Q = np.c_[np.log(q.reads), q.no_adapter / q.reads, q.bad_len / q.reads, np.log(q.kept), q.reads_aligned / q.kept, q.dup_rate, np.log(q.umi_total), q.mirnas_detected]
else: Q = None
Cs = np.logspace(-2, 1, 12)
def nested(X, names):
    aucs, fold_aucs, selc, nsel = [], [], Counter(), []
    for sd in range(SEEDS):
        p = np.zeros(len(y))
        for tr, te in StratifiedKFold(10, shuffle=True, random_state=sd).split(X, y):
            pipe = make_pipeline(StandardScaler(), LogisticRegression(penalty="l1", solver="liblinear", class_weight="balanced", max_iter=5000))
            gs = GridSearchCV(pipe, {"logisticregression__C": Cs}, cv=StratifiedKFold(5, shuffle=True, random_state=100 + sd), scoring="roc_auc", n_jobs=-1)
            gs.fit(X[tr], y[tr]); p[te] = gs.predict_proba(X[te])[:, 1]
            if y[te].min() != y[te].max(): fold_aucs.append(roc_auc_score(y[te], p[te]))
            coef = gs.best_estimator_[-1].coef_[0]; s = [names[i] for i in np.flatnonzero(coef)]; selc.update(s); nsel.append(len(s))
        aucs.append(roc_auc_score(y, p))
    nf = 10 * SEEDS
    return np.array(aucs), np.array(fold_aucs), {k: v / nf for k, v in selc.items()}, np.array(nsel)
variants = [("B1 miRNA, UMI molecules", Xu, list(feats)), ("B2 miRNA, aligned reads", Xr, list(feats))]
if Q is not None: variants.append(("B3 miRNA (UMI) + 8 library-quality metrics", np.c_[Xu, Q], list(feats) + [f"QC{i}" for i in range(8)]))
for nm, X, names in variants:
    t0 = time.time(); a, fa, sf, ns = nested(X, names)
    say(f"{nm}: features {X.shape[1]}; pooled outer AUC mean {a.mean():.3f} (range {a.min():.3f}-{a.max():.3f}); mean per-fold AUC {fa.mean():.3f}; "
        f"median features selected {np.median(ns):.0f}; {time.time()-t0:.0f}s")
    say("   LASSO5 selection frequency: " + ", ".join(f"{m.replace('hsa-','')} {sf.get(m, 0):.2f}" for m in LASSO))
    top = sorted(sf.items(), key=lambda kv: -kv[1])[:12]; say("   most often selected: " + ", ".join(f"{k.replace('hsa-','')} {v:.2f}" for k, v in top))
say("\nEXTRA DONE", time.strftime("%H:%M:%S"))
