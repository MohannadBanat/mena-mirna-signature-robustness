#!/usr/bin/env python3
"""MENA miRNA stress test: female-only DE, depth matching, hemolysis, panel AUCs.
Run in /root/mena after `pod_pipeline.sh matrix`. Writes results/ ."""
import os, sys, json, time, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

W = sys.argv[1] if len(sys.argv) > 1 else "/root/mena"
N_ITER = int(os.environ.get("N_ITER", 100))
os.makedirs(f"{W}/results", exist_ok=True)
R = lambda n: f"{W}/results/{n}"
log = open(R("summary.txt"), "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); log.write(s + "\n"); log.flush()

# ---------- load ----------
umi = pd.read_csv(f"{W}/counts_matrix_umi.tsv", sep="\t", index_col=0)      # miRNA x sample
rds = pd.read_csv(f"{W}/counts_matrix_reads.tsv", sep="\t", index_col=0)
qc = pd.read_csv(f"{W}/qc_table.tsv", sep="\t").set_index("sample")
pub = pd.read_csv(f"{W}/pub_s1.tsv", sep="\t").set_index("mirna")
pub15 = pd.read_csv(f"{W}/pub_s15.tsv", sep="\t").set_index("sample")
sig = pd.read_csv(f"{W}/locked_signature_v1.csv").set_index("miRNA ID")
# harmonise names to miRBase v22 (used by the original study) via stable MIMAT accessions
def fa_map(p):
    m = {}
    for l in open(p):
        if l.startswith(">hsa-"): x = l[1:].split(); m[x[0]] = x[1]
    return m
new_acc = fa_map(f"{W}/ref/mature.fa"); old_by_acc = {v: k for k, v in fa_map(f"{W}/ref/mature22.fa").items()}
ren = {n: old_by_acc[a] for n, a in new_acc.items() if a in old_by_acc and old_by_acc[a] != n}
umi = umi.rename(index=ren); rds = rds.rename(index=ren)
umi = umi.groupby(level=0).sum(); rds = rds.groupby(level=0).sum()
missing = [m for m in sig.index if m not in umi.index]
samples = sorted(umi.columns)
umi, rds = umi[samples], rds.reindex(columns=samples).fillna(0)
grp = pd.Series(["BC" if s.startswith("BC") else "HC" for s in samples], index=samples)
say(f"Harmonised {len(ren)} miRNA names to published nomenclature; panel miRNAs absent from reference: {missing}")
say(f"Samples: {len(samples)}  BC={sum(grp=='BC')}  HC={sum(grp=='HC')}  miRNAs with any count: {(umi.sum(1)>0).sum()}")
# miRBase name harmonisation: mature.fa uses e.g. hsa-miR-451a ; published uses same
LASSO = sig.index[sig.LASSO5].tolist()
UP = sig.index[sig.direction == "up"].tolist(); DOWN = sig.index[sig.direction == "down"].tolist()
HIGH3 = ["hsa-let-7g-5p", "hsa-miR-24-3p", "hsa-miR-142-3p"]

# ---------- 1. QC & depth ----------
qc["group"] = grp
qc["dup_rate"] = 1 - qc["umi_total"] / qc["sample"].map(rds.sum()) if "sample" in qc else np.nan
qc["reads_aligned"] = rds.sum().reindex(qc.index)
qc["dup_rate"] = 1 - qc["umi_total"] / qc["reads_aligned"]
qc["kept_frac"] = qc["kept"] / qc["reads"]
qc.to_csv(R("qc_full.tsv"), sep="\t")
rows = []
for col in ["reads", "kept_frac", "reads_aligned", "umi_total", "mirnas_detected", "dup_rate"]:
    a, b = qc.loc[qc.group == "BC", col], qc.loc[qc.group == "HC", col]
    p = stats.mannwhitneyu(a, b).pvalue
    rows.append([col, a.median(), a.quantile(.25), a.quantile(.75), b.median(), b.quantile(.25), b.quantile(.75), p])
dep = pd.DataFrame(rows, columns=["metric", "BC_median", "BC_q1", "BC_q3", "HC_median", "HC_q1", "HC_q3", "MWU_p"])
dep.to_csv(R("depth_by_group.tsv"), sep="\t", index=False)
say("\n== Depth / QC by group =="); say(dep.round(4).to_string(index=False))

# ---------- 2. female-only DE ----------
def run_deseq(counts_mxs, g, min_count=10, min_frac=0.2):
    keep = ((counts_mxs >= min_count).mean(1) >= min_frac)
    X = counts_mxs.loc[keep].T.astype(int)
    meta = pd.DataFrame({"group": g.loc[X.index]})
    dds = DeseqDataSet(counts=X, metadata=meta, design="~group", quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["group", "BC", "HC"], quiet=True); st.summary()
    return st.results_df, keep

t0 = time.time()
de, keep = run_deseq(umi, grp)
say(f"\n== Female-only DESeq2 (UMI) == tested {keep.sum()} miRNAs in {time.time()-t0:.0f}s")
de = de.join(pub, how="left")
de.to_csv(R("de_female_only_umi.tsv"), sep="\t")
sigde = de[(de.padj < 0.05) & (de.log2FoldChange.abs() >= 1)]
say(f"DE (padj<0.05, |log2FC|>=1): {len(sigde)}  up={sum(sigde.log2FoldChange>0)} down={sum(sigde.log2FoldChange<0)}")
both = de.dropna(subset=["pub_lfc"])
r, p = stats.spearmanr(both.log2FoldChange, both.pub_lfc)
say(f"Concordance with published Table S1 log2FC (n={len(both)} shared): Spearman rho={r:.3f} p={p:.2g}")
pubsig = pub[(pub.pub_padj < 0.05) & (pub.pub_lfc.abs() >= 1)]
say(f"Published DE miRNAs (padj<.05,|FC|>=2): {len(pubsig)}; of these tested here: {pubsig.index.isin(de.index).sum()}; "
    f"replicated same direction & padj<.05: {sum((de.reindex(pubsig.index).padj<0.05) & (np.sign(de.reindex(pubsig.index).log2FoldChange)==np.sign(pubsig.pub_lfc)))}")

# raw-read DE for UMI comparison
de_r, _ = run_deseq(rds, grp)
de_r.to_csv(R("de_female_only_reads.tsv"), sep="\t")

# ---------- 3. panel miRNA fate + detection ----------
norm = umi / umi.sum() * 1e6  # CPM
pr = []
for m in sig.index:
    inC = m in umi.index
    det_bc = (umi.loc[m, grp == "BC"] > 0).mean() if inC else np.nan
    det_hc = (umi.loc[m, grp == "HC"] > 0).mean() if inC else np.nan
    fisher = stats.fisher_exact([[ (umi.loc[m,grp=="BC"]>0).sum(), (umi.loc[m,grp=="BC"]==0).sum()],
                                 [ (umi.loc[m,grp=="HC"]>0).sum(), (umi.loc[m,grp=="HC"]==0).sum()]])[1] if inC else np.nan
    d = de.loc[m] if m in de.index else None
    dr = de_r.loc[m] if m in de_r.index else None
    pr.append(dict(mirna=m, panel=sig.loc[m, "panel"], LASSO5=bool(sig.loc[m, "LASSO5"]), pub_direction=sig.loc[m, "direction"],
                   pub_lfc=sig.loc[m, "log2FoldChange"], pub_baseMean=sig.loc[m, "baseMean"],
                   mean_umi_BC=umi.loc[m, grp=="BC"].mean() if inC else 0, mean_umi_HC=umi.loc[m, grp=="HC"].mean() if inC else 0,
                   detect_BC=det_bc, detect_HC=det_hc, detect_fisher_p=fisher,
                   tested=d is not None, lfc_umi=None if d is None else d.log2FoldChange, padj_umi=None if d is None else d.padj,
                   lfc_reads=None if dr is None else dr.log2FoldChange, padj_reads=None if dr is None else dr.padj))
pt = pd.DataFrame(pr)
pt["replicated"] = (pt.tested) & (pt.padj_umi < 0.05) & (np.sign(pt.lfc_umi) == np.where(pt.pub_direction == "up", 1, -1))
pt.to_csv(R("panel_mirna_fate.tsv"), sep="\t", index=False)
say("\n== Panel miRNA fate (female-only, UMI) =="); say(pt[["mirna","LASSO5","pub_direction","pub_baseMean","mean_umi_BC","mean_umi_HC","detect_BC","detect_HC","tested","lfc_umi","padj_umi","replicated"]].round(3).to_string(index=False))
say(f"Replicated: {pt.replicated.sum()}/{len(pt)} panel miRNAs; LASSO5 replicated: {pt[pt.LASSO5].replicated.sum()}/5")

# ---------- 4. depth matching (binomial thinning) ----------
tot = umi.sum(); target = int(tot.quantile(0.10))
say(f"\n== Depth matching == target UMI depth (10th pct) = {target:,}; {N_ITER} iterations")
rng = np.random.default_rng(20261019)
rep = {m: 0 for m in sig.index}; tested_n = {m: 0 for m in sig.index}
lfcs = {m: [] for m in sig.index}
detB = {m: [] for m in sig.index}; detH = {m: [] for m in sig.index}
t0 = time.time()
for it in range(N_ITER):
    frac = np.minimum(1.0, target / tot.values)
    thin = pd.DataFrame(rng.binomial(umi.values.astype(np.int64), frac[None, :]), index=umi.index, columns=umi.columns)
    for m in sig.index:
        if m in thin.index:
            detB[m].append((thin.loc[m, grp == "BC"] > 0).mean()); detH[m].append((thin.loc[m, grp == "HC"] > 0).mean())
    try:
        d, _ = run_deseq(thin, grp)
    except Exception as e:
        say("iter fail", it, e); continue
    for m in sig.index:
        if m in d.index:
            tested_n[m] += 1; lfcs[m].append(d.loc[m, "log2FoldChange"])
            want = 1 if sig.loc[m, "direction"] == "up" else -1
            if d.loc[m, "padj"] < 0.05 and np.sign(d.loc[m, "log2FoldChange"]) == want: rep[m] += 1
    if it == 0: say(f"  first iteration {time.time()-t0:.0f}s")
dm = pd.DataFrame({"mirna": list(sig.index), "pub_baseMean": sig.baseMean.values,
                   "iters_tested": [tested_n[m] for m in sig.index],
                   "iters_replicated": [rep[m] for m in sig.index],
                   "median_lfc": [np.median(lfcs[m]) if lfcs[m] else np.nan for m in sig.index]})
dm["detect_BC_thinned"] = [np.mean(detB[m]) if detB[m] else np.nan for m in sig.index]
dm["detect_HC_thinned"] = [np.mean(detH[m]) if detH[m] else np.nan for m in sig.index]
dm["frac_replicated"] = dm.iters_replicated / N_ITER
dm["supported"] = dm.frac_replicated >= 0.5
dm.to_csv(R("depth_matching.tsv"), sep="\t", index=False)
say(dm.round(3).to_string(index=False))

# ---------- 5. hemolysis ----------
lc = np.log2(norm + 1)
def ratio(a, b): return (lc.loc[a] - lc.loc[b]) if a in lc.index and b in lc.index else pd.Series(np.nan, index=lc.columns)
hem = pd.DataFrame({"group": grp, "r451_23a": ratio("hsa-miR-451a", "hsa-miR-23a-3p")})
# alternative reference: abundant, not DE, not erythroid
ery = {"hsa-miR-451a","hsa-miR-16-5p","hsa-miR-486-5p","hsa-miR-92a-3p","hsa-miR-144-3p","hsa-miR-144-5p","hsa-miR-15b-5p","hsa-miR-20a-5p"}
cand = de[(de.baseMean > de.baseMean.quantile(0.75)) & (de.padj > 0.5) & (de.log2FoldChange.abs() < 0.2) & (~de.index.isin(ery))]
ref = cand.sort_values("baseMean", ascending=False).index[0] if len(cand) else None
say(f"\n== Hemolysis == alternative reference miRNA: {ref}")
if ref: hem["r451_ref"] = ratio("hsa-miR-451a", ref)
d23 = de.loc["hsa-miR-23a-3p"] if "hsa-miR-23a-3p" in de.index else None
if d23 is not None: say(f"miR-23a-3p in female-only DE: log2FC={d23.log2FoldChange:.3f} padj={d23.padj:.3g}")
for c in [c for c in hem.columns if c != "group"]:
    a, b = hem.loc[hem.group=="BC", c], hem.loc[hem.group=="HC", c]
    say(f"{c}: BC median {a.median():.3f} vs HC {b.median():.3f}  MWU p={stats.mannwhitneyu(a,b).pvalue:.3g}")
pubh = pub15.reindex([s.replace("BC0","BC").replace("BC00","BC") for s in hem.index])
hem.to_csv(R("hemolysis.tsv"), sep="\t")

# ---------- 6. panel scores & AUC (DeLong) ----------
def delong(y, s):
    y = np.asarray(y).astype(bool); s = np.asarray(s, float)
    pos, neg = s[y], s[~y]; m, n = len(pos), len(neg)
    V10 = np.array([((p > neg).sum() + 0.5 * (p == neg).sum()) / n for p in pos])
    V01 = np.array([((pos > q).sum() + 0.5 * (pos == q).sum()) / m for q in neg])
    auc = V10.mean(); se = np.sqrt(V10.var(ddof=1) / m + V01.var(ddof=1) / n)
    return auc, max(0, auc - 1.96 * se), min(1, auc + 1.96 * se)
z = lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, np.nan), axis=0)
def score(up, down=()):
    up = [m for m in up if m in z.index]; down = [m for m in down if m in z.index]
    s = z.loc[up].mean() if up else 0
    if down: s = s - z.loc[down].mean()
    return s, len(up), len(down)
y = (grp == "BC").values
out = []
for name, u, d_ in [("Up8 panel", UP, ()), ("Down11 panel", (), DOWN), ("LASSO5", LASSO, ()), ("High-abundance 3", HIGH3, ()),
                    ("Low-count LASSO (138-5p, 32-3p)", ["hsa-miR-138-5p", "hsa-miR-32-3p"], ())]:
    s, nu, nd = score(u, d_)
    auc, lo, hi = delong(y, s)
    out.append([name, nu, nd, auc, lo, hi])
au = pd.DataFrame(out, columns=["score", "n_up", "n_down", "AUC", "CI_low", "CI_high"])
au.to_csv(R("panel_auc_female_only.tsv"), sep="\t", index=False)
say("\n== Panel AUCs, female-only 152 vs 48 (note: 76 cases overlap original discovery) =="); say(au.round(3).to_string(index=False))
# depth-AUC: does library depth alone separate groups?
for c in ["umi_total", "reads"]:
    auc, lo, hi = delong(y, qc.reindex(samples)[c].values); say(f"AUC of library depth ({c}) alone: {auc:.3f} ({lo:.3f}-{hi:.3f})")
say("\nDONE", time.strftime("%H:%M:%S"))
