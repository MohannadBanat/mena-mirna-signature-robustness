#!/usr/bin/env python3
"""Amendment v5a, A1 item 4: Figure 3 redrawn with counts per million and DESeq2 median-of-ratios side by side."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from sklearn.metrics import roc_curve, roc_auc_score
from PIL import Image
for f in fm.findSystemFonts(["/usr/share/texmf/fonts/opentype/public/tex-gyre"]): fm.fontManager.addfont(f)
plt.rcParams.update({"font.family": "TeX Gyre Pagella", "font.size": 7.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False})
BC, HC, GR, BL, YL, PU = "#D55E00", "#009E73", "#999999", "#0072B2", "#E69F00", "#CC79A7"
W = "/home/claude/rev"; O = f"{W}/revision_2026-11-21"
M = pd.read_csv(f"{W}/main/results/counts_matrix_umi_mirbase22.tsv", sep="\t", index_col=0); S = sorted(M.columns); M = M[S].astype(float)
y = np.array([s.startswith("BC") for s in S])
sig = pd.read_csv(f"{W}/main/results/locked_signature_v1.csv").set_index("miRNA ID")
UP = sig.index[sig.direction == "up"].tolist(); DN = sig.index[sig.direction == "down"].tolist(); L5 = sig.index[sig.LASSO5].tolist()
def mor(M):
    pos = M[(M > 0).all(1)]; lg = np.log(pos); return M / np.exp(lg.sub(lg.mean(1), axis=0).median())
NORM = {"Counts per million": np.log2(M / M.sum() * 1e6 + 1), "DESeq2 median-of-ratios": np.log2(mor(M) + 1)}
def scores(lc):
    z = lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, np.nan), axis=0)
    return {"Up-regulated (8)": z.loc[[m for m in UP if m in z.index]].mean().values,
            "Down-regulated (9)": -z.loc[[m for m in DN if m in z.index]].mean().values,
            "LASSO5": z.loc[L5].mean().values}
SC = {k: scores(v) for k, v in NORM.items()}
fig = plt.figure(figsize=(16.6 / 2.54, 11.0 / 2.54)); gs = fig.add_gridspec(2, 2, width_ratios=[1.65, 1], wspace=0.35, hspace=0.55)
xi = np.arange(len(S))
for r, (nm, sc) in enumerate(SC.items()):
    a = fig.add_subplot(gs[r, 0]); las = sc["LASSO5"]
    a.axvspan(142.5, 151.5, color=YL, alpha=0.22, lw=0)
    a.scatter(xi[y], las[y], s=5, color=BC, lw=0, label="Cancer"); a.scatter(xi[~y], las[~y], s=5, color=HC, lw=0, marker="s", label="Control")
    a.axhline(0, color=GR, lw=0.5, ls=":")
    a.set_xticks([0, 49, 99, 151, 199]); a.set_xticklabels(["BC001", "BC050", "BC100", "BC152", "HC048"], fontsize=6.5)
    a.set_ylabel("LASSO5 score (mean z)"); a.set_title(nm, fontsize=7.5, loc="left")
    if r == 0: a.legend(frameon=False, fontsize=6.5, loc="lower left", ncol=2)
    if r == 1: a.set_xlabel("Library, in deposited sample order")
    a.text(-0.1, 1.08, "AB"[r], transform=a.transAxes, fontweight="bold", fontsize=9)
b = fig.add_subplot(gs[:, 1])
for lab, col in [("Up-regulated (8)", PU), ("LASSO5", BC), ("Down-regulated (9)", BL)]:
    for nm, ls in [("Counts per million", "-"), ("DESeq2 median-of-ratios", "--")]:
        s = SC[nm][lab]; fp, tp, _ = roc_curve(y, s)
        b.plot(fp, tp, color=col, ls=ls, lw=1.0, label=f"{lab}, {'CPM' if ls == '-' else 'DESeq2'}: {roc_auc_score(y, s):.2f}")
b.plot([0, 1], [0, 1], color=GR, lw=0.5, ls=":"); b.set_aspect("equal"); b.set_xlim(0, 1); b.set_ylim(0, 1.01)
b.set_xlabel("1 − specificity"); b.set_ylabel("Sensitivity"); b.legend(frameon=False, fontsize=6.0, loc="lower right", handlelength=1.6)
b.text(-0.3, 1.02, "C", transform=b.transAxes, fontweight="bold", fontsize=9)
fig.subplots_adjust(left=0.07, right=0.975, bottom=0.1, top=0.93)
p = f"{O}/Figure3_revised.png"; fig.savefig(p, dpi=600); im = Image.open(p).convert("RGB"); im.save(p, dpi=(600, 600)); print(p, im.size)
for nm, sc in SC.items(): print(nm, {k: round(roc_auc_score(y, v), 3) for k, v in sc.items()})
