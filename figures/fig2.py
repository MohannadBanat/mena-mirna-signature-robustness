import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
for f in fm.findSystemFonts(["/usr/share/texmf/fonts/opentype/public/tex-gyre"]): fm.fontManager.addfont(f)
plt.rcParams.update({"font.family": "TeX Gyre Pagella", "font.size": 7.5, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False})
# Okabe-Ito colour-blind-safe palette
OI = {"bc": "#D55E00", "ben": "#0072B2", "ctl": "#009E73", "grey": "#999999"}
# LASSO5 direction-only score quantiles (5, 25, 50, 75, 95%) from GSE73002 (z-scored to non-cancer controls)
Q = {"Breast cancer\n(n = 1,280)": ([0.310, 1.838, 2.426, 2.992, 3.871], OI["bc"]),
     "Benign breast\ndisease (n = 54)": ([1.233, 1.609, 2.186, 2.798, 3.793], OI["ben"]),
     "Non-cancer\n(n = 2,686)": ([-1.275, -0.536, 0.000, 0.518, 1.260], OI["ctl"])}
# ROC points from the full curves (gse_roc_points.tsv)
import csv
pts={}
for r in csv.reader(open("/home/claude/ms/gse_roc_points.tsv"), delimiter="\t"):
    if len(r)==3 and not r[0].startswith("SCORE"): pts.setdefault(r[0],([],[])); pts[r[0]][0].append(float(r[1])); pts[r[0]][1].append(float(r[2]))
roc = {"Cancer vs non-cancer, AUC 0.951": (*pts["BC vs control"], OI["bc"], "-"),
       "Benign vs non-cancer, AUC 0.984": (*pts["Benign vs control"], OI["ben"], "-"),
       "Cancer vs benign, AUC 0.558": (*pts["BC vs benign"], "#000000", "--")}
fig, ax = plt.subplots(1, 2, figsize=(16.6 / 2.54, 6.4 / 2.54), gridspec_kw={"width_ratios": [1, 1.05]})
a = ax[0]
for i, (lab, (q, col)) in enumerate(Q.items()):
    a.plot([i, i], [q[0], q[4]], color=col, lw=0.9)
    a.add_patch(plt.Rectangle((i - 0.22, q[1]), 0.44, q[3] - q[1], facecolor=col, alpha=0.25, edgecolor=col, lw=0.8))
    a.plot([i - 0.22, i + 0.22], [q[2], q[2]], color=col, lw=1.4)
a.set_xticks(range(3)); a.set_xticklabels(list(Q.keys())); a.set_xlim(-0.6, 2.6)
a.set_ylabel("LASSO5 score (z units vs non-cancer)"); a.axhline(0, color=OI["grey"], lw=0.5, ls=":")
a.text(-0.17, 1.03, "A", transform=a.transAxes, fontweight="bold", fontsize=9)
b = ax[1]
for lab, (x, y, col, ls) in roc.items(): b.plot(x, y, color=col, ls=ls, lw=1.1, label=lab)
b.plot([0, 1], [0, 1], color=OI["grey"], lw=0.5, ls=":")
b.set_xlabel("1 − specificity"); b.set_ylabel("Sensitivity"); b.set_xlim(0, 1); b.set_ylim(0, 1.01); b.set_aspect("equal")
b.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.04, 0.0), fontsize=6.3, handlelength=1.4, borderaxespad=0.2)
b.text(-0.2, 1.03, "B", transform=b.transAxes, fontweight="bold", fontsize=9)
fig.tight_layout(w_pad=2.2)
fig.savefig("/home/claude/ms/fig5.png", dpi=600)
print("ok")
