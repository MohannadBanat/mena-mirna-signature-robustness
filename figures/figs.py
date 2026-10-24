import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_curve, roc_auc_score
from PIL import Image
for f in fm.findSystemFonts(["/usr/share/texmf/fonts/opentype/public/tex-gyre"]): fm.fontManager.addfont(f)
plt.rcParams.update({"font.family": "TeX Gyre Pagella", "font.size": 7.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False, "mathtext.fontset": "custom",
                     "mathtext.it": "TeX Gyre Pagella:italic", "mathtext.rm": "TeX Gyre Pagella"})
BC, HC, GR, BL, YL, PU, SK = "#D55E00", "#009E73", "#999999", "#0072B2", "#E69F00", "#CC79A7", "#56B4E9"
D = "/home/claude/ms/data/"
def rd(f): return [r.split(",") for r in open(D + f).read().strip().split(";")]
qc = {r[0]: [float(x) for x in r[1:]] for r in rd("qc_full.txt")}
sc = {r[0]: [float(x) for x in r[1:]] for r in rd("scores.txt")}
hem = {r[0]: [float(r[2]), float(r[3])] for r in rd("hem.txt")}
S = sorted(qc); y = np.array([s.startswith("BC") for s in S]); rng = np.random.default_rng(1)
def L(ax, t, x=-0.2): ax.text(x, 1.04, t, transform=ax.transAxes, fontweight="bold", fontsize=9, va="bottom")
def save(fig, n):
    p = f"/home/claude/ms/fig{n}.png"; fig.savefig(p, dpi=600)
    im = Image.open(p).convert("RGB"); im.save(p, dpi=(600, 600)); print(n, im.size)
def strip(ax, v, ylab, log=True):
    for i, (k, col) in enumerate([(True, BC), (False, HC)]):
        w = v[y == k]
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(w)), w, s=3.5, color=col, alpha=0.55, lw=0)
        q1, m, q3 = np.percentile(w, [25, 50, 75])
        ax.plot([i - 0.28, i + 0.28], [m, m], color="k", lw=1.2); ax.plot([i, i], [q1, q3], color="k", lw=0.8)
    if log: ax.set_yscale("log")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Cancer\n(n = 152)", "Control\n(n = 48)"]); ax.set_xlim(-0.6, 1.6); ax.set_ylabel(ylab)
    p = mannwhitneyu(v[y], v[~y]).pvalue
    ax.text(0.5, 1.0, "$p$ < 0.001" if p < 0.001 else f"$p$ = {p:.3f}".replace("0.04", "0.04"), transform=ax.transAxes, ha="center", va="bottom", fontsize=7)

# ---------- Figure 1: the count matrix ----------
reads = np.array([qc[s][0] for s in S]); kept = np.array([qc[s][1] for s in S]); noad = np.array([qc[s][2] for s in S])
badl = np.array([qc[s][3] for s in S]); ali = np.array([qc[s][4] for s in S]); umi = np.array([qc[s][5] for s in S]); det = np.array([qc[s][6] for s in S])
fig = plt.figure(figsize=(16.6 / 2.54, 12.4 / 2.54))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.75, wspace=0.75)
a = fig.add_subplot(gs[0, 0:2])
cats = [("No adapter", noad, GR), ("Insert outside 16–30 nt", badl, "#CCCCCC"), ("Unaligned", kept - ali, YL),
        ("PCR duplicate", ali - umi, SK), ("Unique molecule (UMI)", umi, BL)]
for j, k in enumerate([True, False]):
    left = 0; tot = reads[y == k].sum()
    for lab, v, col in cats:
        w = v[y == k].sum() / tot * 100
        a.barh(j, w, left=left, color=col, height=0.55, edgecolor="white", lw=0.4, label=lab if j == 0 else None)
        if w > 4: a.text(left + w / 2, j, f"{w:.0f}%", ha="center", va="center", fontsize=6.5, color="k")
        left += w
a.set_yticks([0, 1]); a.set_yticklabels(["Cancer", "Control"]); a.invert_yaxis(); a.set_xlim(0, 100)
a.set_xlabel("Share of all raw reads in the group (%)")
a.legend(frameon=False, ncol=3, fontsize=6.3, loc="lower left", bbox_to_anchor=(0, 1.02), handlelength=1, columnspacing=1)
L(a, "A", -0.1)
strip(fig.add_subplot(gs[0, 2]), umi / 1e6, "UMI molecules (millions)"); L(fig.axes[-1], "B", -0.45)
top = rd("top15.txt")
a = fig.add_subplot(gs[1, 0:2])
names = [t[0].replace("hsa-", "") for t in top][:10]
x = np.arange(len(names)); bcv = [float(t[2]) for t in top][:10]; hcv = [float(t[3]) for t in top][:10]
a.bar(x - 0.2, bcv, 0.4, color=BC, label="Cancer"); a.bar(x + 0.2, hcv, 0.4, color=HC, label="Control")
a.set_yscale("log"); a.set_xticks(x); a.set_xticklabels(names, rotation=40, ha="right", fontsize=6.6)
a.set_ylabel("Share of molecules (%)"); a.legend(frameon=False, fontsize=6.5, loc="upper right")
L(a, "C", -0.1)
strip(fig.add_subplot(gs[1, 2]), det, "miRNAs detected per library", log=False); L(fig.axes[-1], "D", -0.45)
fig.subplots_adjust(left=0.08, right=0.99, bottom=0.14, top=0.9)
save(fig, 1)

# ---------- Figure 2: female-only differential expression ----------
de = [r for r in rd("de.txt")]
lfc = np.array([float(r[2]) for r in de]); padj = np.array([float(r[4]) for r in de]); pub = np.array([float(r[5]) if r[5] not in ("", "nan") else np.nan for r in de])
nm = [r[0].replace("hsa-", "") for r in de]
panel = {"let-7g-5p", "miR-24-3p", "miR-142-3p", "miR-223-3p", "miR-4455", "miR-106b-3p", "miR-30a-5p"}
fig, ax = plt.subplots(1, 2, figsize=(16.6 / 2.54, 7.0 / 2.54))
a = ax[0]; lp = -np.log10(np.clip(padj, 1e-40, 1)); sig = (padj < 0.05) & (abs(lfc) >= 1)
a.scatter(lfc[~sig], lp[~sig], s=4, color=GR, alpha=0.6, lw=0, label="Not DE")
a.scatter(lfc[sig & (lfc > 0)], lp[sig & (lfc > 0)], s=6, color=BC, lw=0, label="Up in cancer (25)")
a.scatter(lfc[sig & (lfc < 0)], lp[sig & (lfc < 0)], s=6, color=BL, lw=0, label="Down in cancer (12)")
for i, n in enumerate(nm):
    if n in panel: a.scatter(lfc[i], lp[i], s=16, facecolor="none", edgecolor="k", lw=0.6); a.annotate(n, (lfc[i], lp[i]), xytext=(3, 2), textcoords="offset points", fontsize=5.8)
a.axhline(-np.log10(0.05), color=GR, lw=0.5, ls=":"); [a.axvline(v, color=GR, lw=0.5, ls=":") for v in (-1, 1)]
a.set_xlabel("log$_2$ fold change, cancer vs control"); a.set_ylabel("−log$_{10}$ adjusted $p$")
a.legend(frameon=False, fontsize=6.2, loc="upper left", handletextpad=0.2); L(a, "A", -0.17)
b = ax[1]; ok = ~np.isnan(pub)
b.scatter(pub[ok & ~sig], lfc[ok & ~sig], s=4, color=GR, alpha=0.6, lw=0)
b.scatter(pub[ok & sig], lfc[ok & sig], s=6, color=PU, lw=0, label="DE in female-only analysis")
for i, n in enumerate(nm):
    if n in panel and ok[i]:
        b.scatter(pub[i], lfc[i], s=16, facecolor="none", edgecolor="k", lw=0.6)
        if n in ("miR-4455", "miR-30a-5p", "miR-106b-3p"): b.annotate(n, (pub[i], lfc[i]), xytext=(4, -3 if n != "miR-106b-3p" else -8), textcoords="offset points", fontsize=5.8)
b.annotate("let-7g-5p, miR-24-3p,\nmiR-142-3p, miR-223-3p", (1.05, 1.1), xytext=(1.6, 2.6), fontsize=5.8, arrowprops=dict(arrowstyle="-", lw=0.5, color="k"))
lim = [min(np.nanmin(pub), lfc.min()) - 0.3, max(np.nanmax(pub), lfc.max()) + 0.3]
b.plot(lim, lim, color=GR, lw=0.5, ls=":"); b.axhline(0, color=GR, lw=0.4); b.axvline(0, color=GR, lw=0.4)
b.set_xlabel("Published log$_2$ fold change (mixed-sex discovery)"); b.set_ylabel("log$_2$ fold change (female-only)")
b.text(0.03, 0.97, "Spearman ρ = 0.885, $n$ = 295", transform=b.transAxes, va="top", fontsize=6.8)
b.legend(frameon=False, fontsize=6.2, loc="lower right"); L(b, "B", -0.17)
fig.tight_layout(w_pad=2.0); save(fig, 2)

# ---------- Figure 3: female-only scores ----------
fig = plt.figure(figsize=(16.6 / 2.54, 6.6 / 2.54)); gs = fig.add_gridspec(1, 2, width_ratios=[1.65, 1], wspace=0.35)
a = fig.add_subplot(gs[0]); las = np.array([sc[s][0] for s in S])
xi = np.arange(len(S)); cl = np.array([s in {f"BC{i}" for i in range(144, 153)} for s in S])
a.axvspan(142.5, 151.5, color=YL, alpha=0.22, lw=0)
a.scatter(xi[y], las[y], s=5, color=BC, lw=0, label="Cancer"); a.scatter(xi[~y], las[~y], s=5, color=HC, lw=0, marker="s", label="Control")
a.axhline(0, color=GR, lw=0.5, ls=":")
a.set_xticks([0, 49, 99, 151, 152, 199]); a.set_xticklabels(["BC001", "BC050", "BC100", "BC152", "", "HC048"], fontsize=6.5)
a.text(147, las.max() + 0.05, "BC144–\nBC152", ha="center", va="bottom", fontsize=6.2)
a.set_ylim(las.min() - 0.2, las.max() + 0.55)
a.set_xlabel("Library, in deposited sample order"); a.set_ylabel("LASSO5 score (mean z)")
a.legend(frameon=False, fontsize=6.5, loc="lower left", ncol=2); L(a, "A", -0.1)
b = fig.add_subplot(gs[1])
for j, lab, col, ls in [(1, "Up-regulated 8", PU, "-"), (0, "LASSO5", BC, "-"), (2, "High-abundance 3", BL, "--")]:
    s = np.array([sc[k][j] for k in S]); fp, tp, _ = roc_curve(y, s)
    b.plot(fp, tp, color=col, ls=ls, lw=1.0, label=f"{lab}, AUC {roc_auc_score(y, s):.3f}")
b.plot([0, 1], [0, 1], color=GR, lw=0.5, ls=":"); b.set_aspect("equal"); b.set_xlim(0, 1); b.set_ylim(0, 1.01)
b.set_xlabel("1 − specificity"); b.set_ylabel("Sensitivity"); b.legend(frameon=False, fontsize=6.1, loc="lower right", handlelength=1.4)
L(b, "B", -0.3)
fig.subplots_adjust(left=0.07, right=0.99, bottom=0.17, top=0.9); save(fig, 3)

# ---------- Figure 4: depth matching and haemolysis ----------
fig = plt.figure(figsize=(16.6 / 2.54, 6.4 / 2.54)); gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 0.75, 0.75], wspace=0.8)
a = fig.add_subplot(gs[0])
Dm = [("miR-4455 (up)", .630, .000), ("miR-138-5p (up)*", .585, .122), ("miR-32-3p (up)*", .452, .168), ("miR-766-5p (up)", .452, .179),
      ("miR-382-3p (down)", .473, .640), ("miR-1976 (down)", .448, .356), ("miR-6741-3p (down)", .287, .361), ("miR-500b-5p (down)", .241, .189)]
for i, (lab, bb, h) in enumerate(Dm[::-1]):
    a.plot([h * 100, bb * 100], [i, i], color=GR, lw=0.8, zorder=1)
    a.scatter(bb * 100, i, color=BC, s=14, zorder=2, label="Cancer" if i == 0 else None)
    a.scatter(h * 100, i, color=HC, s=14, zorder=2, marker="s", label="Control" if i == 0 else None)
a.set_yticks(range(len(Dm))); a.set_yticklabels([d[0] for d in Dm[::-1]]); a.set_xlim(-3, 75)
a.set_xlabel("Libraries with detection (%)")
a.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=6.5, handletextpad=0.3, borderaxespad=0.1)
L(a, "A", -0.62)
strip(fig.add_subplot(gs[1]), np.array([hem[s][0] for s in S]), "log$_2$ miR-451a / miR-23a-3p", log=False); L(fig.axes[-1], "B", -0.55)
strip(fig.add_subplot(gs[2]), np.array([hem[s][1] for s in S]), "log$_2$ miR-451a / let-7a-5p", log=False); L(fig.axes[-1], "C", -0.55)
fig.subplots_adjust(left=0.165, right=0.99, bottom=0.2, top=0.88); save(fig, 4)
