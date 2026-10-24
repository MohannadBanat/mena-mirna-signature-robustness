import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from scipy.stats import spearmanr, mannwhitneyu
from PIL import Image
for f in fm.findSystemFonts(["/usr/share/texmf/fonts/opentype/public/tex-gyre"]): fm.fontManager.addfont(f)
plt.rcParams.update({"font.family": "TeX Gyre Pagella", "font.size": 7.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False, "mathtext.fontset": "custom",
                     "mathtext.it": "TeX Gyre Pagella:italic", "mathtext.rm": "TeX Gyre Pagella"})
BC, HC, GR, BL, YL, PU = "#D55E00", "#009E73", "#999999", "#0072B2", "#E69F00", "#CC79A7"
T = {}
for l in open('data/tissue_summary.tsv'):
    p = l.rstrip('\n').split('\t')
    if p[0] in ('UP', 'DN', 'PANEL') and len(p) > 3: T[(p[0], p[1])] = [float(x) for x in p[2:]]
plasma = {"let-7g-5p": 100, "miR-24-3p": 100, "miR-142-3p": 100, "miR-223-3p": 100, "miR-138-5p": 71, "miR-32-3p": 68, "miR-4455": 63, "miR-766-5p": 67,
          "miR-106b-3p": 100, "miR-30a-5p": 100, "miR-382-3p": 60, "miR-1976": 71, "miR-500b-5p": 47, "miR-6741-3p": 46, "miR-660-3p": 16, "miR-3195": 6, "miR-4497": 0, "miR-4639-5p": 0, "miR-106a-3p": 1}
order = ["let-7g-5p", "miR-24-3p", "miR-142-3p", "miR-223-3p", "miR-138-5p", "miR-32-3p", "miR-4455", "miR-766-5p",
         "miR-106b-3p", "miR-30a-5p", "miR-382-3p", "miR-1976", "miR-500b-5p", "miR-6741-3p", "miR-660-3p", "miR-3195", "miR-4497", "miR-4639-5p", "miR-106a-3p"]
lasso = {"let-7g-5p", "miR-24-3p", "miR-142-3p", "miR-138-5p", "miR-32-3p"}
fig = plt.figure(figsize=(16.6 / 2.54, 12.6 / 2.54))
gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], height_ratios=[1.25, 1], hspace=0.5, wspace=0.45)
a = fig.add_subplot(gs[:, 0])
for i, m in enumerate(order[::-1]):
    v = T.get(("PANEL", m), [0, 0, 0])
    tis = v[0] if len(v) > 1 else 0.0
    a.plot([plasma[m], tis], [i + 0.17, i - 0.17], color=GR, lw=0.7, zorder=1)
    a.scatter(plasma[m], i + 0.17, s=15, color=BC, zorder=2, label="Plasma, cancer cases" if i == 0 else None)
    a.scatter(tis, i - 0.17, s=15, color=BL, marker="D", zorder=3, label="Tumour tissue" if i == 0 else None)
a.set_yticks(range(len(order))); a.set_yticklabels([("*" if m in lasso else "") + m + (" (down)" if order.index(m) >= 8 else "") for m in order[::-1]], fontsize=6.6)
a.axhline(10.5, color=GR, lw=0.5, ls=":")
a.set_xlim(-4, 104); a.set_xlabel("Samples with ≥1 molecule (%)")
a.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, fontsize=6.5, handletextpad=0.3, columnspacing=0.8)
a.text(-0.42, 1.03, "A", transform=a.transAxes, fontweight="bold", fontsize=9)
S = [x.split(',') for x in open('data/tissue_samples.txt').read().strip().split(';')]
cols = ["let-7g-5p", "miR-24-3p", "miR-142-3p", "miR-223-3p", "miR-138-5p", "miR-32-3p", "miR-4455", "miR-766-5p", "miR-150-5p", "miR-21-5p", "miR-155-5p"]
X = {c: np.array([float(r[i + 2]) for r in S]) for i, c in enumerate(cols)}
b = fig.add_subplot(gs[0, 1])
b.scatter(X["miR-150-5p"], X["miR-142-3p"], s=6, color=PU, lw=0, label="miR-142-3p")
b.scatter(X["miR-150-5p"], X["miR-223-3p"], s=6, color=YL, lw=0, marker="s", label="miR-223-3p")
r1 = spearmanr(X["miR-150-5p"], X["miR-142-3p"])[0]; r2 = spearmanr(X["miR-150-5p"], X["miR-223-3p"])[0]
b.text(0.03, 0.97, f"ρ = {r1:.2f} (miR-142-3p)\nρ = {r2:.2f} (miR-223-3p)", transform=b.transAxes, va="top", fontsize=6.5)
b.set_xlabel("miR-150-5p, log$_2$(MPM + 1)"); b.set_ylabel("log$_2$(MPM + 1)")
b.legend(frameon=False, fontsize=6.3, loc="lower right", handletextpad=0.2)
b.text(-0.3, 1.03, "B", transform=b.transAxes, fontweight="bold", fontsize=9)
c = fig.add_subplot(gs[1, 1])
up = [v[2] for (g, m), v in T.items() if g == 'UP']; dn = [v[2] for (g, m), v in T.items() if g == 'DN']
rng = np.random.default_rng(3)
for i, (vals, col, lab) in enumerate([(up, BC, f"Up in plasma\n(n = {len(up)})"), (dn, BL, f"Down in plasma\n(n = {len(dn)})")]):
    y = np.log10(np.array(vals) + 0.1)
    c.scatter(i + rng.uniform(-0.15, 0.15, len(y)), y, s=7, color=col, lw=0, alpha=0.8)
    c.plot([i - 0.25, i + 0.25], [np.median(y)] * 2, color="k", lw=1.1)
c.set_xticks([0, 1]); c.set_xticklabels([f"Up in plasma\n(n = {len(up)})", f"Down in plasma\n(n = {len(dn)})"]); c.set_xlim(-0.6, 1.6)
c.set_yticks([-1, 0, 1, 2, 3, 4]); c.set_yticklabels(["0", "1", "10", "100", "10$^3$", "10$^4$"])
c.set_ylabel("Median MPM in tumours")
c.text(0.5, 1.0, f"$p$ = {mannwhitneyu(up, dn).pvalue:.3f}", transform=c.transAxes, ha="center", va="bottom", fontsize=7)
c.text(-0.3, 1.03, "C", transform=c.transAxes, fontweight="bold", fontsize=9)
fig.subplots_adjust(left=0.17, right=0.99, bottom=0.1, top=0.92)
fig.savefig("fig6.png", dpi=600)
im = Image.open("fig6.png").convert("RGB"); im.save("fig6.png", dpi=(600, 600)); print(im.size)
im.resize((im.size[0] * 2 // 7, im.size[1] * 2 // 7)).save("/tmp/claude-0/prev6.png")
