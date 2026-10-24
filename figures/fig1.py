import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
for f in fm.findSystemFonts(["/usr/share/texmf/fonts/opentype/public/tex-gyre"]): fm.fontManager.addfont(f)
plt.rcParams.update({"font.family": "TeX Gyre Pagella", "font.size": 7.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False})
BC, HC, GR = "#D55E00", "#009E73", "#999999"
rows = [r.split(",") for r in open("/home/claude/ms/qc_compact.txt").read().strip().split(";")]
g = np.array([r[0][:2] for r in rows]); reads = np.array([float(r[1]) for r in rows]) / 1e3; umi = np.array([float(r[2]) for r in rows]) / 1e3
rng = np.random.default_rng(1)
fig = plt.figure(figsize=(16.6 / 2.54, 6.2 / 2.54))
gs = fig.add_gridspec(1, 3, width_ratios=[0.75, 0.75, 1.35], wspace=0.95)
def strip(ax, y, ylab, letter, ptxt):
    for i, (k, col, lab) in enumerate([("BC", BC, "Cancer\n(n = 152)"), ("HC", HC, "Control\n(n = 48)")]):
        v = y[g == k]
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=4, color=col, alpha=0.55, lw=0)
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.plot([i - 0.28, i + 0.28], [med, med], color="black", lw=1.2)
        ax.plot([i, i], [q1, q3], color="black", lw=0.8)
    ax.set_yscale("log"); ax.set_xticks([0, 1]); ax.set_xticklabels(["Cancer\n(n = 152)", "Control\n(n = 48)"])
    ax.set_xlim(-0.6, 1.6); ax.set_ylabel(ylab)
    ax.text(0.5, 1.02, ptxt, transform=ax.transAxes, ha="center", va="bottom", fontsize=7)
    ax.text(-0.42, 1.03, letter, transform=ax.transAxes, fontweight="bold", fontsize=9)
strip(fig.add_subplot(gs[0]), reads, "Raw reads (millions)", "A", "$p$ < 0.001")
strip(fig.add_subplot(gs[1]), umi, "UMI-collapsed molecules (millions)", "B", "$p$ = 0.048")
ax = fig.add_subplot(gs[2])
# detection after depth matching (mean of 100 thinned datasets); published direction in brackets
D = [("miR-4455 (up)", .630, .000), ("miR-138-5p (up)*", .585, .122), ("miR-32-3p (up)*", .452, .168), ("miR-766-5p (up)", .452, .179),
     ("miR-382-3p (down)", .473, .640), ("miR-1976 (down)", .448, .356), ("miR-6741-3p (down)", .287, .361), ("miR-500b-5p (down)", .241, .189)]
for i, (lab, b, h) in enumerate(D[::-1]):
    ax.plot([h * 100, b * 100], [i, i], color=GR, lw=0.8, zorder=1)
    ax.scatter(b * 100, i, color=BC, s=14, zorder=2, label="Cancer" if i == 0 else None)
    ax.scatter(h * 100, i, color=HC, s=14, zorder=2, marker="s", label="Control" if i == 0 else None)
ax.set_yticks(range(len(D))); ax.set_yticklabels([d[0] for d in D[::-1]]); ax.set_xlim(-3, 75)
ax.set_xlabel("Samples with detection (%)")
ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=6.8, handletextpad=0.3, borderaxespad=0.1)
ax.text(-0.62, 1.03, "C", transform=ax.transAxes, fontweight="bold", fontsize=9)
fig.subplots_adjust(left=0.08, right=0.99, bottom=0.2, top=0.88)
fig.savefig("/home/claude/ms/fig1.png", dpi=600)
from PIL import Image
im = Image.open("/home/claude/ms/fig1.png").convert("RGB"); im.save("/home/claude/ms/fig1.png", dpi=(600, 600))
im.resize((im.size[0] // 5, im.size[1] // 5)).save("/tmp/claude-0/fig1_preview.png"); print(im.size)
