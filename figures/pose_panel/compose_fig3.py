import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patheffects as pe
from pathlib import Path

HERE = Path(__file__).parent
OUT  = HERE.parent / "fig3_pose_comparison"

a8  = mpimg.imread(HERE / "panel_A8_productive.png")
a12 = mpimg.imread(HERE / "panel_A12_nonproductive.png")

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), dpi=300)
for ax, img, title, dist, prod_label, prod_color in [
    (axes[0], a8,  r"A8 $\times$ $\beta$-carotene aldehyde",  2.27, "productive",     "#2E7D32"),
    (axes[1], a12, r"A12 $\times$ $\beta$-carotene aldehyde", 5.49, "non-productive", "#C62828"),
]:
    ax.imshow(img)
    ax.axis("off")
    ax.set_title(title, fontsize=11, pad=6, fontweight="bold")
    # Distance / verdict box
    txt = "Cys-SG ↔ aldehyde-C: {:.2f} Å".format(dist)
    ax.text(0.03, 0.97, txt, transform=ax.transAxes,
            fontsize=10, va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="grey", lw=0.6, alpha=0.92))
    ax.text(0.97, 0.97, prod_label, transform=ax.transAxes,
            fontsize=10, fontweight="bold", color="white", va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.30", fc=prod_color, ec="none", alpha=0.95))

axes[0].text(-0.02, 1.02, "(a)", transform=axes[0].transAxes,
             fontsize=13, fontweight="bold", va="top", ha="left")
axes[1].text(-0.02, 1.02, "(b)", transform=axes[1].transAxes,
             fontsize=13, fontweight="bold", va="top", ha="left")

# Legend strip at the bottom
legend_y = -0.02
legend_kw = dict(transform=fig.transFigure, fontsize=8, ha="center", va="center")
fig.text(0.5, 0.025,
         "Cartoon: receptor (grey, transparent).  Sticks: catalytic Cys side chain (yellow C, orange S) and "
         "β-carotene-derived polyene-aldehyde ligand.  Spheres: Cys-SG (orange) and aldehyde carbon (firebrick); "
         "dashed line: Cys-SG ↔ aldehyde-C distance.",
         **legend_kw)

plt.tight_layout(rect=[0, 0.05, 1, 1])
fig.savefig(str(OUT) + ".pdf", dpi=300, bbox_inches="tight")
fig.savefig(str(OUT) + ".png", dpi=300, bbox_inches="tight")
print("saved", OUT)
