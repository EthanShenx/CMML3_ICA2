"""
render_panels_py.py
Render the Fig S2 isoform sub-panels (a-g) and Fig S3 panel (d) radar as
individual PNG files in the matplotlib style of 07_make_figures.py.

Inputs are read from data2plot/ CSVs (no master_results.csv required), so
this script reflects whatever threshold the CSVs were exported at.

Run from the project root or scripts/:
    python scripts/render_panels_py.py
Outputs land in figures/panels_py/.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT       = SCRIPT_DIR.parent
DATA2PLOT  = ROOT / "data2plot"
OUT        = ROOT / "figures" / "panels_py"
OUT.mkdir(parents=True, exist_ok=True)

# Constants mirroring 07_make_figures.py
PROTEINS = ["A8", "A9", "A10", "A11", "A12", "A13", "A14"]
CAT_ORDER = [
    "Short-chain", "Medium-chain", "α,β-Unsaturated",
    "Aromatic", "Retinoid", "Semialdehyde", "Other",
]
CAT_COLORS = {
    "Short-chain":    "#4CAF50",
    "Medium-chain":   "#2196F3",
    "α,β-Unsaturated": "#FF9800",
    "Aromatic":       "#9C27B0",
    "Retinoid":       "#F44336",
    "Semialdehyde":   "#00BCD4",
    "Other":          "#9E9E9E",
}
PROT_COLORS = {
    "A8":  "#a6cee3", "A9":  "#1f78b4", "A10": "#b2df8a",
    "A11": "#33a02c", "A12": "#fb9a99", "A13": "#e31a1c",
    "A14": "#fdbf6f",
}

DG_THRESHOLD = -7.0

matplotlib.rcParams["font.family"]  = "Arial"
matplotlib.rcParams["pdf.fonttype"] = 42


def short_name(ligand, maxlen=22):
    name = ligand.split("_", 1)[-1] if "_" in ligand else ligand
    return name if len(name) <= maxlen else name[: maxlen - 1] + "…"


def lbl_si(ax, text, fontsize=18):
    ax.text(-0.08, 1.04, f"({text})", transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="bottom", ha="left",
            family="Arial")


def render_figS2_panel(prot, panel_letter, df):
    sub = df[df["protein"] == prot].sort_values("DG").reset_index(drop=True)
    colors_s = [CAT_COLORS.get(c, "#9E9E9E") for c in sub["category"]]

    fig, ax = plt.subplots(figsize=(5.5, 6))
    ax.barh(range(len(sub)), sub["DG"], color=colors_s, height=0.75, linewidth=0)

    for j, row in sub.iterrows():
        if pd.notna(row.get("distance")) and row["distance"] <= 4.0:
            ax.text(row["DG"] - 0.04, j, "★",
                    ha="right", va="center", fontsize=6, color="#FFD700")

    if len(sub):
        top = sub.iloc[0]
        ax.text(0.03, 0.97,
                f"top: {short_name(top['ligand'], 24)}",
                transform=ax.transAxes,
                ha="left", va="top", fontsize=10, color="#222")

    ax.axvline(DG_THRESHOLD, ls="--", color="grey", lw=0.8, alpha=0.7)
    ax.set_yticks([])
    ax.set_xlabel("ΔG (kcal mol$^{-1}$)", fontsize=13)
    ax.tick_params(axis="x", labelsize=11)
    ax.set_title(f"{prot}  (median ΔG = {sub['DG'].median():.2f})",
                 fontsize=13, pad=6, color=PROT_COLORS[prot])
    lbl_si(ax, panel_letter)

    name = f"figS2{panel_letter}_{prot}.png"
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ {name}")


def render_figS2():
    csv = DATA2PLOT / "figS2_data.csv"
    df = pd.read_csv(csv)
    for prot, letter in zip(PROTEINS, "abcdefg"):
        render_figS2_panel(prot, letter, df)


def render_figS3d():
    csv = DATA2PLOT / "figS3d_radar.csv"
    df = pd.read_csv(csv)
    cats_r = [c for c in CAT_ORDER if c != "Other"]
    N = len(cats_r)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, polar=True)

    for prot in PROTEINS:
        vals = []
        for cat in cats_r:
            row = df[(df["protein"] == prot) & (df["category"] == cat)]
            vals.append(float(row["frac"].iloc[0]) if len(row) else 0.0)
        vals += vals[:1]
        ax.plot(angles, vals, color=PROT_COLORS[prot],
                lw=1.5, label=prot, alpha=0.85)
        ax.fill(angles, vals, color=PROT_COLORS[prot], alpha=0.08)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(cats_r, size=11)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["25%", "50%", "75%", "100%"], size=10)
    ax.set_ylim(0, 1)
    thr_txt = f"$\\leq {DG_THRESHOLD:g}$"
    ax.set_title(f"Selectivity profile\n(fraction with $\\Delta G$ {thr_txt} kcal mol$^{{-1}}$)",
                 fontsize=13, pad=34)
    ax.legend(loc="center left", bbox_to_anchor=(1.18, 0.5),
              fontsize=10, framealpha=0.85)
    lbl_si(ax, "d")

    name = "figS3d_radar.pdf"
    fig.savefig(OUT / name, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ {name}")


if __name__ == "__main__":
    render_figS2()
    render_figS3d()
