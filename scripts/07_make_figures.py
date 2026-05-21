"""
07_make_figures.py
Generate all manuscript figures from master_results.csv for the ALDH-aldehyde
docking study.  Produces fig1_landscape, fig2_geometry, figS1_qc, figS2_isoforms,
figS3_chemspace (PDF + PNG, 300 dpi) and table_S1_dual_winners.tex.

Run from the scripts/ directory or project root:
    python scripts/07_make_figures.py
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
import seaborn as sns
from scipy.stats import spearmanr, mannwhitneyu, kruskal, pearsonr
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT       = SCRIPT_DIR.parent
SDF_DIR    = ROOT.parent / "Data" / "substrate"          # ICA2/Data/substrate/
RESULTS    = ROOT / "results"
FIGS       = ROOT / "figures"
REPORT     = ROOT / "report"
FIGS.mkdir(exist_ok=True)
REPORT.mkdir(exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
df      = pd.read_csv(RESULTS / "master_results.csv")
df_ald  = df[df["is_aldehyde"] == True].copy()

PROTEINS = ["A8", "A9", "A10", "A11", "A12", "A13", "A14"]

# ── Substrate categories ───────────────────────────────────────────────────────
CATEGORY_MAP = {
    # Short-chain saturated (C1–C5)
    "1_formaldehyde":              "Short-chain",
    "2_acetaldehyde":              "Short-chain",
    "3_propylaldehyde":            "Short-chain",
    "4_malonaldehyde":             "Short-chain",
    "5_butyraldehyde":             "Short-chain",
    "6_amylaldehyde":              "Short-chain",
    "9_2-hydroxypropionaldehyde":  "Short-chain",
    # Medium-chain aliphatic (C6–C16)
    "7_glutaraldehyde":            "Medium-chain",
    "8_hexanal":                   "Medium-chain",
    "15_heptanal":                 "Medium-chain",
    "16_octanaldehyde":            "Medium-chain",
    "17_nonanal":                  "Medium-chain",
    "18_decanal":                  "Medium-chain",
    "19_dodecanal":                "Medium-chain",
    "24_Tetradecanal":             "Medium-chain",
    "25_hexadecanal":              "Medium-chain",
    # α,β-Unsaturated / lipid-peroxidation
    "11_prop-2-enal":              "α,β-Unsaturated",
    "12_prop-2-ynal":              "α,β-Unsaturated",
    "13_(2E)-but-2-enal":          "α,β-Unsaturated",
    "14_trans-2-Hexenal":          "α,β-Unsaturated",
    "20_Citral":                   "α,β-Unsaturated",
    "21_cis-2-Octenal":            "α,β-Unsaturated",
    "22_trans-2-Octenal":          "α,β-Unsaturated",
    "23_4-Hydroxynonenal":         "α,β-Unsaturated",
    "26_trans-2-hexadecenal":      "α,β-Unsaturated",
    "34_4-hydroxy-2-nonenal":      "α,β-Unsaturated",
    # Aromatic / heteroaromatic
    "27_Benzaldehyde":             "Aromatic",
    "28_Phenylacetaldehyde":       "Aromatic",
    "29_Salicylaldehyde":          "Aromatic",
    "30_(E)-3-(3,4-Dihydroxyphenyl)prop-2-enal": "Aromatic",
    "35_Cinnamaldehyde":           "Aromatic",
    "36_3-(4-Isopropylphenyl)-2-methylpropanal":  "Aromatic",
    "37_thiobenzaldehyde":         "Aromatic",
    "38_(2R)-2-[(2S,3R,4S)-3,4-Dihydroxy-5-oxo -tetrahydrofuran-2-yl]-2-hydroxy-acetaldehyde": "Aromatic",
    "39_4-nitrobenzaldehyde":      "Aromatic",
    "40_Furan-2-carbaldehyde":     "Aromatic",
    "41_Indole-3-acetaldehyde":    "Aromatic",
    # Retinoids / polyenes
    "31_(2E,4E,6E,8E)-3,7-dimethyl-9-(2,6,6-trimethylcyclohexen-1-yl)nona-2,4,6,8-tetraenal": "Retinoid",
    "32_9-cis-Retinal":            "Retinoid",
    "33_11-cis-Retinol":           "Retinoid",
    # Amino-acid semialdehydes
    "10_4-guanidinobutyraldehyde": "Semialdehyde",
    "42_Glutamate-5-semialdehyde": "Semialdehyde",
    "43_Glutamate-1-semialdehyde": "Semialdehyde",
    "44_glutarate semialdehyde":   "Semialdehyde",
    "45_adipate semialdehyde":     "Semialdehyde",
    "46_Succinic semialdehyde":    "Semialdehyde",
    "47_methylmalonate semialdehyde": "Semialdehyde",
    "48_alpha-aminoadipate semialdehyde": "Semialdehyde",
    "49_2-Aminomuconic semialdehyde": "Semialdehyde",
}
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
# RColorBrewer "Paired" — fixed palette used consistently across all figures
PROT_COLORS = {
    "A8":  "#a6cee3", "A9":  "#1f78b4", "A10": "#b2df8a",
    "A11": "#33a02c", "A12": "#fb9a99", "A13": "#e31a1c",
    "A14": "#fdbf6f",
}

# Annotate category
df["category"]     = df["ligand"].map(CATEGORY_MAP).fillna("Other")
df_ald["category"] = df_ald["ligand"].map(CATEGORY_MAP).fillna("Other")


# ── Global typography (Arial for everything, including SI figures) ───────────
matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["pdf.fonttype"] = 42  # embed Arial in the PDF


# ── Helpers ───────────────────────────────────────────────────────────────────
def lbl(ax, text, fontsize=11):
    """Bold panel label at top-left (used by main Figures 1 and 2)."""
    ax.text(-0.08, 1.06, text, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="bottom", ha="left",
            family="Arial")


def lbl_si(ax, text, fontsize=18):
    """SI panel label (Arial, bold, no brackets, enlarged)."""
    ax.text(-0.08, 1.06, text, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="bottom", ha="left",
            family="Arial")


# Fixed abbreviations for long substrate names, defined once and reused
# everywhere (axes labels, lollipop legends, dual-winner table). Substrates
# whose stripped name is already short enough fall through to the default.
SUBSTRATE_ABBR = {
    "31_(2E,4E,6E,8E)-3,7-dimethyl-9-(2,6,6-trimethylcyclohexen-1-yl)nona-2,4,6,8-tetraenal": "β-carotene ald.",
    "32_9-cis-Retinal":            "9-cis-retinal",
    "33_11-cis-Retinol":           "11-cis-retinol",
    "23_4-Hydroxynonenal":         "4-HNE",
    "34_4-hydroxy-2-nonenal":      "4-OH-2-nonenal",
    "36_3-(4-Isopropylphenyl)-2-methylpropanal": "iPr-Ph methylpropanal",
    "30_(E)-3-(3,4-Dihydroxyphenyl)prop-2-enal": "3,4-(OH)₂-cinnamald.",
    "38_(2R)-2-[(2S,3R,4S)-3,4-Dihydroxy-5-oxo -tetrahydrofuran-2-yl]-2-hydroxy-acetaldehyde": "DHF-OH-acetald.",
    "10_4-guanidinobutyraldehyde": "γ-guanidinobutyr-ald.",
    "48_alpha-aminoadipate semialdehyde": "α-AASA",
    "42_Glutamate-5-semialdehyde": "Glu-5-semiald.",
    "43_Glutamate-1-semialdehyde": "Glu-1-semiald.",
    "47_methylmalonate semialdehyde": "MeMal-semiald.",
    "49_2-Aminomuconic semialdehyde": "2-AM-semiald.",
    "50_10-formyltetrahydrofolate": "10-formyl-THF",
    "57_S-3-oxopropylglutathione": "S-3-oxopropyl-GSH",
    "58_S-(4-oxobutan-2-yl)glutathione": "S-4-oxobutyl-GSH",
    "60_S-Formylglutathione":      "S-formyl-GSH",
    "26_trans-2-hexadecenal":      "t-2-hexadecenal",
    "40_Furan-2-carbaldehyde":     "furan-2-carbald.",
    "41_Indole-3-acetaldehyde":    "indole-3-acetald.",
    "44_glutarate semialdehyde":   "glutarate-semiald.",
    "45_adipate semialdehyde":     "adipate-semiald.",
    "46_Succinic semialdehyde":    "succinic-semiald.",
    "29_Salicylaldehyde":          "salicylald.",
    "39_4-nitrobenzaldehyde":      "4-NO₂-benzald.",
    "37_thiobenzaldehyde":         "thiobenzald.",
    "27_Benzaldehyde":             "benzald.",
    "28_Phenylacetaldehyde":       "phenylacetald.",
    "35_Cinnamaldehyde":           "cinnamald.",
}


def short_name(ligand, maxlen=22):
    """Return a short, human-readable substrate label.

    Order of preference: fixed abbreviation > numeric-prefix-stripped name
    (truncated as a last resort)."""
    if ligand in SUBSTRATE_ABBR:
        return SUBSTRATE_ABBR[ligand]
    parts = ligand.split("_", 1)
    name = parts[1] if len(parts) > 1 else parts[0]
    return name if len(name) <= maxlen else name[:maxlen - 1] + "…"


def wilson_ci(k, n, z=1.96):
    """Wilson score 95% CI for a proportion k/n."""
    if n == 0:
        return 0.0, 0.0
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    margin = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return max(0, centre - margin), min(1, centre + margin)


def format_pval(p):
    if p < 0.001:
        return f"p = {p:.2e}"
    return f"p = {p:.3f}"


# ── Pre-compute stats ─────────────────────────────────────────────────────────
rho_cd, p_cd = spearmanr(
    df_ald["CYS_SG_to_aldC_distance_A"].dropna(),
    df_ald.loc[df_ald["CYS_SG_to_aldC_distance_A"].notna(), "DG_prediction_kcalmol"],
)
rho_conf, p_conf = spearmanr(
    df["rank1_confidence"].dropna(),
    df.loc[df["rank1_confidence"].notna(), "DG_prediction_kcalmol"],
)
mask_d = df_ald["CYS_SG_to_aldC_distance_A"].notna()
prod_dg  = df_ald.loc[mask_d & (df_ald["CYS_SG_to_aldC_distance_A"] <= 4.0), "DG_prediction_kcalmol"]
nprod_dg = df_ald.loc[mask_d & (df_ald["CYS_SG_to_aldC_distance_A"] >  4.0), "DG_prediction_kcalmol"]
stat_mwu, p_mwu = mannwhitneyu(prod_dg, nprod_dg, alternative="two-sided")

print("=== STATISTICS FOR REPORT ===")
print(f"Total pairs:             {len(df)}")
print(f"Aldehyde pairs:          {len(df_ald)}")
print(f"DG range:                {df.DG_prediction_kcalmol.min():.2f} to {df.DG_prediction_kcalmol.max():.2f} kcal/mol")
print(f"DG median:               {df.DG_prediction_kcalmol.median():.2f} kcal/mol")
print(f"Productive (<=4A):       {(df_ald.CYS_SG_to_aldC_distance_A<=4).sum()}/{len(df_ald)}")
print(f"rho(confidence,DG):      {rho_conf:.3f}  {format_pval(p_conf)}")
print(f"rho(distance,DG):        {rho_cd:.3f}  {format_pval(p_cd)}")
print(f"MWU prod vs nonprod:     prod median={prod_dg.median():.2f}, nonprod={nprod_dg.median():.2f}, {format_pval(p_mwu)}")
dual = df_ald[(df_ald.DG_prediction_kcalmol<=-7)&(df_ald.CYS_SG_to_aldC_distance_A<=4.0)]
print(f"Dual winners (DG<=-7, dist<=4): {len(dual)}")
print()


# ── Figure 1 ──────────────────────────────────────────────────────────────────
def make_fig1():
    sns.set_theme(style="ticks", font_scale=0.85)
    fig = plt.figure(figsize=(20, 15))
    gs  = GridSpec(3, 12, figure=fig,
                   height_ratios=[1.1, 1.4, 1.2],
                   hspace=0.52, wspace=0.55)

    ax_a = fig.add_subplot(gs[0, 0:4])
    ax_b = fig.add_subplot(gs[0, 4:8])
    ax_c = fig.add_subplot(gs[0, 8:12])
    ax_d = fig.add_subplot(gs[1, :])
    ax_e = fig.add_subplot(gs[2, 0:3])
    ax_f = fig.add_subplot(gs[2, 3:7])
    ax_g = fig.add_subplot(gs[2, 7:12])

    # ── Panel a: workflow schematic ──
    ax_a.set_xlim(0, 10); ax_a.set_ylim(0, 4.5); ax_a.axis("off")
    steps = [
        (1.0, 3.2, "Protein PDB\n+ Ligand SDF\n(420 pairs)", "#A5D6A7"),
        (4.0, 3.2, "DiffDock-L 1.1\n10 poses/pair\n2× A16 GPU",  "#90CAF9"),
        (7.0, 3.2, "PRODIGY-LIG\nΔG scoring",                     "#FFCC80"),
        (9.5, 3.2, "Cys-SG ↔\naldehyde-C\ndistance",              "#F48FB1"),
    ]
    for (x, y, text, col) in steps:
        box = FancyBboxPatch((x-0.9, y-0.7), 1.8, 1.55,
                              boxstyle="round,pad=0.1",
                              facecolor=col, edgecolor="#555", lw=1.2, zorder=2)
        ax_a.add_patch(box)
        ax_a.text(x, y+0.08, text, ha="center", va="center",
                  fontsize=7.5, fontweight="bold", zorder=3)
    for i in range(len(steps)-1):
        x0 = steps[i][0]+0.9; x1 = steps[i+1][0]-0.9
        ax_a.annotate("", xy=(x1, steps[i+1][1]+0.08), xytext=(x0, steps[i][1]+0.08),
                      arrowprops=dict(arrowstyle="->", color="#333", lw=1.8), zorder=4)
    ax_a.text(5.25, 0.6, "master_results.csv\n(ΔG + confidence + distance)", ha="center",
              fontsize=8, color="#555", style="italic",
              bbox=dict(boxstyle="round,pad=0.25", facecolor="#FFF9C4", edgecolor="#bbb", lw=0.8))
    ax_a.text(5.25, 0.15, "→ affinity map + geometry filter", ha="center", fontsize=7.5, color="#666")
    lbl(ax_a, "(a)")

    # ── Panel b: chemical space (RDKit MW vs logP) ──
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        mw_d, lp_d = {}, {}
        for lig in df["ligand"].unique():
            sdf = SDF_DIR / f"{lig}.sdf"
            if not sdf.exists(): continue
            try:
                mol = Chem.MolFromMolFile(str(sdf), sanitize=True)
                if mol:
                    mw_d[lig] = Descriptors.ExactMolWt(mol)
                    lp_d[lig] = Descriptors.MolLogP(mol)
            except Exception:
                pass
        prop_df = pd.DataFrame({"ligand": list(mw_d), "MW": list(mw_d.values()),
                                 "logP": [lp_d[k] for k in mw_d]})
        prop_df["category"] = prop_df["ligand"].map(CATEGORY_MAP).fillna("Other")
        for cat in CAT_ORDER:
            sub = prop_df[prop_df["category"] == cat]
            ax_b.scatter(sub["logP"], sub["MW"], c=CAT_COLORS[cat],
                         label=cat, s=40, alpha=0.85, edgecolors="white", lw=0.4, zorder=3)
        ax_b.set_xlabel("Calculated log P", fontsize=8)
        ax_b.set_ylabel("Exact MW (Da)", fontsize=8)
        ax_b.legend(loc="upper left", fontsize=6, framealpha=0.7,
                    handlelength=1.0, markerscale=0.8)
        rdkit_ok = True
    except ImportError:
        ax_b.text(0.5, 0.5, "RDKit not available\n(skipped)", transform=ax_b.transAxes,
                  ha="center", va="center", fontsize=10, color="grey")
        rdkit_ok = False
    ax_b.set_title("Chemical space of 59 substrates", fontsize=8, pad=4)
    lbl(ax_b, "(b)")

    # ── Panel c: DiffDock rank-1 confidence distributions ──
    order_c = sorted(PROTEINS, key=lambda p: df[df.protein==p].rank1_confidence.median())
    data_c  = [df[df.protein==p].rank1_confidence.dropna().values for p in order_c]
    parts   = ax_c.violinplot(data_c, positions=range(len(order_c)),
                               vert=False, widths=0.7, showmedians=True,
                               showextrema=True)
    for pc, prot in zip(parts["bodies"], order_c):
        pc.set_facecolor(PROT_COLORS[prot]); pc.set_alpha(0.75)
    parts["cmedians"].set_color("black"); parts["cmedians"].set_lw(1.5)
    for spine in ("cmins","cmaxes","cbars"):
        parts[spine].set_color("grey"); parts[spine].set_lw(0.8)
    ax_c.set_yticks(range(len(order_c))); ax_c.set_yticklabels(order_c, fontsize=8)
    ax_c.set_xlabel("DiffDock rank-1 confidence", fontsize=8)
    ax_c.set_title("Confidence per isoform", fontsize=8, pad=4)
    lbl(ax_c, "(c)")

    # ── Panel d: affinity heatmap ──
    pivot = df.pivot_table(index="protein", columns="ligand", values="DG_prediction_kcalmol")
    pivot = pivot.reindex(PROTEINS)
    filled = pivot.fillna(pivot.mean())
    col_link = linkage(pdist(filled.T.values, metric="euclidean"), method="ward")
    col_ord  = dendrogram(col_link, no_plot=True)["leaves"]
    pivot_s  = pivot.iloc[:, col_ord]

    im = ax_d.imshow(pivot_s.values, aspect="auto", cmap="RdYlBu_r",
                     vmin=-10.5, vmax=-5.0, interpolation="nearest")
    cb = fig.colorbar(im, ax=ax_d, fraction=0.015, pad=0.005,
                      label="ΔG$_{noelec}$ (kcal mol$^{-1}$)")
    cb.ax.tick_params(labelsize=7)

    ax_d.set_yticks(range(len(PROTEINS))); ax_d.set_yticklabels(PROTEINS, fontsize=8)
    # x-axis: substrate number colored by class
    lig_ord = list(pivot_s.columns)
    ax_d.set_xticks(range(len(lig_ord)))
    ax_d.set_xticklabels(
        [l.split("_")[0] for l in lig_ord],
        rotation=90, fontsize=5.5
    )
    for tick, lig in zip(ax_d.get_xticklabels(), lig_ord):
        tick.set_color(CAT_COLORS.get(CATEGORY_MAP.get(lig, "Other"), "#333"))
    ax_d.set_xlabel("Substrate (number, coloured by class; columns hierarchically clustered)",
                    fontsize=7.5)
    ax_d.set_ylabel("ALDH isoform", fontsize=8)

    # Mark top-1 per row with white star
    for i, prot in enumerate(PROTEINS):
        row = pivot_s.loc[prot]
        best_j = row.dropna().values.argmin() if row.dropna().size else None
        if best_j is not None:
            j_idx = list(pivot_s.columns).index(row.dropna().index[best_j])
            ax_d.text(j_idx, i, "★", ha="center", va="center",
                      fontsize=8, color="white", fontweight="bold")

    # Category legend
    legend_handles = [mpatches.Patch(facecolor=CAT_COLORS[c], label=c) for c in CAT_ORDER]
    ax_d.legend(handles=legend_handles, loc="lower right",
                bbox_to_anchor=(1.0, -0.4), ncol=7,
                fontsize=6.5, framealpha=0.8, handlelength=0.9)
    lbl(ax_d, "(d)")

    # ── Panel e: per-protein ΔG violin ──
    order_e = sorted(PROTEINS, key=lambda p: df_ald[df_ald.protein==p].DG_prediction_kcalmol.median())
    sns.violinplot(data=df_ald, x="DG_prediction_kcalmol", y="protein",
                   order=order_e, palette=PROT_COLORS, ax=ax_e,
                   inner="box", linewidth=0.8, scale="width", orient="h")
    ax_e.axvline(-7.0, ls="--", color="grey", lw=0.8, alpha=0.7)
    ax_e.set_xlabel("ΔG (kcal mol$^{-1}$)", fontsize=8)
    ax_e.set_ylabel(""); ax_e.tick_params(labelsize=8)
    ax_e.set_title("Affinity distributions\nper isoform", fontsize=8, pad=4)
    lbl(ax_e, "(e)")

    # ── Panel f: top-5 per protein lollipop ──
    top5 = (df_ald.groupby("protein", group_keys=False)
                  .apply(lambda g: g.nsmallest(5, "DG_prediction_kcalmol")))
    y_pos, yticks, ytick_labels = 0, [], []
    for prot in PROTEINS:
        sub = top5[top5.protein==prot].sort_values("DG_prediction_kcalmol")
        for _, row in sub.iterrows():
            ax_f.hlines(y_pos, row.DG_prediction_kcalmol, -5.1,
                        color=PROT_COLORS[prot], lw=1.5, alpha=0.7)
            ax_f.plot(row.DG_prediction_kcalmol, y_pos,
                      "o", color=PROT_COLORS[prot], ms=6, zorder=3)
            yticks.append(y_pos)
            ytick_labels.append(short_name(row.ligand, 20))
            y_pos += 1
        y_pos += 0.5  # gap between proteins
    ax_f.set_yticks(yticks)
    ax_f.set_yticklabels(ytick_labels, fontsize=5.5)
    ax_f.set_xlabel("ΔG (kcal mol$^{-1}$)", fontsize=8)
    ax_f.axvline(-7.0, ls="--", color="grey", lw=0.8, alpha=0.6)
    ax_f.set_title("Top-5 substrates per isoform", fontsize=8, pad=4)
    prot_handles = [mpatches.Patch(facecolor=PROT_COLORS[p], label=p) for p in PROTEINS]
    ax_f.legend(handles=prot_handles, loc="lower right",
                fontsize=6, handlelength=0.8, framealpha=0.7)
    lbl(ax_f, "(f)")

    # ── Panel g: substrate selectivity bar ──
    thr = -7.0
    sel = (df_ald[df_ald.DG_prediction_kcalmol <= thr]
           .groupby("ligand")["protein"].nunique().reset_index())
    sel.columns = ["ligand", "n_prots"]
    all_ligs = pd.DataFrame({"ligand": df_ald["ligand"].unique()})
    sel = all_ligs.merge(sel, on="ligand", how="left").fillna(0)
    sel["lig_num"] = sel["ligand"].str.split("_", n=1).str[0].astype(int)
    sel = sel.sort_values("lig_num")
    colors_g = [CAT_COLORS.get(CATEGORY_MAP.get(l, "Other"), "#9E9E9E")
                for l in sel["ligand"]]
    ax_g.bar(range(len(sel)), sel["n_prots"], color=colors_g, width=0.8, linewidth=0)
    ax_g.axhline(4, ls="--", color="grey", lw=0.8, alpha=0.7)
    ax_g.set_xticks([])
    ax_g.set_ylabel("No. of isoforms\nbinding ≤ −7 kcal mol$^{-1}$", fontsize=7.5)
    ax_g.set_xlabel("Substrate (sorted by ID)", fontsize=7.5)
    ax_g.set_title("Pan-ALDH binding at ΔG ≤ −7 kcal mol$^{-1}$", fontsize=8, pad=4)
    for i, row in sel[sel["n_prots"] >= 4].iterrows():
        idx = sel.index.get_loc(i)
        ax_g.text(idx, row.n_prots + 0.05, row["ligand"].split("_")[0],
                  ha="center", va="bottom", fontsize=5.5, rotation=90)
    lbl(ax_g, "(g)")

    fig.suptitle("Systematic blind docking of 413 ALDH–aldehyde pairs reveals "
                 "isoform-specific affinity landscapes",
                 fontsize=10, y=0.98, fontweight="bold")
    fig.savefig(FIGS / "fig1_landscape.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(FIGS / "fig1_landscape.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("✓ Figure 1 saved")


# ── Figure 2 ──────────────────────────────────────────────────────────────────
def make_fig2():
    sns.set_theme(style="ticks", font_scale=0.85)
    fig = plt.figure(figsize=(20, 17))
    gs  = GridSpec(5, 12, figure=fig,
                   height_ratios=[1.1, 1.5, 1.1, 1.3, 1.0],
                   hspace=0.55, wspace=0.55)

    ax_a  = fig.add_subplot(gs[0, 0:4])
    ax_b  = fig.add_subplot(gs[0, 4:8])
    ax_c  = fig.add_subplot(gs[0, 8:12])
    ax_c1 = fig.add_subplot(gs[1, 0:6])
    ax_c2 = fig.add_subplot(gs[1, 6:12])
    ax_d  = fig.add_subplot(gs[2, 0:3])
    ax_e  = fig.add_subplot(gs[2, 3:6])
    ax_f  = fig.add_subplot(gs[2, 6:9])
    ax_g  = fig.add_subplot(gs[2, 9:12])
    ax_h  = fig.add_subplot(gs[3, :])
    ax_i  = fig.add_subplot(gs[4, :])      # panel (i) now spans full row
    # Panel (j) class composition has been relocated to Supplementary Fig S5

    # ── Panels c1 / c2: structural single-pair comparison ──
    import matplotlib.image as mpimg
    POSE_DIR = FIGS / "pose_panel"
    a8_img  = mpimg.imread(POSE_DIR / "panel_A8_productive.png")
    a12_img = mpimg.imread(POSE_DIR / "panel_A12_nonproductive.png")
    for ax, img, hdr, dist, verdict, vcol in [
        (ax_c1, a8_img,  r"A8 $\times$ $\beta$-carotene ald.",  2.27, "productive",     "#2E7D32"),
        (ax_c2, a12_img, r"A12 $\times$ $\beta$-carotene ald.", 5.49, "non-productive", "#C62828"),
    ]:
        ax.imshow(img); ax.axis("off")
        ax.set_title(hdr, fontsize=9, pad=4, fontweight="bold")
        ax.text(0.03, 0.97, f"Cys-SG ↔ aldehyde-C: {dist:.2f} Å",
                transform=ax.transAxes, fontsize=8, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="grey", lw=0.6, alpha=0.9))
        ax.text(0.97, 0.97, verdict,
                transform=ax.transAxes, fontsize=8, fontweight="bold",
                color="white", va="top", ha="right",
                bbox=dict(boxstyle="round,pad=0.28", fc=vcol, ec="none", alpha=0.95))
    lbl(ax_c1, "(c1)")
    lbl(ax_c2, "(c2)")

    dist_vals = df_ald["CYS_SG_to_aldC_distance_A"].dropna()
    n_prod = (dist_vals <= 4.0).sum()
    n_tot  = len(dist_vals)

    # ── Panel a: distance histogram ──
    ax_a.hist(dist_vals, bins=40, color="#5C85D6", edgecolor="white", lw=0.3, alpha=0.85)
    ax_a.axvline(4.0, color="#D62728", ls="--", lw=1.8, label="4 Å threshold")
    ax_a.fill_betweenx([0, ax_a.get_ylim()[1] if ax_a.get_ylim()[1]>0 else 50],
                        0, 4.0, color="#A8D5A2", alpha=0.15, zorder=0)
    ax_a.set_xlabel("Cys-SG ↔ aldehyde-C distance (Å)", fontsize=8)
    ax_a.set_ylabel("Number of complexes", fontsize=8)
    ax_a.set_title("Distance distribution\n(all aldehyde pairs)", fontsize=8, pad=4)
    ax_a.text(0.97, 0.95, f"{n_prod}/{n_tot}\nproductive\n(≤ 4 Å)",
              transform=ax_a.transAxes, ha="right", va="top",
              fontsize=7.5, color="#D62728",
              bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#D62728", lw=0.8))
    ax_a.legend(fontsize=7.5, loc="upper right", framealpha=0.7)
    lbl(ax_a, "(a)")

    # ── Panel b: per-protein distance violin ──
    order_b = sorted(PROTEINS,
                     key=lambda p: df_ald[df_ald.protein==p]
                                       .CYS_SG_to_aldC_distance_A.median())
    ax_b.axhspan(0, 4, alpha=0.12, color="#A8D5A2", zorder=0)
    sns.violinplot(data=df_ald, x="protein", y="CYS_SG_to_aldC_distance_A",
                   order=order_b, palette=PROT_COLORS, ax=ax_b,
                   inner="box", linewidth=0.8, scale="width")
    ax_b.axhline(4.0, color="#D62728", ls="--", lw=1.2, alpha=0.8, label="4 Å")
    ax_b.set_xlabel("ALDH isoform", fontsize=8)
    ax_b.set_ylabel("Distance (Å)", fontsize=8)
    ax_b.set_title("Catalytic-Cys distance\nper isoform", fontsize=8, pad=4)
    ax_b.legend(fontsize=7, framealpha=0.7)
    lbl(ax_b, "(b)")

    # ── Panel c: ΔG vs distance scatter ──
    mask_c = df_ald["CYS_SG_to_aldC_distance_A"].notna()
    sc = ax_c.scatter(
        df_ald.loc[mask_c, "CYS_SG_to_aldC_distance_A"],
        df_ald.loc[mask_c, "DG_prediction_kcalmol"],
        c=df_ald.loc[mask_c, "rank1_confidence"].fillna(-1),
        cmap="viridis", s=18, alpha=0.75, edgecolor="none",
    )
    ax_c.axvline(4.0, color="grey", ls="--", lw=0.8, alpha=0.8)
    cb_c = fig.colorbar(sc, ax=ax_c, fraction=0.04, pad=0.03)
    cb_c.set_label("DiffDock confidence", fontsize=7)
    cb_c.ax.tick_params(labelsize=6)
    rho_c, p_c = spearmanr(
        df_ald.loc[mask_c, "CYS_SG_to_aldC_distance_A"],
        df_ald.loc[mask_c, "DG_prediction_kcalmol"],
    )
    ax_c.text(0.97, 0.97,
              f"Spearman ρ = {rho_c:.2f}\n{format_pval(p_c)}",
              transform=ax_c.transAxes, ha="right", va="top",
              fontsize=7,
              bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="grey", lw=0.6))
    ax_c.set_xlabel("Cys-SG distance (Å)", fontsize=8)
    ax_c.set_ylabel("ΔG (kcal mol$^{-1}$)", fontsize=8)
    ax_c.set_title("Affinity vs catalytic geometry", fontsize=8, pad=4)
    lbl(ax_c, "(c)")

    # ── Panel d: four-quadrant scatter ──
    XTHRESH, YTHRESH = 4.0, -7.0
    QUAD_COLS = {
        "Productive\nhigh-affinity":     "#FFD700",
        "High-affinity\nnon-productive":  "#E07B39",
        "Productive\nlow-affinity":       "#2196F3",
        "Neither":                        "#BDBDBD",
    }

    def get_quad(row):
        if pd.isna(row["CYS_SG_to_aldC_distance_A"]): return "Neither"
        prod = row["CYS_SG_to_aldC_distance_A"] <= XTHRESH
        strg = row["DG_prediction_kcalmol"]       <= YTHRESH
        if prod and strg:  return "Productive\nhigh-affinity"
        if not prod and strg: return "High-affinity\nnon-productive"
        if prod and not strg: return "Productive\nlow-affinity"
        return "Neither"

    df_ald["quad"] = df_ald.apply(get_quad, axis=1)
    for q, col in QUAD_COLS.items():
        sub = df_ald[df_ald["quad"] == q]
        ax_d.scatter(sub["CYS_SG_to_aldC_distance_A"], sub["DG_prediction_kcalmol"],
                     c=col, label=f"{q} (n={len(sub)})", s=14, alpha=0.72,
                     edgecolor="white", lw=0.2, zorder=3)
    ax_d.axvline(XTHRESH, color="grey", ls="--", lw=0.8)
    ax_d.axhline(YTHRESH, color="grey", ls="--", lw=0.8)
    ax_d.set_xlabel("Distance (Å)", fontsize=8)
    ax_d.set_ylabel("ΔG (kcal mol$^{-1}$)", fontsize=8)
    ax_d.set_title("Two-metric quadrant\nanalysis", fontsize=8, pad=4)
    ax_d.legend(loc="lower right", fontsize=5.5, framealpha=0.8, handlelength=0.9)
    lbl(ax_d, "(d)")

    # ── Panel e: productive fraction per protein ──
    bars_e, errs_lo, errs_hi = [], [], []
    for prot in PROTEINS:
        sub = df_ald[df_ald.protein==prot]["CYS_SG_to_aldC_distance_A"].dropna()
        k   = (sub <= 4.0).sum(); n = len(sub)
        frac = k / n if n else 0
        lo, hi = wilson_ci(k, n)
        bars_e.append(frac); errs_lo.append(frac - lo); errs_hi.append(hi - frac)
    x_e = np.arange(len(PROTEINS))
    ax_e.bar(x_e, bars_e, color=[PROT_COLORS[p] for p in PROTEINS],
             yerr=[errs_lo, errs_hi], capsize=3, width=0.65,
             error_kw={"lw": 1.2, "capthick": 1.2})
    ax_e.axhline(n_prod / n_tot, color="black", ls="--", lw=0.9, alpha=0.6,
                 label=f"Overall ({n_prod/n_tot:.0%})")
    ax_e.set_xticks(x_e); ax_e.set_xticklabels(PROTEINS, fontsize=8)
    ax_e.set_ylabel("Productive fraction", fontsize=8)
    ax_e.set_ylim(0, 1)
    ax_e.set_title("Productive poses\nper isoform (95% CI)", fontsize=8, pad=4)
    ax_e.legend(fontsize=7, framealpha=0.7)
    lbl(ax_e, "(e)")

    # ── Panel f: productive fraction per substrate class ──
    prod_frac = {}
    for cat in CAT_ORDER:
        sub = df_ald[df_ald.category==cat]["CYS_SG_to_aldC_distance_A"].dropna()
        if len(sub) == 0: prod_frac[cat] = 0; continue
        prod_frac[cat] = (sub <= 4.0).sum() / len(sub)
    cats_sorted = sorted(CAT_ORDER[:-1], key=lambda c: -prod_frac[c])
    colors_f = [CAT_COLORS[c] for c in cats_sorted]
    fracs_f  = [prod_frac[c]  for c in cats_sorted]
    bars_f   = ax_f.barh(cats_sorted, fracs_f, color=colors_f, height=0.6)
    ax_f.axvline(n_prod / n_tot, color="black", ls="--", lw=0.9, alpha=0.6)
    for bar, frac in zip(bars_f, fracs_f):
        ax_f.text(frac + 0.01, bar.get_y() + bar.get_height()/2,
                  f"{frac:.0%}", va="center", ha="left", fontsize=7.5)
    ax_f.set_xlim(0, 1.05)
    ax_f.set_xlabel("Productive fraction (≤ 4 Å)", fontsize=8)
    ax_f.set_title("Productive poses per\nsubstrate class", fontsize=8, pad=4)
    lbl(ax_f, "(f)")

    # ── Panel g: ΔG violin productive vs non-productive ──
    df_vg = df_ald[mask_c].copy()
    df_vg["pose"] = np.where(df_vg["CYS_SG_to_aldC_distance_A"] <= 4, "Productive\n(≤ 4 Å)", "Non-productive\n(> 4 Å)")
    sns.violinplot(data=df_vg, x="pose", y="DG_prediction_kcalmol",
                   palette={"Productive\n(≤ 4 Å)": "#5C85D6", "Non-productive\n(> 4 Å)": "#D62728"},
                   ax=ax_g, inner="box", linewidth=0.8)
    stat_g, p_g = mannwhitneyu(prod_dg, nprod_dg, alternative="two-sided")
    ax_g.text(0.5, 0.02, f"Mann–Whitney\n{format_pval(p_g)}",
              transform=ax_g.transAxes, ha="center", va="bottom", fontsize=7,
              bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="grey", lw=0.6))
    ax_g.set_xlabel("")
    ax_g.set_ylabel("ΔG (kcal mol$^{-1}$)", fontsize=8)
    ax_g.set_title("ΔG: productive vs\nnon-productive poses", fontsize=8, pad=4)
    ax_g.tick_params(axis="x", labelsize=7.5)
    lbl(ax_g, "(g)")

    # ── Panel h: dual-winner lollipop ──
    dw = (df_ald[(df_ald.DG_prediction_kcalmol <= -7.0) &
                 (df_ald.CYS_SG_to_aldC_distance_A <= 4.0)]
             .nsmallest(15, "DG_prediction_kcalmol"))
    dw = dw.sort_values("DG_prediction_kcalmol", ascending=False)
    y_h = np.arange(len(dw))
    ax_h.hlines(y_h, dw["DG_prediction_kcalmol"].values, -6.8,
                color=[PROT_COLORS.get(p, "grey") for p in dw["protein"]],
                lw=2.0, alpha=0.8)
    sc_h = ax_h.scatter(dw["DG_prediction_kcalmol"], y_h,
                         s=120 / (dw["CYS_SG_to_aldC_distance_A"] + 0.5),
                         c=[PROT_COLORS.get(p, "grey") for p in dw["protein"]],
                         edgecolors="black", lw=0.5, zorder=4)
    ax_h.set_yticks(y_h)
    ax_h.set_yticklabels(
        [f"{row.protein}  {short_name(row.ligand, 30)}" for _, row in dw.iterrows()],
        fontsize=7.5
    )
    ax_h.set_xlabel("ΔG (kcal mol$^{-1}$)", fontsize=8)
    ax_h.set_title("Top-15 dual-winner pairs (ΔG ≤ −7 kcal mol⁻¹ AND Cys-SG ≤ 4 Å)\n"
                   "Dot size ∝ proximity (1/distance)",
                   fontsize=8, pad=4)
    ax_h.axvline(-7.0, ls="--", color="grey", lw=0.8, alpha=0.7)
    ax_h.invert_xaxis()
    prot_handles2 = [Line2D([0],[0], marker="o", color="w",
                             markerfacecolor=PROT_COLORS[p], ms=8, label=p) for p in PROTEINS]
    ax_h.legend(handles=prot_handles2, loc="lower right",
                fontsize=7, ncol=7, framealpha=0.8)
    lbl(ax_h, "(h)")

    # ── Panel i: three-gate filter funnel ──
    n_start = len(df_ald)
    n_aff   = int(((df_ald["DG_prediction_kcalmol"] <= -7.0)).sum())
    n_dual  = int(((df_ald["DG_prediction_kcalmol"] <= -7.0) &
                   (df_ald["CYS_SG_to_aldC_distance_A"] <= 4.0)).sum())
    # Ligand efficiency uses the pre-computed per-ligand heavy-atom counts
    # from results/interaction_features.csv (produced by 09_interaction_analysis.py).
    feat = pd.read_csv(RESULTS / "interaction_features.csv")
    heavy_map = dict(zip(feat["ligand"], feat["HeavyAtoms"]))
    df_le = df_ald.copy()
    df_le["heavy"] = df_le["ligand"].map(heavy_map)
    df_le["LE"]    = df_le["DG_prediction_kcalmol"].abs() / df_le["heavy"]
    le_med = df_le["LE"].median()
    n_triple = int(((df_le["DG_prediction_kcalmol"] <= -7.0) &
                    (df_le["CYS_SG_to_aldC_distance_A"] <= 4.0) &
                    (df_le["LE"] >= le_med)).sum())
    stages = ["All aldehyde\npairs",
              f"ΔG ≤ −7\nkcal mol$^{{-1}}$",
              "+ Cys-SG\n≤ 4 Å",
              f"+ LE ≥ {le_med:.2f}\nkcal mol$^{{-1}}$\natom$^{{-1}}$"]
    counts = [n_start, n_aff, n_dual, n_triple]
    cols_funnel = ["#9E9E9E", "#5C85D6", "#FFB300", "#2E7D32"]
    y_pos = np.arange(len(stages))[::-1]
    ax_i.barh(y_pos, counts, color=cols_funnel, edgecolor="black", lw=0.6, height=0.6)
    for y, n, prev in zip(y_pos, counts, [None] + counts[:-1]):
        retain = "" if prev is None else f"  ({100*n/prev:.0f}% retained)"
        ax_i.text(n + 8, y, f"n = {n}{retain}", va="center", fontsize=8)
    ax_i.set_yticks(y_pos); ax_i.set_yticklabels(stages, fontsize=8)
    ax_i.set_xlabel("Pairs retained", fontsize=8)
    ax_i.set_xlim(0, n_start * 1.35)
    ax_i.set_title("Three-gate filter funnel: affinity × geometry × ligand efficiency",
                   fontsize=8, pad=4)
    ax_i.invert_yaxis()
    lbl(ax_i, "(i)")

    # Panel (j) relocated to Supplementary Fig S5 (see make_figS5).

    fig.suptitle("Catalytic-cysteine geometry stratifies docking poses "
                 "and identifies productive substrate–enzyme pairs",
                 fontsize=10, y=0.98, fontweight="bold")
    fig.savefig(FIGS / "fig2_geometry.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(FIGS / "fig2_geometry.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("✓ Figure 2 saved")


# ── SI Figure S1 ──────────────────────────────────────────────────────────────
def make_figS1():
    sns.set_theme(style="ticks", font_scale=1.0)
    fig = plt.figure(figsize=(16, 5.5))
    gs  = GridSpec(1, 3, figure=fig, wspace=0.40,
                   bottom=0.28, top=0.92, left=0.05, right=0.84)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    # Panel (d) "Binding affinity by substrate class" removed because S3(b)
    # already shows the same data as a violin.

    # ── Panel a: PRODIGY return code success ──
    rc_counts = (df.groupby(["protein", "prodigy_rc"])
                   .size().reset_index(name="n"))
    success = (rc_counts[rc_counts["prodigy_rc"] == 0]
                .set_index("protein")["n"].reindex(PROTEINS, fill_value=0))
    total   = df.groupby("protein").size().reindex(PROTEINS, fill_value=0)
    fail    = total - success
    x_a = np.arange(len(PROTEINS))
    ax_a.bar(x_a, success.values / total.values * 100,
             color="#4CAF50", label="Success (rc=0)", width=0.6)
    ax_a.bar(x_a, fail.values / total.values * 100,
             bottom=success.values / total.values * 100,
             color="#F44336", label="Failure (rc≠0)", width=0.6)
    ax_a.set_xticks(x_a); ax_a.set_xticklabels(PROTEINS, fontsize=12)
    ax_a.set_ylabel("PRODIGY-LIG run outcome (%)", fontsize=13)
    ax_a.set_title("PRODIGY-LIG pipeline success rates", fontsize=13, pad=6)
    ax_a.set_ylim(0, 110)
    ax_a.tick_params(axis="y", labelsize=11)
    ax_a.legend(fontsize=10, loc="upper center", ncol=2,
                bbox_to_anchor=(0.5, -0.18), frameon=False)
    lbl_si(ax_a, "a")

    # ── Panel b: SMARTS aldehyde detection ──
    ald_counts = df.groupby(["protein", "is_aldehyde"]).size().unstack(fill_value=0)
    ald_t = ald_counts.get(True,  pd.Series(0, index=PROTEINS))
    ald_f = ald_counts.get(False, pd.Series(0, index=PROTEINS))
    x_b   = np.arange(len(PROTEINS))
    ax_b.bar(x_b, ald_t.reindex(PROTEINS, fill_value=0),
             color="#5C85D6", label="Aldehyde (SMARTS ✓)", width=0.6)
    ax_b.bar(x_b, ald_f.reindex(PROTEINS, fill_value=0),
             bottom=ald_t.reindex(PROTEINS, fill_value=0),
             color="#FF8C00", label="Non-aldehyde", width=0.6)
    ax_b.legend(fontsize=10, loc="upper center", ncol=2,
                bbox_to_anchor=(0.5, -0.18), frameon=False)
    ax_b.set_xticks(x_b); ax_b.set_xticklabels(PROTEINS, rotation=0, fontsize=12)
    ax_b.set_ylabel("Number of substrates", fontsize=13)
    ax_b.set_title("SMARTS aldehyde identification per isoform",
                   fontsize=13, pad=6)
    ax_b.set_xlabel("")
    ax_b.tick_params(axis="y", labelsize=11)
    lbl_si(ax_b, "b")

    # ── Panel c: confidence vs ΔG by substrate class ──
    mask_c = df["rank1_confidence"].notna()
    for cat in CAT_ORDER:
        sub = df[mask_c & (df["category"] == cat)]
        if len(sub) < 3: continue
        rho_s, p_s = spearmanr(sub["rank1_confidence"], sub["DG_prediction_kcalmol"])
        ax_c.scatter(sub["rank1_confidence"], sub["DG_prediction_kcalmol"],
                     c=CAT_COLORS[cat], label=f"{cat} (ρ={rho_s:.2f})",
                     s=22, alpha=0.7, edgecolor="none")
    ax_c.set_xlabel("DiffDock rank-1 confidence", fontsize=13)
    ax_c.set_ylabel("ΔG (kcal mol$^{-1}$)", fontsize=13)
    ax_c.set_title("Confidence–affinity correlation by substrate class",
                   fontsize=13, pad=6)
    ax_c.tick_params(labelsize=11)
    ax_c.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
                ncol=1, fontsize=9, frameon=False, handlelength=1.0)
    lbl_si(ax_c, "c")

    fig.savefig(FIGS / "figS1_qc.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(FIGS / "figS1_qc.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("✓ SI Figure S1 saved")


# ── SI Figure S2 ──────────────────────────────────────────────────────────────
def make_figS2():
    sns.set_theme(style="ticks", font_scale=1.0)
    fig, axes = plt.subplots(4, 2, figsize=(13, 16))
    axes = axes.flatten()

    for idx, prot in enumerate(PROTEINS):
        ax = axes[idx]
        sub = df_ald[df_ald.protein == prot].sort_values("DG_prediction_kcalmol")
        colors_s = [CAT_COLORS.get(CATEGORY_MAP.get(l, "Other"), "#9E9E9E")
                    for l in sub["ligand"]]
        ax.barh(range(len(sub)), sub["DG_prediction_kcalmol"],
                color=colors_s, height=0.75, linewidth=0)

        # Mark productive poses (dist <= 4 Å) with a gold star
        for j, (_, row) in enumerate(sub.iterrows()):
            if pd.notna(row["CYS_SG_to_aldC_distance_A"]) and row["CYS_SG_to_aldC_distance_A"] <= 4.0:
                ax.text(row["DG_prediction_kcalmol"] - 0.04, j, "★",
                        ha="right", va="center", fontsize=6,
                        color="#FFD700")

        # Label the single top substrate in the top-left corner of the panel
        if len(sub):
            top = sub.iloc[0]
            ax.text(0.03, 0.97, f"top: {short_name(top['ligand'], 24)}",
                    transform=ax.transAxes,
                    ha="left", va="top", fontsize=10, color="#222")

        ax.axvline(-7.0, ls="--", color="grey", lw=0.8, alpha=0.7)
        ax.set_yticks([])
        ax.set_xlabel("ΔG (kcal mol$^{-1}$)", fontsize=13)
        ax.tick_params(axis="x", labelsize=11)
        ax.set_title(
            f"{prot}  (median ΔG = {sub.DG_prediction_kcalmol.median():.2f})",
            fontsize=13, pad=6, color=PROT_COLORS[prot],
        )
        lbl_si(ax, "abcdefg"[idx])

    # ── Panel h: inter-isoform Pearson correlation (moved from Fig. 1) ──
    ax_h = axes[7]
    piv_prot = df.pivot_table(index="protein", columns="ligand",
                               values="DG_prediction_kcalmol").fillna(
                    df.pivot_table(index="protein", columns="ligand",
                                   values="DG_prediction_kcalmol").mean())
    corr = piv_prot.T.corr(method="pearson").loc[PROTEINS, PROTEINS]
    sns.heatmap(corr, cmap="RdYlBu_r", center=0, vmin=0, vmax=1,
                annot=True, fmt=".2f", annot_kws={"size": 11},
                ax=ax_h, cbar_kws={"label": "Pearson r", "shrink": 0.8},
                square=True, linewidths=0.4, linecolor="white")
    ax_h.tick_params(labelsize=11)
    ax_h.set_title("Inter-isoform correlation",
                   fontsize=13, pad=6, color="#333")
    lbl_si(ax_h, "h")

    # Category legend for the whole figure
    handles_leg = [mpatches.Patch(facecolor=CAT_COLORS[c], label=c) for c in CAT_ORDER]
    star_handle = Line2D([0],[0], marker="*", color="w",
                          markerfacecolor="#FFD700", ms=12, label="★ productive (≤ 4 Å)")
    fig.legend(handles=handles_leg + [star_handle],
               loc="lower center", ncol=4, fontsize=11, framealpha=0.85,
               bbox_to_anchor=(0.5, 0.005))
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(FIGS / "figS2_isoforms.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(FIGS / "figS2_isoforms.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("✓ SI Figure S2 saved")


# ── SI Figure S3 ──────────────────────────────────────────────────────────────
def make_figS3():
    sns.set_theme(style="ticks", font_scale=1.0)
    fig = plt.figure(figsize=(14, 10))
    gs  = GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.45)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1], polar=True)

    # ── Panel a: Tanimoto similarity dendrogram ──
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import DataStructs

        fps_list, fp_names, fp_cats = [], [], []
        for lig in df["ligand"].unique():
            sdf = SDF_DIR / f"{lig}.sdf"
            if not sdf.exists(): continue
            try:
                mol = Chem.MolFromMolFile(str(sdf), sanitize=True)
                if mol:
                    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
                    fps_list.append(fp)
                    fp_names.append(lig.split("_")[0])
                    fp_cats.append(CATEGORY_MAP.get(lig, "Other"))
            except Exception:
                pass

        n = len(fps_list)
        if n > 4:
            dist_mat = np.zeros((n, n))
            for i in range(n):
                sims = DataStructs.BulkTanimotoSimilarity(fps_list[i], fps_list)
                dist_mat[i] = 1 - np.array(sims)
            condensed = squareform(dist_mat)
            Z = linkage(condensed, method="ward")
            dend = dendrogram(Z, labels=fp_names,
                              color_threshold=0, above_threshold_color="grey",
                              ax=ax_a, leaf_rotation=90, leaf_font_size=9,
                              no_labels=False)
            # Colour leaf labels by category
            leaf_order = dend["leaves"]
            cats_ord   = [fp_cats[i] for i in leaf_order]
            for tick, cat in zip(ax_a.get_xticklabels(), cats_ord):
                tick.set_color(CAT_COLORS.get(cat, "#333"))
            ax_a.set_xlabel("Substrate number (colour = class)", fontsize=13)
            ax_a.set_ylabel("Ward linkage distance", fontsize=13)
            ax_a.set_title("Tanimoto similarity dendrogram (Morgan FP, r=2)",
                           fontsize=13, pad=6)
            ax_a.tick_params(axis="y", labelsize=11)
    except Exception as e:
        ax_a.text(0.5, 0.5, f"Skipped\n({e})", transform=ax_a.transAxes,
                  ha="center", va="center", fontsize=11, color="grey")
    lbl_si(ax_a, "a")

    # ── Panel b: ΔG by substrate class ──
    cats_b = CAT_ORDER[:-1]
    df_b   = df_ald[df_ald["category"].isin(cats_b)].copy()
    sns.violinplot(data=df_b, x="category", y="DG_prediction_kcalmol",
                   order=cats_b, palette=CAT_COLORS,
                   ax=ax_b, inner="box", linewidth=0.9)
    groups_b = [df_ald[df_ald.category==c].DG_prediction_kcalmol.dropna().values
                for c in cats_b if len(df_ald[df_ald.category==c]) >= 2]
    if groups_b:
        stat_b, p_b = kruskal(*groups_b)
        ax_b.text(0.97, 0.97, f"Kruskal–Wallis\n{format_pval(p_b)}",
                  transform=ax_b.transAxes, ha="right", va="top", fontsize=11,
                  bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                            edgecolor="grey", lw=0.7))
    ax_b.set_xticklabels(ax_b.get_xticklabels(), rotation=30, ha="right",
                          fontsize=11)
    ax_b.tick_params(axis="y", labelsize=11)
    ax_b.set_xlabel("")
    ax_b.set_ylabel("ΔG (kcal mol$^{-1}$)", fontsize=13)
    ax_b.set_title("Affinity distribution by substrate class",
                   fontsize=13, pad=6)
    lbl_si(ax_b, "b")

    # ── Panel c: median distance by class per protein ──
    med_dist = (df_ald.groupby(["protein", "category"])
                      ["CYS_SG_to_aldC_distance_A"].median().reset_index())
    for prot in PROTEINS:
        sub_c = med_dist[med_dist.protein==prot]
        cat_order_c = [c for c in CAT_ORDER if c in sub_c.category.values]
        x_vals = [cat_order_c.index(c) if c in cat_order_c else np.nan
                  for c in sub_c.category]
        ax_c.plot(sub_c.category.values,
                  sub_c.CYS_SG_to_aldC_distance_A.values,
                  "o-", color=PROT_COLORS[prot], label=prot,
                  ms=7, lw=1.6, alpha=0.85)
    ax_c.axhline(4.0, color="#D62728", ls="--", lw=1.0, alpha=0.7, label="4 Å")
    ax_c.set_xticklabels(ax_c.get_xticklabels(), rotation=30, ha="right",
                          fontsize=11)
    ax_c.tick_params(axis="y", labelsize=11)
    ax_c.set_ylabel("Median Cys-SG distance (Å)", fontsize=13)
    ax_c.set_title("Catalytic geometry by substrate class and isoform",
                   fontsize=13, pad=6)
    ax_c.legend(fontsize=9, loc="upper center", ncol=4, frameon=False,
                bbox_to_anchor=(0.5, -0.28))
    lbl_si(ax_c, "c")

    # ── Panel d: radar chart — isoform selectivity across classes ──
    cats_r = CAT_ORDER[:-1]  # 6 categories
    N = len(cats_r)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    for prot in PROTEINS:
        vals = []
        for cat in cats_r:
            sub_r = df_ald[(df_ald.protein==prot) & (df_ald.category==cat)]
            if len(sub_r) == 0:
                vals.append(0)
            else:
                vals.append((sub_r.DG_prediction_kcalmol <= -7.0).sum() / len(sub_r))
        vals += vals[:1]
        ax_d.plot(angles, vals, color=PROT_COLORS[prot], lw=1.5, label=prot, alpha=0.85)
        ax_d.fill(angles, vals, color=PROT_COLORS[prot], alpha=0.08)

    ax_d.set_xticks(angles[:-1])
    ax_d.set_xticklabels(cats_r, size=11)
    ax_d.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax_d.set_yticklabels(["25%", "50%", "75%", "100%"], size=10)
    ax_d.set_title("Selectivity profile (fraction with ΔG ≤ −7 kcal mol⁻¹)",
                   fontsize=13, pad=34)
    ax_d.legend(loc="center left", bbox_to_anchor=(1.18, 0.5),
                fontsize=10, framealpha=0.85)
    lbl_si(ax_d, "d")

    fig.savefig(FIGS / "figS3_chemspace.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(FIGS / "figS3_chemspace.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("✓ SI Figure S3 saved")


# ── SI Table S1 ───────────────────────────────────────────────────────────────
def make_table_S1():
    dw = (df_ald[(df_ald.DG_prediction_kcalmol <= -7.0) &
                 (df_ald.CYS_SG_to_aldC_distance_A <= 4.0)]
             .sort_values("DG_prediction_kcalmol")
             .head(20)
             .copy())
    dw["substrate"] = dw["ligand"].apply(lambda l: short_name(l, 35).replace("_", " "))
    dw["category"]  = dw["ligand"].map(CATEGORY_MAP).fillna("Other")
    dw["near"] = dw["CYS_SG_to_aldC_distance_A"] <= 2.5

    lines = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\small",
        # Title above the table
        r"{\noindent\textbf{Table S1: Top-20 dual-winner protein--substrate pairs.}\par}",
        r"\vspace{0.4em}",
        r"\begin{tabular}{@{}llrrrl@{}}",
        r"\toprule",
        r"\textbf{Protein} & \textbf{Substrate} & "
        r"$\boldsymbol{\Delta G}$ & \textbf{Dist.} & "
        r"\textbf{Confidence} & \textbf{Class} \\",
        r" & & \textbf{(kcal mol$^{-1}$)} & \textbf{(\AA)} & & \\",
        r"\midrule",
    ]
    def tex_safe(s):
        return (str(s)
                .replace("α", r"$\alpha$")
                .replace("β", r"$\beta$")
                .replace("→", r"$\rightarrow$")
                .replace("↔", r"$\leftrightarrow$")
                .replace("★", r"$\star$")
                .replace("≤", r"$\leq$")
                .replace("≥", r"$\geq$"))

    for _, row in dw.iterrows():
        conf_str = f"{row.rank1_confidence:.2f}" if pd.notna(row.rank1_confidence) else "---"
        dist_str = f"{row.CYS_SG_to_aldC_distance_A:.2f}"
        cells = [tex_safe(row.protein), tex_safe(row.substrate),
                 f"{row.DG_prediction_kcalmol:.2f}", dist_str,
                 conf_str, tex_safe(row.category)]
        shade = r"\rowcolor{nearattack}" if row.near else r"\rowcolor{otherrow}"
        lines.append("  " + shade + " " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    # Legend below the table
    lines += [
        r"\vspace{0.4em}",
        r"\caption*{\small Dual-winner: PRODIGY-LIG "
        r"$\Delta G_{\text{noelec}} \leq -7.0$\,kcal\,mol$^{-1}$ "
        r"\textit{and} catalytic Cys-SG $\leq 4$\,\AA. "
        r"Blue-shaded rows satisfy the near-attack criterion "
        r"($\leq 2.5$\,\AA); grey-shaded rows do not. "
        r"Substrate class colour coding matches Figs.~1--2.}",
        r"\label{tab:S1}",
        r"\end{table*}",
    ]

    out = REPORT / "table_S1_dual_winners.tex"
    out.write_text("\n".join(lines))
    print(f"✓ SI Table S1 saved → {out}")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print(f"Data: {len(df)} pairs, {len(df_ald)} aldehyde pairs")
    make_fig1()
    make_fig2()
    make_figS1()
    make_figS2()
    make_figS3()
    make_table_S1()
    print("\nAll outputs written to:", FIGS)
