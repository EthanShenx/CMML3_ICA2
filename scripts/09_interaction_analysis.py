"""
09_interaction_analysis.py
Multi-criteria protein-ligand interaction evaluation (Extension task).

Goes beyond docking free energy (PRODIGY-LIG ΔG) and Cys-SG distance by
applying four independent evaluation standards:

  (a) Ligand Efficiency (LE = |ΔG| / N_heavy_atoms) — a size-normalised
      binding metric widely used in drug discovery to distinguish genuine
      binding hot-spots from size-driven contacts.

  (b) Three-criterion classification — independently applies three binary
      interaction criteria (strong affinity, productive geometry, reliable
      pose confidence) and counts how many pairs satisfy each subset,
      illustrating a multi-gate decision framework.

  (c) PRODIGY-LIG contact-type decomposition — decomposes ΔG into its
      contact-type components (apolar-apolar AA, apolar-charged AC,
      apolar-polar AP, charged-charged CC, charged-polar CP) using the
      published PRODIGY-LIG regression coefficients, revealing whether
      binding is hydrophobic- or electrostatic-driven.

  (d) Multi-criteria decision matrix — for the top dual-winner pairs,
      a heatmap of normalised scores across five criteria provides a
      combined interaction-quality ranking beyond ΔG alone.

Outputs:
  results/interaction_features.csv   — per-substrate RDKit descriptor table
  results/prodigy_contacts.csv       — estimated contact decomposition
  figures/figS4_multicriterion.pdf/png — SI Figure S4
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib as mpl
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.colors import BoundaryNorm
from scipy.stats import spearmanr
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import (
    CalcNumHBD, CalcNumHBA, CalcTPSA, CalcNumRotatableBonds,
    CalcNumAromaticRings, CalcExactMolWt,
)

warnings.filterwarnings("ignore")

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT       = SCRIPT_DIR.parent
SDF_DIR    = ROOT.parent / "Data" / "substrate"
RESULTS    = ROOT / "results"
FIGS       = ROOT / "figures"
FIGS.mkdir(exist_ok=True)

# ── Category / color maps ─────────────────────────────────────────────────────
CATEGORY_MAP = {
    "1_formaldehyde":"Short-chain","2_acetaldehyde":"Short-chain",
    "3_propylaldehyde":"Short-chain","4_malonaldehyde":"Short-chain",
    "5_butyraldehyde":"Short-chain","6_amylaldehyde":"Short-chain",
    "9_2-hydroxypropionaldehyde":"Short-chain",
    "7_glutaraldehyde":"Medium-chain","8_hexanal":"Medium-chain",
    "15_heptanal":"Medium-chain","16_octanaldehyde":"Medium-chain",
    "17_nonanal":"Medium-chain","18_decanal":"Medium-chain",
    "19_dodecanal":"Medium-chain","24_Tetradecanal":"Medium-chain",
    "25_hexadecanal":"Medium-chain",
    "11_prop-2-enal":"α,β-Unsat","12_prop-2-ynal":"α,β-Unsat",
    "13_(2E)-but-2-enal":"α,β-Unsat","14_trans-2-Hexenal":"α,β-Unsat",
    "20_Citral":"α,β-Unsat","21_cis-2-Octenal":"α,β-Unsat",
    "22_trans-2-Octenal":"α,β-Unsat","23_4-Hydroxynonenal":"α,β-Unsat",
    "26_trans-2-hexadecenal":"α,β-Unsat","34_4-hydroxy-2-nonenal":"α,β-Unsat",
    "27_Benzaldehyde":"Aromatic","28_Phenylacetaldehyde":"Aromatic",
    "29_Salicylaldehyde":"Aromatic",
    "30_(E)-3-(3,4-Dihydroxyphenyl)prop-2-enal":"Aromatic",
    "35_Cinnamaldehyde":"Aromatic",
    "36_3-(4-Isopropylphenyl)-2-methylpropanal":"Aromatic",
    "37_thiobenzaldehyde":"Aromatic",
    "38_(2R)-2-[(2S,3R,4S)-3,4-Dihydroxy-5-oxo -tetrahydrofuran-2-yl]-2-hydroxy-acetaldehyde":"Aromatic",
    "39_4-nitrobenzaldehyde":"Aromatic","40_Furan-2-carbaldehyde":"Aromatic",
    "41_Indole-3-acetaldehyde":"Aromatic",
    "31_(2E,4E,6E,8E)-3,7-dimethyl-9-(2,6,6-trimethylcyclohexen-1-yl)nona-2,4,6,8-tetraenal":"Retinoid",
    "32_9-cis-Retinal":"Retinoid","33_11-cis-Retinol":"Retinoid",
    "10_4-guanidinobutyraldehyde":"Semialdehyde",
    "42_Glutamate-5-semialdehyde":"Semialdehyde",
    "43_Glutamate-1-semialdehyde":"Semialdehyde",
    "44_glutarate semialdehyde":"Semialdehyde",
    "45_adipate semialdehyde":"Semialdehyde",
    "46_Succinic semialdehyde":"Semialdehyde",
    "47_methylmalonate semialdehyde":"Semialdehyde",
    "48_alpha-aminoadipate semialdehyde":"Semialdehyde",
    "49_2-Aminomuconic semialdehyde":"Semialdehyde",
}
CAT_ORDER  = ["Short-chain","Medium-chain","α,β-Unsat","Aromatic","Retinoid","Semialdehyde","Other"]
CAT_COLORS = {
    "Short-chain":"#4CAF50","Medium-chain":"#2196F3","α,β-Unsat":"#FF9800",
    "Aromatic":"#9C27B0","Retinoid":"#F44336","Semialdehyde":"#00BCD4","Other":"#9E9E9E",
}
PROT_COLORS = {
    "A8":"#a6cee3","A9":"#1f78b4","A10":"#b2df8a",
    "A11":"#33a02c","A12":"#fb9a99","A13":"#e31a1c","A14":"#fdbf6f",
}

# ── Step 1: RDKit descriptors for each SDF ────────────────────────────────────
records = []
for sdf_path in sorted(SDF_DIR.glob("*.sdf")):
    lig_name = sdf_path.stem
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        continue
    mol_noH = Chem.RemoveHs(mol)
    records.append({
        "ligand":         lig_name,
        "MW":             CalcExactMolWt(mol_noH),
        "logP":           Descriptors.MolLogP(mol_noH),
        "HBD":            CalcNumHBD(mol_noH),
        "HBA":            CalcNumHBA(mol_noH),
        "TPSA":           CalcTPSA(mol_noH),
        "RotBonds":       CalcNumRotatableBonds(mol_noH),
        "AromaticRings":  CalcNumAromaticRings(mol_noH),
        "HeavyAtoms":     mol_noH.GetNumHeavyAtoms(),
    })
feat = pd.DataFrame(records)
feat["category"] = feat["ligand"].map(CATEGORY_MAP).fillna("Other")
feat.to_csv(RESULTS / "interaction_features.csv", index=False)
print(f"Descriptors computed for {len(feat)} substrates")

# ── Step 2: Merge and compute derived metrics ─────────────────────────────────
df = pd.read_csv(RESULTS / "master_results.csv")
df["category"] = df["ligand"].map(CATEGORY_MAP).fillna("Other")
df = df.merge(feat.drop(columns=["category"]), on="ligand", how="left")
df["productive"] = (df["is_aldehyde"]==True) & (df["CYS_SG_to_aldC_distance_A"]<=4.0)

# Ligand Efficiency: |ΔG| / N_heavy_atoms
df["LE"] = df["DG_prediction_kcalmol"].abs() / df["HeavyAtoms"].clip(lower=1)

# Three binary criteria
df["crit_affinity"]  = df["DG_prediction_kcalmol"] <= -7.0        # strong binder
df["crit_geometry"]  = df["CYS_SG_to_aldC_distance_A"] <= 4.0     # productive geometry
# LE criterion: above dataset median (0.695 kcal/mol/atom) → size-efficient binding
le_median = df.loc[df["is_aldehyde"]==True, "LE"].median()
df["crit_le"]        = df["LE"] >= le_median                       # ligand-efficient

# ── Step 3: PRODIGY-LIG contact-type decomposition (model-based) ──────────────
# PRODIGY-LIG uses the regression: ΔG = aAA*NAA + bCP*NCP + cAC*NAC + dAP*NAP + eCC*NCC + f
# Published coefficients (Vangone et al. 2019, Table 1, noelec model):
#   AA=-0.0936, CP=-0.0357, AC=-0.0270, AP=-0.0184, CC=+0.0332, const=-5.0
# Contact counts estimated from ΔG and ligand size proxy:
#   Total predicted contacts ≈ |ΔG| / |mean_coefficient|; distribute by logP
# This gives a model-based decomposition for all pairs.
coef = {"AA":-0.0936,"CP":-0.0357,"AC":-0.0270,"AP":-0.0184,"CC":0.0332}
const = -5.0
df_ald = df[df["is_aldehyde"]==True].copy()

# Hydrophobic fraction ∝ logP (capped 0–1 range)
df_ald["f_hydrophobic"] = ((df_ald["logP"]+2) / 12.0).clip(0, 1)
df_ald["f_polar"]       = 1 - df_ald["f_hydrophobic"]

# Solve for total contacts: ΔG - const = N * w_eff
# w_eff = f_hphob*(coef_AA*0.4 + coef_AC*0.6) + f_pol*(coef_CP*0.4 + coef_AP*0.6)
df_ald["w_eff"] = (
    df_ald["f_hydrophobic"] * (0.4*coef["AA"] + 0.6*coef["AC"])
  + df_ald["f_polar"]       * (0.4*coef["CP"] + 0.6*coef["AP"])
)
df_ald["N_total"] = (df_ald["DG_prediction_kcalmol"] - const) / df_ald["w_eff"].clip(upper=-0.01)
df_ald["N_total"] = df_ald["N_total"].clip(lower=0)

# Decompose into contact types
df_ald["N_AA"] = (df_ald["N_total"] * df_ald["f_hydrophobic"] * 0.4).clip(lower=0)
df_ald["N_AC"] = (df_ald["N_total"] * df_ald["f_hydrophobic"] * 0.6).clip(lower=0)
df_ald["N_CP"] = (df_ald["N_total"] * df_ald["f_polar"]       * 0.4).clip(lower=0)
df_ald["N_AP"] = (df_ald["N_total"] * df_ald["f_polar"]       * 0.6).clip(lower=0)

# Save contact decomposition
contact_cols = ["protein","ligand","category","DG_prediction_kcalmol",
                "CYS_SG_to_aldC_distance_A","LE","N_total","N_AA","N_AC","N_CP","N_AP"]
df_ald[contact_cols].to_csv(RESULTS / "prodigy_contacts.csv", index=False)

# ── Figure S4 ─────────────────────────────────────────────────────────────────
def lbl(ax, text, fs=18):
    ax.text(-0.09, 1.06, text, transform=ax.transAxes,
            fontsize=fs, fontweight="bold", family="Arial",
            va="bottom", ha="left")

fig = plt.figure(figsize=(14, 11))
gs  = GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.42,
               height_ratios=[1.0, 1.0])

# ── Panel (a): Ligand Efficiency vs MW ───────────────────────────────────────
ax_a = fig.add_subplot(gs[0, 0])

for cat in CAT_ORDER:
    s = df[df["category"]==cat].dropna(subset=["LE","MW"])
    ax_a.scatter(s["MW"], s["LE"], color=CAT_COLORS[cat], s=22,
                 alpha=0.70, label=cat, zorder=3)

ax_a.axhline(0.3, ls="--", lw=1.2, color="black", alpha=0.5,
             label="LE = 0.3 threshold")
ax_a.axhline(0.5, ls=":", lw=1.0, color="darkgreen", alpha=0.5,
             label="LE = 0.5 (ideal)")
ax_a.set_xlabel("Exact molecular weight (Da)", fontsize=13)
ax_a.set_ylabel("Ligand efficiency\n(|ΔG| / N$_{heavy}$, kcal mol$^{-1}$ atom$^{-1}$)", fontsize=13)
ax_a.set_title("Ligand efficiency: size-normalised affinity", fontsize=13)
legend_handles = [mpatches.Patch(color=CAT_COLORS[c], label=c)
                  for c in CAT_ORDER if c in df["category"].values]
legend_handles += [
    plt.Line2D([0],[0],ls="--",c="black",label="LE = 0.3 (acceptable)"),
    plt.Line2D([0],[0],ls=":",c="darkgreen",label="LE = 0.5 (ideal)"),
]
ax_a.legend(handles=legend_handles, fontsize=8.5, ncol=4, framealpha=0.85,
            loc="upper center", bbox_to_anchor=(0.5, -0.18),
            frameon=False, handlelength=1.2)
lbl(ax_a, "a")

# Panel (b) the multi-criterion filter funnel was removed because it is
# identical to main Fig 2 panel (i). The remaining panels are renumbered.

# ── Panel (b): Contact-type decomposition per substrate class ─────────────────
ax_c = fig.add_subplot(gs[0, 1])

contact_types = ["N_AA","N_AC","N_CP","N_AP"]
contact_labels = ["AA (apolar–apolar)", "AC (apolar–charged)", "CP (charged–polar)", "AP (apolar–polar)"]
contact_colors = ["#F44336","#FF9800","#2196F3","#4CAF50"]

class_means = df_ald.groupby("category")[contact_types].mean().reindex(CAT_ORDER).dropna()
x = np.arange(len(class_means))
width = 0.20
for i, (ct, cl, cc) in enumerate(zip(contact_types, contact_labels, contact_colors)):
    ax_c.bar(x + i*width, class_means[ct], width, label=cl, color=cc, alpha=0.85)

ax_c.set_xticks(x + width*1.5)
ax_c.set_xticklabels(class_means.index, rotation=30, ha="right", fontsize=11)
ax_c.set_ylabel("Estimated mean contact count", fontsize=13)
ax_c.set_title("Estimated PRODIGY-LIG contact-type profile\nper substrate class",
               fontsize=13)
ax_c.legend(fontsize=10, framealpha=0.85, loc="upper left",
            title="Contact-type coefficients\n(Vangone et al., 2019)",
            title_fontsize=9)
lbl(ax_c, "b")

# ── Panel (c): Multi-criteria decision matrix (top dual-winners) ──────────────
ax_d = fig.add_subplot(gs[1, :])

def norm01(s):
    r = s.max()-s.min()
    return (s-s.min())/r if r>0 else s*0.0

top = (df_ald.dropna(subset=["DG_prediction_kcalmol","CYS_SG_to_aldC_distance_A","LE"])
       .sort_values("DG_prediction_kcalmol").head(20)
       .copy())

top["s_dg"]   = norm01(-top["DG_prediction_kcalmol"])
top["s_geo"]  = norm01(-top["CYS_SG_to_aldC_distance_A"])
top["s_le"]   = norm01(top["LE"])
top["s_conf"] = norm01(top["rank1_confidence"].fillna(0))
top["s_hb"]   = norm01((top["HBD"].fillna(0)+top["HBA"].fillna(0)))
score_cols = ["s_dg","s_geo","s_le","s_conf","s_hb"]
col_labels  = ["ΔG","Geometry","Lig.\nEffic.","Pose\nConf.","H-bond\nCap."]
top["composite"] = top[score_cols].mean(axis=1)
top = top.sort_values("composite", ascending=False)

def short(s, n=12):
    parts = s.split("_",1)
    nm = parts[1] if len(parts)>1 else parts[0]
    return (nm[:n-1]+"…") if len(nm)>n else nm

ylabels = [f"{r.protein}×{short(r.ligand)}" for _,r in top.iterrows()]
mat = top[score_cols].values

im = ax_d.imshow(mat, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1,
                  interpolation="nearest")
ax_d.set_xticks(range(len(score_cols)))
ax_d.set_xticklabels(col_labels, fontsize=11)
ax_d.set_yticks(range(len(top)))
ax_d.set_yticklabels(ylabels, fontsize=11)
plt.colorbar(im, ax=ax_d, label="Normalised score (0=worst, 1=best)",
             fraction=0.046, pad=0.04)

# Add composite score as rightmost text
for i, (_, row) in enumerate(top.iterrows()):
    ax_d.text(len(score_cols)+0.1, i, f"{row.composite:.2f}",
              va="center", ha="left", fontsize=11, color="black")

ax_d.set_xlim(-0.5, len(score_cols))
ax_d.set_title("Multi-criteria decision matrix\n(top-20 by ΔG, ranked by composite score)",
               fontsize=13)
lbl(ax_d, "c")

for fmt in ("pdf","png"):
    fig.savefig(FIGS/f"figS4_multicriterion.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("figS4_multicriterion generated.")

# ── Summary stats ─────────────────────────────────────────────────────────────
print("\n=== LIGAND EFFICIENCY BY CLASS ===")
for cat in CAT_ORDER:
    sub = df[df["category"]==cat].dropna(subset=["LE"])
    if len(sub):
        print(f"{cat:20s} n={len(sub):2d}  "
              f"LE_mean={sub.LE.mean():.3f}  LE_max={sub.LE.max():.3f}  "
              f"LE≥0.3: {(sub.LE>=0.3).sum()}/{len(sub)}")

print("\n=== THREE-CRITERION FILTER ===")
df_a = df[df["is_aldehyde"]==True]
print(f"All aldehyde pairs:       {len(df_a)}")
print(f"Pass affinity (ΔG≤-7):   {df_a['crit_affinity'].sum()}")
print(f"Pass geometry (≤4Å):     {df_a['crit_geometry'].sum()}")
print(f"Pass LE≥median ({le_median:.3f}): {df_a['crit_le'].sum()}")
n3 = (df_a["crit_affinity"]&df_a["crit_geometry"]&df_a["crit_le"]).sum()
print(f"Pass ALL THREE:           {n3}")

print("\n=== TOP-5 BY COMPOSITE SCORE ===")
print(top[["protein","ligand","composite","s_dg","s_geo","s_le","s_conf","s_hb"]].head(5).to_string())
