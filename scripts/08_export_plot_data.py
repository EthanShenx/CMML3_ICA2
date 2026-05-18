"""
08_export_plot_data.py
Mirror the data-preparation logic of 07_make_figures.py and write one CSV
per panel into ShenYuchen/data2plot/. Performs no plotting; rendering is
delegated to the R scripts in ShenYuchen/visualization/.

Run from the project root:
    python scripts/08_export_plot_data.py
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, kruskal
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT       = SCRIPT_DIR.parent
SDF_DIR    = ROOT.parent / "Data" / "substrate"
RESULTS    = ROOT / "results"
OUT        = ROOT / "data2plot"
OUT.mkdir(exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
df     = pd.read_csv(RESULTS / "master_results.csv")
df_ald = df[df["is_aldehyde"] == True].copy()

PROTEINS = ["A8", "A9", "A10", "A11", "A12", "A13", "A14"]

# ── Substrate categories (copied verbatim from 07_make_figures.py) ────────────
CATEGORY_MAP = {
    "1_formaldehyde":              "Short-chain",
    "2_acetaldehyde":              "Short-chain",
    "3_propylaldehyde":            "Short-chain",
    "4_malonaldehyde":             "Short-chain",
    "5_butyraldehyde":             "Short-chain",
    "6_amylaldehyde":              "Short-chain",
    "9_2-hydroxypropionaldehyde":  "Short-chain",
    "7_glutaraldehyde":            "Medium-chain",
    "8_hexanal":                   "Medium-chain",
    "15_heptanal":                 "Medium-chain",
    "16_octanaldehyde":            "Medium-chain",
    "17_nonanal":                  "Medium-chain",
    "18_decanal":                  "Medium-chain",
    "19_dodecanal":                "Medium-chain",
    "24_Tetradecanal":             "Medium-chain",
    "25_hexadecanal":              "Medium-chain",
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
    "31_(2E,4E,6E,8E)-3,7-dimethyl-9-(2,6,6-trimethylcyclohexen-1-yl)nona-2,4,6,8-tetraenal": "Retinoid",
    "32_9-cis-Retinal":            "Retinoid",
    "33_11-cis-Retinol":           "Retinoid",
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
    if ligand in SUBSTRATE_ABBR:
        return SUBSTRATE_ABBR[ligand]
    parts = ligand.split("_", 1)
    name = parts[1] if len(parts) > 1 else parts[0]
    return name if len(name) <= maxlen else name[:maxlen - 1] + "…"


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    margin = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return max(0, centre - margin), min(1, centre + margin)


df["category"]     = df["ligand"].map(CATEGORY_MAP).fillna("Other")
df_ald["category"] = df_ald["ligand"].map(CATEGORY_MAP).fillna("Other")
df["short_name"]   = df["ligand"].apply(short_name)
df_ald["short_name"] = df_ald["ligand"].apply(short_name)

mask_d   = df_ald["CYS_SG_to_aldC_distance_A"].notna()
prod_dg  = df_ald.loc[mask_d & (df_ald["CYS_SG_to_aldC_distance_A"] <= 4.0), "DG_prediction_kcalmol"]
nprod_dg = df_ald.loc[mask_d & (df_ald["CYS_SG_to_aldC_distance_A"] >  4.0), "DG_prediction_kcalmol"]
n_prod   = (df_ald["CYS_SG_to_aldC_distance_A"] <= 4.0).sum()
n_tot    = mask_d.sum()


# =====================================================================
#  FIGURE 1
# =====================================================================
# 1b — chemical space (RDKit MW / logP per ligand)
def export_fig1b():
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        print("⚠ RDKit missing; skipping fig1b")
        return
    rows = []
    for lig in df["ligand"].unique():
        sdf = SDF_DIR / f"{lig}.sdf"
        if not sdf.exists():
            continue
        mol = Chem.MolFromMolFile(str(sdf), sanitize=True)
        if mol is None:
            continue
        rows.append({
            "ligand":   lig,
            "MW":       Descriptors.ExactMolWt(mol),
            "logP":     Descriptors.MolLogP(mol),
            "category": CATEGORY_MAP.get(lig, "Other"),
        })
    pd.DataFrame(rows).to_csv(OUT / "fig1b_chemspace.csv", index=False)


# 1c — DiffDock rank-1 confidence per pair
def export_fig1c():
    df[["protein", "rank1_confidence"]].dropna().to_csv(
        OUT / "fig1c_confidence.csv", index=False)


# 1d — affinity heatmap pivot, Ward-clustered column order
def export_fig1d():
    pivot   = df.pivot_table(index="protein", columns="ligand",
                              values="DG_prediction_kcalmol").reindex(PROTEINS)
    filled  = pivot.fillna(pivot.mean())
    col_link = linkage(pdist(filled.T.values, metric="euclidean"), method="ward")
    col_ord  = dendrogram(col_link, no_plot=True)["leaves"]
    pivot_s  = pivot.iloc[:, col_ord]

    long = pivot_s.reset_index().melt(id_vars="protein", var_name="ligand",
                                       value_name="DG")
    long["col_index"]   = long["ligand"].map({l: i for i, l in enumerate(pivot_s.columns)})
    long["lig_num"]     = long["ligand"].str.split("_", n=1).str[0]
    long["category"]    = long["ligand"].map(CATEGORY_MAP).fillna("Other")
    long.to_csv(OUT / "fig1d_heatmap_long.csv", index=False)

    # Dendrogram segment coordinates from the SAME linkage used to reorder
    # the columns. scipy renders each merge as a "U"; we emit one row per
    # straight segment so R can plot with geom_segment.
    dd = dendrogram(col_link, no_plot=True)
    # scipy's leaf x positions are 5, 15, 25, ... (i.e. 10*i + 5). Rescale to
    # the heatmap's col_index axis (0..N-1).
    leaves = dd["leaves"]                     # column order (length N)
    n_leaves = len(leaves)
    rescale = lambda x: (x - 5.0) / 10.0       # 5,15,25,... → 0,1,2,...
    seg_rows = []
    for icoord, dcoord in zip(dd["icoord"], dd["dcoord"]):
        # Each "U" has 4 points giving 3 segments: left-up, top, right-down.
        xs = [rescale(v) for v in icoord]
        ys = list(dcoord)
        seg_rows.append({"x": xs[0], "xend": xs[1], "y": ys[0], "yend": ys[1]})
        seg_rows.append({"x": xs[1], "xend": xs[2], "y": ys[1], "yend": ys[2]})
        seg_rows.append({"x": xs[2], "xend": xs[3], "y": ys[2], "yend": ys[3]})
    pd.DataFrame(seg_rows).to_csv(OUT / "fig1d_dendro_segments.csv", index=False)

    # Top-1 substrate per row (★)
    top1 = []
    for prot in PROTEINS:
        row = pivot_s.loc[prot].dropna()
        if not row.empty:
            best = row.idxmin()
            top1.append({"protein": prot, "ligand": best,
                         "col_index": list(pivot_s.columns).index(best)})
    pd.DataFrame(top1).to_csv(OUT / "fig1d_top1.csv", index=False)


# 1e — ΔG violin per isoform (aldehydes only, as in source script)
def export_fig1e():
    df_ald[["protein", "DG_prediction_kcalmol"]].to_csv(
        OUT / "fig1e_dg_violin.csv", index=False)


# 1f — top-5 substrates per protein
def export_fig1f():
    top5 = (df_ald.groupby("protein", group_keys=False)
                  .apply(lambda g: g.nsmallest(5, "DG_prediction_kcalmol")))
    top5["short_name"] = top5["ligand"].apply(lambda l: short_name(l, 20))
    top5[["protein", "ligand", "short_name", "DG_prediction_kcalmol"]] \
        .to_csv(OUT / "fig1f_top5.csv", index=False)


# 1g — pan-ALDH binding (legacy bar input + new substrate × isoform matrix)
def export_fig1g():
    thr = -6.0
    sel = (df_ald[df_ald.DG_prediction_kcalmol <= thr]
           .groupby("ligand")["protein"].nunique().reset_index())
    sel.columns = ["ligand", "n_prots"]
    all_ligs = pd.DataFrame({"ligand": df_ald["ligand"].unique()})
    sel = all_ligs.merge(sel, on="ligand", how="left").fillna(0)
    sel["lig_num"]  = sel["ligand"].str.split("_", n=1).str[0].astype(int)
    sel["category"] = sel["ligand"].map(CATEGORY_MAP).fillna("Other")
    sel = sel.sort_values("lig_num").reset_index(drop=True)
    sel.to_csv(OUT / "fig1g_pan_aldh.csv", index=False)

    # Substrate × isoform binding matrix for the dot-matrix viz: one row per
    # (protein, ligand) with is_strong = DG <= -6.
    mat = df_ald[["protein", "ligand", "DG_prediction_kcalmol"]].copy()
    mat["is_strong"] = (mat["DG_prediction_kcalmol"] <= thr).astype(int)
    mat["lig_num"]   = mat["ligand"].str.split("_", n=1).str[0].astype(int)
    mat["category"]  = mat["ligand"].map(CATEGORY_MAP).fillna("Other")
    mat["short_name"] = mat["ligand"].apply(short_name)
    # Attach the per-ligand promiscuity score (n isoforms with DG <= -6)
    mat = mat.merge(sel[["ligand", "n_prots"]], on="ligand", how="left")
    mat.rename(columns={"DG_prediction_kcalmol": "DG"}, inplace=True)
    mat.to_csv(OUT / "fig1g_matrix.csv", index=False)


# 1h — inter-isoform Pearson correlation (full matrix)
def export_fig1h():
    piv = df.pivot_table(index="protein", columns="ligand",
                          values="DG_prediction_kcalmol")
    piv = piv.fillna(piv.mean())
    corr = piv.T.corr(method="pearson").loc[PROTEINS, PROTEINS]
    corr.reset_index().to_csv(OUT / "fig1h_corr.csv", index=False)


# =====================================================================
#  FIGURE 2
# =====================================================================
def export_fig2a():
    df_ald.loc[mask_d, ["CYS_SG_to_aldC_distance_A"]] \
        .rename(columns={"CYS_SG_to_aldC_distance_A": "distance"}) \
        .to_csv(OUT / "fig2a_distance.csv", index=False)


def export_fig2b():
    df_ald[["protein", "CYS_SG_to_aldC_distance_A"]] \
        .dropna() \
        .rename(columns={"CYS_SG_to_aldC_distance_A": "distance"}) \
        .to_csv(OUT / "fig2b_dist_violin.csv", index=False)


def export_fig2c():
    sub = df_ald.loc[mask_d, [
        "protein", "ligand", "CYS_SG_to_aldC_distance_A",
        "DG_prediction_kcalmol", "rank1_confidence"]] \
        .rename(columns={"CYS_SG_to_aldC_distance_A": "distance",
                         "DG_prediction_kcalmol": "DG",
                         "rank1_confidence": "confidence"})
    sub.to_csv(OUT / "fig2c_scatter.csv", index=False)
    rho_c, p_c = spearmanr(sub["distance"], sub["DG"])
    pd.DataFrame([{"rho": rho_c, "p": p_c}]).to_csv(
        OUT / "fig2c_spearman.csv", index=False)


def export_fig2d():
    XTHRESH, YTHRESH = 4.0, -7.0

    def get_quad(row):
        if pd.isna(row["CYS_SG_to_aldC_distance_A"]):
            return "Neither"
        prod = row["CYS_SG_to_aldC_distance_A"] <= XTHRESH
        strg = row["DG_prediction_kcalmol"]       <= YTHRESH
        if prod and strg:    return "Productive high-affinity"
        if not prod and strg:return "High-affinity non-productive"
        if prod and not strg:return "Productive low-affinity"
        return "Neither"

    sub = df_ald.copy()
    sub["quad"] = sub.apply(get_quad, axis=1)
    out = sub[["protein", "ligand", "CYS_SG_to_aldC_distance_A",
               "DG_prediction_kcalmol", "quad"]] \
        .rename(columns={"CYS_SG_to_aldC_distance_A": "distance",
                         "DG_prediction_kcalmol": "DG"})
    out.to_csv(OUT / "fig2d_quadrant.csv", index=False)


def export_fig2e():
    rows = []
    for prot in PROTEINS:
        sub = df_ald[df_ald.protein == prot]["CYS_SG_to_aldC_distance_A"].dropna()
        k   = (sub <= 4.0).sum(); n = len(sub)
        frac = k / n if n else 0
        lo, hi = wilson_ci(k, n)
        rows.append({"protein": prot, "k": int(k), "n": int(n),
                     "frac": frac, "ci_lo": lo, "ci_hi": hi})
    pd.DataFrame(rows).to_csv(OUT / "fig2e_productive.csv", index=False)
    pd.DataFrame([{"overall_frac": n_prod / n_tot}]) \
        .to_csv(OUT / "fig2e_overall.csv", index=False)


def export_fig2f():
    rows = []
    for cat in CAT_ORDER:
        sub = df_ald[df_ald.category == cat]["CYS_SG_to_aldC_distance_A"].dropna()
        if len(sub) == 0:
            frac = 0.0
        else:
            frac = (sub <= 4.0).sum() / len(sub)
        rows.append({"category": cat, "frac": frac, "n": int(len(sub))})
    pd.DataFrame(rows).to_csv(OUT / "fig2f_class_prod.csv", index=False)
    pd.DataFrame([{"overall_frac": n_prod / n_tot}]) \
        .to_csv(OUT / "fig2f_overall.csv", index=False)


def export_fig2g():
    df_vg = df_ald[mask_d].copy()
    df_vg["pose"] = np.where(df_vg["CYS_SG_to_aldC_distance_A"] <= 4,
                              "Productive (≤ 4 Å)", "Non-productive (> 4 Å)")
    df_vg[["pose", "DG_prediction_kcalmol"]] \
        .rename(columns={"DG_prediction_kcalmol": "DG"}) \
        .to_csv(OUT / "fig2g_pose_violin.csv", index=False)
    stat_g, p_g = mannwhitneyu(prod_dg, nprod_dg, alternative="two-sided")
    pd.DataFrame([{"stat": stat_g, "p": p_g}]).to_csv(
        OUT / "fig2g_mwu.csv", index=False)


def export_fig2h():
    dw = (df_ald[(df_ald.DG_prediction_kcalmol <= -7.0) &
                 (df_ald.CYS_SG_to_aldC_distance_A <= 4.0)]
             .nsmallest(15, "DG_prediction_kcalmol"))
    dw = dw.sort_values("DG_prediction_kcalmol", ascending=False)
    dw["short_name"] = dw["ligand"].apply(lambda l: short_name(l, 30))
    dw["label"]      = dw["protein"] + "  " + dw["short_name"]
    dw[["protein", "ligand", "short_name", "label",
        "DG_prediction_kcalmol", "CYS_SG_to_aldC_distance_A"]] \
        .rename(columns={"DG_prediction_kcalmol": "DG",
                         "CYS_SG_to_aldC_distance_A": "distance"}) \
        .to_csv(OUT / "fig2h_dual_lollipop.csv", index=False)


def export_fig2i():
    feat = pd.read_csv(RESULTS / "interaction_features.csv")
    heavy_map = dict(zip(feat["ligand"], feat["HeavyAtoms"]))
    le_df = df_ald.copy()
    le_df["heavy"] = le_df["ligand"].map(heavy_map)
    le_df["LE"]    = le_df["DG_prediction_kcalmol"].abs() / le_df["heavy"]
    le_med = le_df["LE"].median()

    n_start = len(df_ald)
    n_aff   = int((df_ald["DG_prediction_kcalmol"] <= -7.0).sum())
    n_dual  = int(((df_ald["DG_prediction_kcalmol"] <= -7.0) &
                   (df_ald["CYS_SG_to_aldC_distance_A"] <= 4.0)).sum())
    n_triple = int(((le_df["DG_prediction_kcalmol"] <= -7.0) &
                    (le_df["CYS_SG_to_aldC_distance_A"] <= 4.0) &
                    (le_df["LE"] >= le_med)).sum())

    stages = [
        ("All aldehyde\npairs",                                       n_start),
        ("ΔG ≤ −7\nkcal mol⁻¹",                                       n_aff),
        ("+ Cys-SG\n≤ 4 Å",                                          n_dual),
        (f"+ LE ≥ {le_med:.2f}\nkcal mol⁻¹\natom⁻¹",                  n_triple),
    ]
    pd.DataFrame([{"stage": s, "n": n, "order": i}
                  for i, (s, n) in enumerate(stages)]) \
      .to_csv(OUT / "fig2i_funnel.csv", index=False)
    pd.DataFrame([{"le_med": le_med}]).to_csv(
        OUT / "fig2i_le_median.csv", index=False)


# =====================================================================
#  FIGURE S1
# =====================================================================
def export_figS1a():
    rc_counts = (df.groupby(["protein", "prodigy_rc"]).size()
                   .reset_index(name="n"))
    success = (rc_counts[rc_counts["prodigy_rc"] == 0]
                .set_index("protein")["n"].reindex(PROTEINS, fill_value=0))
    total   = df.groupby("protein").size().reindex(PROTEINS, fill_value=0)
    fail    = total - success
    pd.DataFrame({
        "protein":      PROTEINS,
        "success_pct":  (success.values / total.values * 100),
        "fail_pct":     (fail.values    / total.values * 100),
        "success_n":    success.values,
        "fail_n":       fail.values,
        "total":        total.values,
    }).to_csv(OUT / "figS1a_prodigy.csv", index=False)


def export_figS1b():
    ald_counts = df.groupby(["protein", "is_aldehyde"]).size().unstack(fill_value=0)
    ald_t = ald_counts.get(True,  pd.Series(0, index=PROTEINS))
    ald_f = ald_counts.get(False, pd.Series(0, index=PROTEINS))
    pd.DataFrame({
        "protein":       PROTEINS,
        "n_aldehyde":    ald_t.reindex(PROTEINS, fill_value=0).values,
        "n_nonaldehyde": ald_f.reindex(PROTEINS, fill_value=0).values,
    }).to_csv(OUT / "figS1b_smarts.csv", index=False)


def export_figS1c():
    sub = df.dropna(subset=["rank1_confidence"])[
        ["protein", "ligand", "rank1_confidence",
         "DG_prediction_kcalmol", "category"]]
    sub = sub.rename(columns={"DG_prediction_kcalmol": "DG",
                              "rank1_confidence": "confidence"})
    sub.to_csv(OUT / "figS1c_conf_dg.csv", index=False)
    rho_list = []
    for cat in CAT_ORDER:
        s = sub[sub["category"] == cat]
        if len(s) < 3:
            continue
        rho, p = spearmanr(s["confidence"], s["DG"])
        rho_list.append({"category": cat, "rho": rho, "p": p, "n": len(s)})
    pd.DataFrame(rho_list).to_csv(OUT / "figS1c_spearman.csv", index=False)


def export_figS1d():
    cats_d = CAT_ORDER[:-1]
    sub = df_ald[df_ald["category"].isin(cats_d)][
        ["category", "DG_prediction_kcalmol"]] \
        .rename(columns={"DG_prediction_kcalmol": "DG"})
    sub.to_csv(OUT / "figS1d_dg_class.csv", index=False)
    groups = [df_ald[df_ald.category == c].DG_prediction_kcalmol.dropna().values
              for c in cats_d if len(df_ald[df_ald.category == c]) >= 2]
    stat_kw, p_kw = kruskal(*groups)
    pd.DataFrame([{"stat": stat_kw, "p": p_kw}]).to_csv(
        OUT / "figS1d_kruskal.csv", index=False)


# =====================================================================
#  FIGURE S2  (one CSV; each per-isoform R script filters by protein)
# =====================================================================
def export_figS2():
    sub = df_ald[["protein", "ligand", "short_name", "category",
                  "DG_prediction_kcalmol", "CYS_SG_to_aldC_distance_A"]] \
        .rename(columns={"DG_prediction_kcalmol": "DG",
                         "CYS_SG_to_aldC_distance_A": "distance"})
    sub["is_productive"] = (sub["distance"] <= 4.0).fillna(False)
    sub.to_csv(OUT / "figS2_data.csv", index=False)


# =====================================================================
#  FIGURE S3
# =====================================================================
def export_figS3a():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import DataStructs
    except ImportError:
        print("⚠ RDKit missing; skipping figS3a")
        return
    fps_list, fp_names, fp_cats, fp_ligs = [], [], [], []
    for lig in df["ligand"].unique():
        sdf = SDF_DIR / f"{lig}.sdf"
        if not sdf.exists():
            continue
        mol = Chem.MolFromMolFile(str(sdf), sanitize=True)
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
        fps_list.append(fp)
        fp_names.append(lig.split("_")[0])
        fp_cats.append(CATEGORY_MAP.get(lig, "Other"))
        fp_ligs.append(lig)

    n = len(fps_list)
    dist_mat = np.zeros((n, n))
    for i in range(n):
        sims = DataStructs.BulkTanimotoSimilarity(fps_list[i], fps_list)
        dist_mat[i] = 1 - np.array(sims)
    pd.DataFrame(dist_mat, index=fp_names, columns=fp_names) \
        .to_csv(OUT / "figS3a_distmat.csv")
    pd.DataFrame({"label": fp_names, "ligand": fp_ligs,
                  "category": fp_cats}) \
        .to_csv(OUT / "figS3a_meta.csv", index=False)


def export_figS3b():
    cats_b = CAT_ORDER[:-1]
    sub = df_ald[df_ald["category"].isin(cats_b)][
        ["category", "DG_prediction_kcalmol"]] \
        .rename(columns={"DG_prediction_kcalmol": "DG"})
    sub.to_csv(OUT / "figS3b_dg_class.csv", index=False)
    groups_b = [df_ald[df_ald.category == c].DG_prediction_kcalmol.dropna().values
                for c in cats_b if len(df_ald[df_ald.category == c]) >= 2]
    stat_b, p_b = kruskal(*groups_b)
    pd.DataFrame([{"stat": stat_b, "p": p_b}]).to_csv(
        OUT / "figS3b_kruskal.csv", index=False)


def export_figS3c():
    med = (df_ald.groupby(["protein", "category"])
                ["CYS_SG_to_aldC_distance_A"].median().reset_index()
                .rename(columns={"CYS_SG_to_aldC_distance_A": "median_distance"}))
    med.to_csv(OUT / "figS3c_med_dist.csv", index=False)


def export_figS3d():
    cats_r = CAT_ORDER[:-1]
    rows = []
    for prot in PROTEINS:
        for cat in cats_r:
            sub = df_ald[(df_ald.protein == prot) & (df_ald.category == cat)]
            if len(sub) == 0:
                frac = 0.0
            else:
                frac = (sub.DG_prediction_kcalmol <= -6.0).sum() / len(sub)
            rows.append({"protein": prot, "category": cat, "frac": frac})
    pd.DataFrame(rows).to_csv(OUT / "figS3d_radar.csv", index=False)


# =====================================================================
#  TABLE S1
# =====================================================================
def export_tableS1():
    dw = (df_ald[(df_ald.DG_prediction_kcalmol <= -7.0) &
                 (df_ald.CYS_SG_to_aldC_distance_A <= 4.0)]
             .sort_values("DG_prediction_kcalmol")
             .head(20)
             .copy())
    dw["substrate"] = dw["ligand"].apply(lambda l: short_name(l, 35).replace("_", " "))
    dw["category"]  = dw["ligand"].map(CATEGORY_MAP).fillna("Other")
    dw["near"]      = dw["CYS_SG_to_aldC_distance_A"] <= 2.5
    dw[["protein", "ligand", "substrate", "DG_prediction_kcalmol",
        "CYS_SG_to_aldC_distance_A", "rank1_confidence", "category", "near"]] \
        .rename(columns={"DG_prediction_kcalmol": "DG",
                         "CYS_SG_to_aldC_distance_A": "distance",
                         "rank1_confidence": "confidence"}) \
        .to_csv(OUT / "tableS1_dual_winners.csv", index=False)


# =====================================================================
#  Main
# =====================================================================
if __name__ == "__main__":
    print(f"Source: {len(df)} pairs ({len(df_ald)} aldehydes)")
    for fn in [export_fig1b, export_fig1c, export_fig1d, export_fig1e,
               export_fig1f, export_fig1g, export_fig1h,
               export_fig2a, export_fig2b, export_fig2c, export_fig2d,
               export_fig2e, export_fig2f, export_fig2g, export_fig2h,
               export_fig2i,
               export_figS1a, export_figS1b, export_figS1c, export_figS1d,
               export_figS2,
               export_figS3a, export_figS3b, export_figS3c, export_figS3d,
               export_tableS1]:
        try:
            fn()
            print(f"  ✓ {fn.__name__}")
        except Exception as e:
            print(f"  ✗ {fn.__name__}: {e}")
    print(f"\nAll CSVs written to {OUT}")
