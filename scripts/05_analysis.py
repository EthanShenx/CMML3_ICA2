"""Merge prodigy + distances + DiffDock confidence into one master CSV
and produce the report figures.
"""
import csv
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

ROOT = Path("/root/ShenYuchen")
prodigy = pd.read_csv(ROOT / "prodigy_results.csv")
dists = pd.read_csv(ROOT / "distances.csv")
df = prodigy.merge(dists, on=["protein", "ligand"], how="left")
df.to_csv(ROOT / "master_results.csv", index=False)
print(f"merged rows: {len(df)}")

# Pretty ligand label = number only (sort key)
df["ligand_num"] = df["ligand"].str.split("_", n=1).str[0].astype(int)

# ---- Figure 1: heat-map of binding affinity (DG_noelec) ----
pivot = df.pivot_table(
    index="protein", columns="ligand_num", values="DG_prediction_kcalmol"
)
# sort index by ALDH letter then number
def _key(p):
    return (p[0], int(p[1:]))
pivot = pivot.reindex(sorted(pivot.index, key=_key))
pivot = pivot.reindex(sorted(pivot.columns), axis=1)

fig, ax = plt.subplots(figsize=(11, 5.5))
im = ax.imshow(pivot.values, cmap="RdYlBu", aspect="auto")
ax.set_xticks(range(len(pivot.columns)))
ax.set_xticklabels(pivot.columns, rotation=90, fontsize=6)
ax.set_yticks(range(len(pivot.index)))
ax.set_yticklabels(pivot.index, fontsize=8)
ax.set_xlabel("Substrate ID")
ax.set_ylabel("ALDH protein")
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cb.set_label("ΔG_noelec (kcal/mol)\n(more negative = stronger binding)")
fig.tight_layout()
fig.savefig(ROOT / "fig1_affinity_heatmap.png", dpi=300)
plt.close(fig)

# ---- Figure 2: catalytic-Cys ↔ aldehyde-C distance vs ΔG ----
fig, ax = plt.subplots(figsize=(6, 4.5))
sc = ax.scatter(
    df["CYS_SG_to_aldC_distance_A"],
    df["DG_prediction_kcalmol"],
    c=df["rank1_confidence"],
    cmap="viridis",
    s=18,
    alpha=0.75,
    edgecolor="none",
)
ax.set_xlabel("Cys SG ↔ aldehyde-C distance (Å)")
ax.set_ylabel("ΔG_noelec (kcal/mol)")
ax.axvline(4.0, ls="--", color="grey", lw=0.7,
           label="catalytic threshold (~4 Å)")
cb = fig.colorbar(sc, ax=ax)
cb.set_label("DiffDock rank-1 confidence")
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout()
fig.savefig(ROOT / "fig2_dg_vs_distance.png", dpi=300)
plt.close(fig)

# ---- summary statistics ----
print("\n=== summary ===")
print(df.select_dtypes(include="number").describe().round(2).T)

n_close = (df["CYS_SG_to_aldC_distance_A"] <= 4.0).sum()
print(f"\nPairs with catalytic-Cys ↔ aldehyde-C ≤ 4 Å: {n_close}/{len(df)}")
strong = df[df["DG_prediction_kcalmol"] < -7].sort_values("DG_prediction_kcalmol").head(10)
print("\nTop 10 strongest binders:")
print(strong[["protein", "ligand", "DG_prediction_kcalmol",
              "CYS_SG_to_aldC_distance_A", "rank1_confidence"]].to_string(index=False))
