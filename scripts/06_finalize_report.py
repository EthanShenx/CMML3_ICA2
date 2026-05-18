"""Fill the report.md placeholders with numbers from master_results.csv.

Idempotent: re-run after the pipeline finishes to update the figures.
"""
import math
import re
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

ROOT = Path("/root/ShenYuchen")
TPL = Path(__file__).parent.parent / "report" / "report.md"
OUT = Path(__file__).parent.parent / "report" / "report_filled.md"
df = pd.read_csv(ROOT / "master_results.csv")

n = len(df)
n_productive = int((df["CYS_SG_to_aldC_distance_A"] <= 4.0).sum())
mean_dg = df["DG_prediction_kcalmol"].mean()
median_dg = df["DG_prediction_kcalmol"].median()
top10 = df.nsmallest(10, "DG_prediction_kcalmol")[
    ["protein", "ligand", "DG_prediction_kcalmol", "CYS_SG_to_aldC_distance_A"]
]

prod = df[df["CYS_SG_to_aldC_distance_A"] <= 4.0]
rho, p = spearmanr(prod["DG_prediction_kcalmol"], prod["rank1_confidence"])

text = TPL.read_text()
text = text.replace("**N=… of 1560**", f"**N={n_productive} of {n}**")
text = re.sub(r"Spearman ρ = …", f"Spearman ρ = {rho:.2f} (p={p:.2g})", text)
text = text.replace("1560 ALDH–aldehyde pairs",
                    f"{n} ALDH–aldehyde pairs")

# Append a results table appendix
text += "\n\n## Supplementary Table 1 — top 10 strongest binders\n\n"
text += top10.to_markdown(index=False, floatfmt=".2f")

OUT.write_text(text)
print(f"wrote {OUT}")
print(f"n={n}, productive={n_productive}, median_DG={median_dg:.2f}")
