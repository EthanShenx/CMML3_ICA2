# visualization/

R rewrites of every panel produced by `scripts/07_make_figures.py`. The
script-level logic is unchanged — Python now only computes and exports
the numeric inputs, R does the rendering.

## Pipeline

```
scripts/08_export_plot_data.py
        │
        ▼
ShenYuchen/data2plot/*.csv     ← numeric inputs, one CSV per panel
        │
        ▼
ShenYuchen/visualization/NN_*.R ← one panel per script
        │
        ▼
ShenYuchen/figures/panels_R/NN_*.png
```

## Layout

| # | File | Panel |
|---|---|---|
| 01 | `01_fig1a_workflow.R` | Fig 1 (a) — pipeline schematic |
| 02 | `02_fig1b_chemspace.R` | Fig 1 (b) — MW vs logP |
| 03 | `03_fig1c_confidence_violin.R` | Fig 1 (c) — DiffDock confidence per isoform |
| 04 | `04_fig1d_affinity_heatmap.R` | Fig 1 (d) — ALDH × substrate ΔG heatmap |
| 05 | `05_fig1e_dg_violin.R` | Fig 1 (e) — ΔG violin per isoform |
| 06 | `06_fig1f_top5_lollipop.R` | Fig 1 (f) — top-5 substrates per isoform |
| 07 | `07_fig1g_pan_aldh_matrix.R` | Fig 1 (g) — substrate × isoform binding dot matrix |
| 08 | `08_fig1h_corr_heatmap.R` | Fig 1 (h) — inter-isoform Pearson r |
| 09 | `09_fig2a_distance_hist.R` | Fig 2 (a) — Cys-SG distance histogram |
| 10 | `10_fig2b_dist_violin.R` | Fig 2 (b) — distance violin per isoform |
| 11 | `11_fig2c_dg_vs_dist.R` | Fig 2 (c) — ΔG vs distance scatter |
| 12 | `12_fig2c1_pose_A8.R` | Fig 2 (c1) — A8 productive pose |
| 13 | `13_fig2c2_pose_A12.R` | Fig 2 (c2) — A12 non-productive pose |
| 14 | `14_fig2d_quadrant.R` | Fig 2 (d) — two-metric quadrant |
| 15 | `15_fig2e_productive_fraction.R` | Fig 2 (e) — productive fraction per isoform |
| 16 | `16_fig2f_class_productive.R` | Fig 2 (f) — productive fraction per class |
| 17 | `17_fig2g_pose_violin.R` | Fig 2 (g) — ΔG productive vs non-productive |
| 18 | `18_fig2h_dual_lollipop.R` | Fig 2 (h) — top-15 dual-winner pairs |
| 19 | `19_fig2i_filter_funnel.R` | Fig 2 (i) — three-gate filter funnel |
| 21 | `21_figS1a_prodigy_success.R` | Fig S1 (a) — PRODIGY-LIG success rates |
| 22 | `22_figS1b_smarts.R` | Fig S1 (b) — SMARTS aldehyde identification |
| 23 | `23_figS1c_conf_dg_class.R` | Fig S1 (c) — confidence-affinity per class |
| 24 | `24_figS1d_dg_boxplot.R` | Fig S1 (d) — ΔG distribution per class |
| 25–31 | `25_figS2a_A8.R` … `31_figS2g_A14.R` | Fig S2 (a)–(g) — per-isoform ranked substrates |
| 32 | `32_figS3a_dendrogram.R` | Fig S3 (a) — Tanimoto dendrogram |
| 33 | `33_figS3b_dg_class_violin.R` | Fig S3 (b) — ΔG violin per class |
| 34 | `34_figS3c_median_dist.R` | Fig S3 (c) — median distance per class |
| 35 | `35_figS3d_radar.R` | Fig S3 (d) — selectivity radar |
| 36 | `36_tableS1_dual_winners.R` | Table S1 — top-20 dual winners (PNG + LaTeX) |

`_palettes.R` defines the shared isoform colours (RColorBrewer "Paired"),
substrate-class colours, and path helpers. `_figS2_panel.R` is the shared
renderer used by 25–31.

## Reproducing

```bash
# 1. Compute panel inputs (writes data2plot/*.csv)
cd ShenYuchen
python scripts/08_export_plot_data.py

# 2. Render every panel (writes figures/panels_R/*.png + report/table_S1*.tex)
cd visualization
Rscript run_all.R
```

To render a single panel, just `Rscript NN_*.R` from this directory; the
helper resolves `data2plot/` and `figures/panels_R/` relative to itself.

## R dependencies

```r
install.packages(c("ggplot2", "dplyr", "tidyr",
                   "ggdendro", "gridExtra", "png", "scales"))
```
