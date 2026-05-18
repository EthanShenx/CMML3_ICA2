<p align="center">
  <img src="./assets/ZJE_logo.png" alt="ZJU-UoE Institute logo" width="780">
</p>

---

<h1 align="center">CMML3 ICA2: Systematic Blind Docking of Human ALDH Isoforms × Aldehyde Substrates</h1>

<p align="center">
  <img src="https://img.shields.io/badge/Python-3.x-3776AB?logo=python&logoColor=white" alt="Python badge">
  <img src="https://img.shields.io/badge/R-4.x-276DC3?logo=r&logoColor=white" alt="R badge">
  <img src="https://img.shields.io/badge/Docking-DiffDock--L%201.1-1f883d" alt="DiffDock badge">
  <img src="https://img.shields.io/badge/Scoring-PRODIGY--LIG-6f42c1" alt="PRODIGY badge">
  <img src="https://img.shields.io/badge/Outputs-Data%20%7C%20Figures-orange" alt="Outputs badge">
</p>

---

## Overview

This folder is a compact, reviewer-facing reproduction package for ICA2. It contains:

- the full analysis pipeline (DiffDock-L blind docking → PRODIGY-LIG ΔG scoring → catalytic-Cys geometry → master CSV) under [`scripts/`](./scripts/)
- the per-panel R rendering pipeline under [`visualization/`](./visualization/), driven by numeric inputs in [`data2plot/`](./data2plot/)
- the precomputed master tables in [`results/`](./results/) and rendered panels / figures in [`figures/`](./figures/)

No absolute paths are required at the visualization stage. The analysis stage is configured for a GPU cluster (`/root/ShenYuchen/`) and is documented for reference; all downstream steps run locally from the precomputed `results/master_results.csv`.

## Repository Layout

```text
ShenYuchen_github_repo/
  README.md                           # this file
  README_server.md                    # cluster-side install / run notes (Section 2 of the brief)
  notebook.ipynb                      # reviewer-facing exploratory notebook
  manual/
    Website_Tasks_Manual.md           # browser walkthrough (Section 1 of the brief)
  webpage_test/
    diffdock_ui_practice/             # screenshots of the DiffDock UI exercise

  scripts/                            # analysis pipeline (Python + shell)
    01_make_diffdock_csv.py           # build 7 × 60 batch CSV
    02_combine.py                     # protein + rank-1 pose -> combined PDB (pure Python)
    02_combine_pymol.py               # PyMOL alternative (kept for reference)
    03_run_prodigy.py                 # batch PRODIGY-LIG ΔG, parses to CSV
    04_distances.py                   # catalytic Cys-SG ↔ aldehyde-C geometry
    05_analysis.py                    # merge PRODIGY + distance + confidence → master CSV
    06_finalize_report.py             # injects numbers into the report body
    07_make_figures.py                # all main + SI figures (matplotlib)
    08_export_plot_data.py            # writes per-panel CSVs into data2plot/ for R
    09_interaction_analysis.py        # contact-type fingerprints (Fig S4)
    parallel_dl.py                    # 8-thread HTTP-Range downloader (ESM2 weights)
    launch.sh                         # waits for ESM2, then runs master_pipeline.sh
    master_pipeline.sh                # 2-GPU DiffDock + downstream chain
    run_all.sh                        # serial-mode reference script

  visualization/                      # one R script per panel (ggplot2)
    _palettes.R                       # shared colours + helpers
    _figS2_panel.R                    # shared renderer for the per-isoform set
    01_fig1a_workflow.R … 36_tableS1_dual_winners.R
    run_all.R                         # render every panel
    README.md                         # panel-by-panel index

  data2plot/                          # numeric inputs for R (one CSV per panel)
  results/                            # merged pipeline tables (master CSV, etc.)
  figures/
    fig1_landscape.pdf                # final composed main figure 1
    fig2_geometry.pdf                 # final composed main figure 2
    figS1_qc.pdf … figS4_multicriterion.pdf   # final composed SI figures
    panels_R/*.png                    # R-rendered panel PNGs used to compose finals
    pose_panel/                       # PyMOL pose-render assets (Fig 2c1 / 2c2)

  assets/ZJE_logo.png                 # logo
```

## Core Code

| File | Purpose |
| --- | --- |
| [`scripts/01_make_diffdock_csv.py`](./scripts/01_make_diffdock_csv.py) | Build the DiffDock batch CSV for every protein × ligand pair |
| [`scripts/02_combine.py`](./scripts/02_combine.py) | Merge rank-1 ligand pose with its protein into a PRODIGY-ready PDB |
| [`scripts/03_run_prodigy.py`](./scripts/03_run_prodigy.py) | Batch `prodigy_lig` over every combined PDB; parse ΔG into a CSV |
| [`scripts/04_distances.py`](./scripts/04_distances.py) | Compute catalytic Cys-SG ↔ aldehyde-C distance per pose |
| [`scripts/05_analysis.py`](./scripts/05_analysis.py) | Merge PRODIGY + distance + DiffDock confidence into `master_results.csv` |
| [`scripts/07_make_figures.py`](./scripts/07_make_figures.py) | Compose every final main + SI figure (PDF + PNG) |
| [`scripts/08_export_plot_data.py`](./scripts/08_export_plot_data.py) | Export per-panel numeric CSVs into `data2plot/` |
| [`scripts/09_interaction_analysis.py`](./scripts/09_interaction_analysis.py) | Contact-type fingerprints and decision matrix (Fig S4) |
| [`visualization/_palettes.R`](./visualization/_palettes.R) | Shared isoform / substrate-class colours and rendering helpers |
| [`visualization/run_all.R`](./visualization/run_all.R) | Render every numbered panel script into `figures/panels_R/` |

## Output Locations

| Location | Contents |
| --- | --- |
| [`results/`](./results/) | Merged pipeline tables (`master_results.csv`, `distances.csv`, `prodigy_results.csv`, `prodigy_contacts.csv`, `interaction_features.csv`) |
| [`data2plot/`](./data2plot/) | One CSV per panel — the numeric inputs consumed by every R script |
| [`figures/`](./figures/) | Final composed main + SI figure PDFs and the rendered per-panel PNGs in [`figures/panels_R/`](./figures/panels_R/) |

## What Changed Between Stages

| Stage | Tool | Input | Output |
| --- | --- | --- | --- |
| 1. Pair list | `01_make_diffdock_csv.py` | `Data/proteinPDB/`, `Data/substrate/` | DiffDock batch CSV |
| 2. Blind docking | DiffDock-L 1.1 (external) | batch CSV | rank-1 ligand SDF per pair |
| 3. Combine | `02_combine.py` | protein PDB + rank-1 SDF | PRODIGY-ready combined PDB |
| 4. Affinity | `03_run_prodigy.py` | combined PDBs | `results/prodigy_results.csv` (ΔG, contacts) |
| 5. Geometry | `04_distances.py` | combined PDBs | `results/distances.csv` (Cys-SG ↔ C=O distance) |
| 6. Merge | `05_analysis.py` | prodigy + distance + DiffDock confidence | `results/master_results.csv` |
| 7. Figures | `07_make_figures.py` | `master_results.csv` | `figures/fig*.pdf` + SI |
| 8. Plot data | `08_export_plot_data.py` | `master_results.csv` | `data2plot/*.csv` |
| 9. Render | `visualization/*.R` | `data2plot/*.csv` | `figures/panels_R/*.png` |
| 10. Interactions | `09_interaction_analysis.py` | `prodigy_contacts.csv` | `interaction_features.csv` + Fig S4 panels |

## File Relationships (Call Graph)

### Analysis pipeline (cluster-side, GPU required for step 2)

```text
01_make_diffdock_csv.py
  -> diffdock_input.csv

[external] DiffDock-L 1.1
  -> rank1 ligand SDF per protein × substrate pair

02_combine.py
  -> combined PDBs (chain A + chain L)

03_run_prodigy.py
  -> results/prodigy_results.csv

04_distances.py
  -> results/distances.csv

05_analysis.py
  -> merges prodigy + distances + diffdock confidence
  -> results/master_results.csv
```

### Visualization pipeline (local)

```text
07_make_figures.py
  -> figures/fig1_landscape.{pdf,png}
  -> figures/fig2_geometry.{pdf,png}
  -> figures/figS1_qc.{pdf,png}
  -> figures/figS2_isoforms.{pdf,png}
  -> figures/figS3_chemspace.{pdf,png}
  -> report/table_S1_dual_winners.tex

08_export_plot_data.py
  -> data2plot/fig1b_chemspace.csv … data2plot/tableS1_dual_winners.csv

visualization/run_all.R
  -> sources every NN_*.R panel script
  -> reads data2plot/*.csv
  -> writes figures/panels_R/*.png
```

### Interaction analysis (Fig S4)

```text
09_interaction_analysis.py
  -> results/interaction_features.csv
  -> figures/figS4_multicriterion.{pdf,png}
```

## Quick Start

### 1. Visualization only (precomputed `results/master_results.csv`)

This is the fastest path: it reuses the docking/scoring results already
checked into [`results/`](./results/) and re-renders every panel locally.

```bash
cd ShenYuchen_github_repo

# (a) Compose finals via matplotlib + export per-panel CSVs for R
python scripts/07_make_figures.py
python scripts/08_export_plot_data.py
python scripts/09_interaction_analysis.py

# (b) Render every R panel (writes figures/panels_R/*.png)
cd visualization
Rscript run_all.R
```

### 2. Full analysis pipeline (GPU cluster)

The DiffDock + PRODIGY stages are configured for a GPU cluster with
DiffDock-L 1.1 weights at `/root/DiffDock/workdir/v1.1/` and the ESM2-650M
cache at `/root/.cache/torch/hub/checkpoints/`. From the cluster:

```bash
# one-shot — runs DiffDock on both GPUs in parallel, then everything downstream
nohup bash scripts/launch.sh </dev/null >/dev/null 2>&1 &
disown

# inspect progress
tail -F launch.log master.log run_a.log run_b.log
```

The default scope is **7 proteins (ALDH1A1–1A3, 1A8–1A14, equivalents) ×
60 aldehydes** = 420 protein–ligand pairs.

## Key Parameters (As Used For The Precomputed Outputs)

### DiffDock-L 1.1

- top-k pose retention: rank-1 only (used for downstream geometry / scoring)
- batch CSV is produced once by `01_make_diffdock_csv.py` and split into two
  GPU shards (`csv_a.csv`, `csv_b.csv`) by `master_pipeline.sh`

### PRODIGY-LIG (`03_run_prodigy.py`)

- protein chain: `A`, ligand chain: `L`, ligand residue: `UNK`
- reported metric: `ΔG_noelec` (electrostatic-free, used as the headline ΔG)
- dual-winner threshold: `ΔG ≤ -7.0 kcal mol⁻¹` **and** Cys-SG ≤ 4 Å

### Geometry (`04_distances.py`)

- catalytic Cys per isoform is hard-coded from the canonical FNxxGQ[TS]C
  motif; A10 uses Cys207 because the spreadsheet Cys244 is absent from the
  deposited PDB
- aldehyde carbon detected by RDKit SMARTS on the rank-1 ligand SDF
- Bürgi–Dunitz threshold for the *near-attack* sub-class: **≤ 2.5 Å**;
  general productive threshold: **≤ 4.0 Å**

### Figure composition (`07_make_figures.py`)

- font: Arial (`matplotlib.rcParams["font.family"] = "Arial"`), `pdf.fonttype = 42`
- isoform palette: `RColorBrewer "Paired"` (A8 – A14)
- substrate-class palette: 6 classes (Aromatic, Medium-chain, Retinoid,
  Short-chain, α,β-Unsaturated, Other)
- dendrogram: Ward linkage on Tanimoto distances of Morgan fingerprints (r = 2)

## ICA Figures

- Main figures (final, as submitted): [`figures/fig1_landscape.pdf`](./figures/fig1_landscape.pdf), [`figures/fig2_geometry.pdf`](./figures/fig2_geometry.pdf)
- Supplementary figures: [`figures/figS1_qc.pdf`](./figures/figS1_qc.pdf), [`figures/figS2_isoforms.pdf`](./figures/figS2_isoforms.pdf), [`figures/figS3_chemspace.pdf`](./figures/figS3_chemspace.pdf), [`figures/figS4_multicriterion.pdf`](./figures/figS4_multicriterion.pdf)
- Panel PNGs used to build the finals: [`figures/panels_R/*.png`](./figures/panels_R/)
- Pose-render assets: [`figures/pose_panel/`](./figures/pose_panel/) (A8 productive + A12 non-productive)

## Dependencies

### Python (analysis + figure composition)

```bash
pip install numpy pandas scipy matplotlib seaborn rdkit biopython prodigy-lig
```

### R (per-panel visualization)

```r
install.packages(c("ggplot2", "dplyr", "tidyr",
                   "ggdendro", "patchwork", "gridExtra",
                   "png", "scales", "RColorBrewer", "ragg"))
```

Fonts: Arial must be available system-side; the R helper uses
`ragg::agg_png` to render Arial reliably on macOS / Linux.
