# Server-side pipeline (Section 2 of the ICA)

The whole pipeline lives in `/mnt/ShenYuchen/`. It is split into five
self-contained Python scripts; each is idempotent (skip files that
already exist) so you can re-run any step.

## 0. Prerequisites (one-off setup)

The conda environment `diffdock` provides DiffDock + RDKit + Biopython +
prodigy-lig. PyMOL lives in its own env (`pymol`) because Schrödinger's
binary build does not coexist easily with PyTorch.

```bash
# 1. Activate base, create envs
source /root/miniconda3/etc/profile.d/conda.sh

# 2. Build the diffdock env (already done; refer to install log)
#    conda create -n diffdock python=3.9 setuptools=69.5.1 pip
#    pip install torch==1.13.1+cu117 --extra-index-url https://download.pytorch.org/whl/cu117
#    pip install torch-{scatter,sparse,cluster,spline-conv} \
#                -f https://pytorch-geometric.com/whl/torch-1.13.1+cu117.html
#    pip install e3nn==0.5.1 fair-esm==2.0.0 networkx==2.8.4 pandas==1.5.1 \
#                pybind11==2.11.1 pytorch-lightning==1.9.5 rdkit==2022.03.3 \
#                scikit-learn==1.1.0 torch-geometric==2.2.0 torchmetrics==0.11.0
#    pip install /root/dllogger
#    pip install --no-build-isolation /root/openfold
#    pip install dm-tree ml-collections
#    pip install prody==2.4.1                                # newer = compiles
#    pip install prodigy-lig prodigy-prot biopython
#
# 3. Build the pymol env
#    conda create -n pymol -c conda-forge -y pymol-open-source rdkit
```

## 1. Build the DiffDock CSV

```bash
source /root/miniconda3/etc/profile.d/conda.sh && conda activate diffdock
python /mnt/ShenYuchen/scripts/01_make_diffdock_csv.py
# -> /mnt/ShenYuchen/diffdock_input.csv  (1560 rows)
```

## 2. Run DiffDock

```bash
cd /root/DiffDock
nohup python -m inference \
    --config default_inference_args.yaml \
    --protein_ligand_csv /mnt/ShenYuchen/diffdock_input.csv \
    --out_dir /mnt/ShenYuchen/diffdock_results \
    --samples_per_complex 10 \
    > /mnt/ShenYuchen/diffdock_run.log 2>&1 &
```

For 1560 pairs on a single A16 (~30 s/pair) this is ~13 h. To split
across both A16 GPUs:

```bash
# Split CSV into halves
head -n 781 /mnt/ShenYuchen/diffdock_input.csv  > /tmp/csv_a.csv
(head -n 1 /mnt/ShenYuchen/diffdock_input.csv;
 tail -n +782 /mnt/ShenYuchen/diffdock_input.csv) > /tmp/csv_b.csv

cd /root/DiffDock
CUDA_VISIBLE_DEVICES=0 nohup python -m inference \
    --config default_inference_args.yaml --protein_ligand_csv /tmp/csv_a.csv \
    --out_dir /mnt/ShenYuchen/diffdock_results > /tmp/run_a.log 2>&1 &
CUDA_VISIBLE_DEVICES=1 nohup python -m inference \
    --config default_inference_args.yaml --protein_ligand_csv /tmp/csv_b.csv \
    --out_dir /mnt/ShenYuchen/diffdock_results > /tmp/run_b.log 2>&1 &
```

## 3. Combine pose + protein with PyMOL

```bash
conda activate pymol
pymol -cq /mnt/ShenYuchen/scripts/02_combine_pymol.py
# -> /mnt/ShenYuchen/complexes/<protein>__<ligand>.pdb
```

## 4. PRODIGY-LIG binding affinity

```bash
conda activate diffdock         # has prodigy-lig
python /mnt/ShenYuchen/scripts/03_run_prodigy.py
# -> /mnt/ShenYuchen/prodigy_results.csv
```

## 5. Distance between catalytic Cys and aldehyde C

```bash
python /mnt/ShenYuchen/scripts/04_distances.py
# -> /mnt/ShenYuchen/distances.csv
```

## 6. Aggregate + plot

```bash
python /mnt/ShenYuchen/scripts/05_analysis.py
# -> /mnt/ShenYuchen/master_results.csv
# -> /mnt/ShenYuchen/fig1_affinity_heatmap.png
# -> /mnt/ShenYuchen/fig2_dg_vs_distance.png
```
