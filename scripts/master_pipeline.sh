#!/usr/bin/env bash
# Master pipeline runner. Launched under nohup; needs no further attention.
# Runs DiffDock on both GPUs in parallel, then combine -> prodigy -> distance
# -> analysis. Writes everything to /root/ShenYuchen/ and copies the small
# artefacts to /mnt/ShenYuchen/ at the end so they survive a reboot.

set -u
ROOT=/root/ShenYuchen
LOG=$ROOT/master.log
exec >>"$LOG" 2>&1
date

export PATH=/root/miniconda3/envs/diffdock/bin:$PATH
export TORCH_HOME=/root/.cache/torch

cd /root/DiffDock

echo "[$(date)] launching DiffDock GPU 0"
CUDA_VISIBLE_DEVICES=0 python -u -m inference \
    --config default_inference_args.yaml \
    --protein_ligand_csv "$ROOT/csv_a.csv" \
    --out_dir "$ROOT/diffdock_results" \
    > "$ROOT/run_a.log" 2>&1 &
PID_A=$!

echo "[$(date)] launching DiffDock GPU 1"
CUDA_VISIBLE_DEVICES=1 python -u -m inference \
    --config default_inference_args.yaml \
    --protein_ligand_csv "$ROOT/csv_b.csv" \
    --out_dir "$ROOT/diffdock_results" \
    > "$ROOT/run_b.log" 2>&1 &
PID_B=$!

echo "[$(date)] waiting for PIDs $PID_A and $PID_B"
wait $PID_A
echo "[$(date)] PID_A finished with exit $?"
wait $PID_B
echo "[$(date)] PID_B finished with exit $?"

cd "$ROOT"

echo "[$(date)] combining poses"
python "$ROOT/scripts/02_combine.py"

echo "[$(date)] running PRODIGY-LIG"
python "$ROOT/scripts/03_run_prodigy.py"

echo "[$(date)] computing distances"
python "$ROOT/scripts/04_distances.py"

echo "[$(date)] aggregating + plotting"
python "$ROOT/scripts/05_analysis.py"

echo "[$(date)] filling report"
pip install --quiet scipy tabulate >/dev/null 2>&1 || true
python "$ROOT/scripts/06_finalize_report.py" || true

# Copy small artefacts to /mnt for persistence (skip the 1560-pair raw poses)
echo "[$(date)] syncing artefacts to /mnt/ShenYuchen"
mkdir -p /mnt/ShenYuchen/results
cp -f "$ROOT"/*.csv "$ROOT"/*.png /mnt/ShenYuchen/results/ 2>/dev/null || true
cp -rf "$ROOT/scripts" "$ROOT/manual" "$ROOT/report" "$ROOT/notebook.ipynb" \
       "$ROOT/README_server.md" /mnt/ShenYuchen/ 2>/dev/null || true

echo "[$(date)] DONE"
