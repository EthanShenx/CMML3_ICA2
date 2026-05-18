#!/usr/bin/env bash
# End-to-end pipeline. Run from /mnt/ShenYuchen on the GPU node.
set -euo pipefail
source /root/miniconda3/etc/profile.d/conda.sh

ROOT=/root/ShenYuchen
SCRIPTS=$ROOT/scripts

echo "=== 1. Build DiffDock CSV ==="
conda activate diffdock
python "$SCRIPTS/01_make_diffdock_csv.py"

echo "=== 2. DiffDock inference (split across 2 GPUs) ==="
mkdir -p "$ROOT/diffdock_results"
total=$(($(wc -l < "$ROOT/diffdock_input.csv") - 1))
half=$(( (total + 1) / 2 ))
head -n 1 "$ROOT/diffdock_input.csv" > /tmp/csv_a.csv
head -n $((half + 1)) "$ROOT/diffdock_input.csv" | tail -n +2 >> /tmp/csv_a.csv
head -n 1 "$ROOT/diffdock_input.csv" > /tmp/csv_b.csv
tail -n +$((half + 2)) "$ROOT/diffdock_input.csv" >> /tmp/csv_b.csv
cd /root/DiffDock
CUDA_VISIBLE_DEVICES=0 python -m inference \
    --config default_inference_args.yaml --protein_ligand_csv /tmp/csv_a.csv \
    --out_dir "$ROOT/diffdock_results" > "$ROOT/run_a.log" 2>&1 &
PID_A=$!
CUDA_VISIBLE_DEVICES=1 python -m inference \
    --config default_inference_args.yaml --protein_ligand_csv /tmp/csv_b.csv \
    --out_dir "$ROOT/diffdock_results" > "$ROOT/run_b.log" 2>&1 &
PID_B=$!
wait $PID_A $PID_B
echo "DiffDock done"

echo "=== 3. Combine pose+protein (Biopython/RDKit) ==="
python "$SCRIPTS/02_combine.py"

echo "=== 4. PRODIGY-LIG ==="
python "$SCRIPTS/03_run_prodigy.py"

echo "=== 5. Distances ==="
python "$SCRIPTS/04_distances.py"

echo "=== 6. Analysis & plots ==="
python "$SCRIPTS/05_analysis.py"

echo "=== ALL DONE ==="
