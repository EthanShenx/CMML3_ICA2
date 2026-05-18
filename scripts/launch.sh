#!/usr/bin/env bash
# Wrapper that ensures ESM2 weights are present (with resumable download)
# and then runs the master DiffDock pipeline.
LOG=/root/ShenYuchen/launch.log
exec >>"$LOG" 2>&1
set -u

ESM_DIR=/root/.cache/torch/hub/checkpoints
ESM=$ESM_DIR/esm2_t33_650M_UR50D.pt
ESM_REG=$ESM_DIR/esm2_t33_650M_UR50D-contact-regression.pt
URL=https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt
URL_REG=https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t33_650M_UR50D-contact-regression.pt
EXPECTED=2604537549   # bytes — verified via HEAD

mkdir -p "$ESM_DIR"

echo "[$(date)] launch.sh started"

# Resumable download with retries: parallel HTTP Range download (8 chunks)
download_until_complete () {
    local url=$1 out=$2 want=$3
    while :; do
        local have=0
        [ -f "$out" ] && have=$(stat -c %s "$out")
        if [ "$have" -ge "$want" ]; then
            echo "[$(date)] $out is complete ($have bytes)"
            return 0
        fi
        echo "[$(date)] downloading $url -> $out (have=$have, want=$want)"
        /root/miniconda3/envs/diffdock/bin/python /root/ShenYuchen/scripts/parallel_dl.py "$url" "$out" "$want" || sleep 10
    done
}

# Regression file is small and may already be present
[ ! -f "$ESM_REG" ] && wget -q -c -t 0 -O "$ESM_REG" "$URL_REG"
download_until_complete "$URL" "$ESM" "$EXPECTED"

bash /root/ShenYuchen/scripts/master_pipeline.sh
echo "[$(date)] launch.sh DONE"
