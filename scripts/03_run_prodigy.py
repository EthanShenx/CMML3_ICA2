"""Batch run prodigy_lig over every combined PDB and parse DG into a CSV."""
import csv
import re
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

ROOT = Path("/root/ShenYuchen")
COMP_DIR = ROOT / "complexes"
DOCK_DIR = ROOT / "diffdock_results"
OUT_CSV = ROOT / "prodigy_results.csv"

DG_RE = re.compile(r"-?\d+\.\d+")


def confidence_from_pair_dir(d: Path):
    """Extract DiffDock rank-1 confidence from the dedicated SDF filename."""
    for sdf in d.glob("rank1_confidence*.sdf"):
        m = re.search(r"rank1_confidence(-?\d+\.\d+)", sdf.name)
        if m:
            return float(m.group(1))
    return None


def run_one(pdb: Path):
    name = pdb.stem
    try:
        cp = subprocess.run(
            ["prodigy_lig", "-c", "A", "L:UNK", "-i", str(pdb)],
            capture_output=True, text=True, timeout=120,
        )
        out = cp.stdout + cp.stderr
        # Output format is two lines:
        #   "Job name\tDGprediction (low refinement) (Kcal/mol)"
        #   "<job-name>\t<dg-value>"
        # The data row starts with the input file's stem (== `name`).
        dg = None
        for line in out.splitlines():
            line_strip = line.strip()
            if line_strip.startswith(name):
                m = DG_RE.search(line_strip.split(maxsplit=1)[-1])
                if m:
                    dg = float(m.group(0))
                    break
        return name, dg, cp.returncode, out
    except Exception as e:
        return name, None, -1, str(e)


def main():
    pdbs = sorted(COMP_DIR.glob("*.pdb"))
    print(f"{len(pdbs)} complexes to score", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        for fut in as_completed([ex.submit(run_one, p) for p in pdbs]):
            name, dg, rc, _out = fut.result()
            prot, lig = name.split("__", 1)
            conf = confidence_from_pair_dir(DOCK_DIR / name)
            rows.append({"protein": prot, "ligand": lig, "rank1_confidence": conf,
                         "DG_prediction_kcalmol": dg, "prodigy_rc": rc})
            if len(rows) % 50 == 0:
                print(f"scored {len(rows)}/{len(pdbs)}", flush=True)
    rows.sort(key=lambda r: (r["protein"], r["ligand"]))
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT_CSV}", flush=True)


if __name__ == "__main__":
    main()
