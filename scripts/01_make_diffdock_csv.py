"""Build the DiffDock batch CSV: every protein x every ligand pair."""
import csv
import re
from pathlib import Path

ROOT = Path("/root/ShenYuchen")
PROT_DIR = ROOT / "Data" / "proteinPDB"
LIG_DIR = ROOT / "Data" / "substrate"
OUT_CSV = ROOT / "diffdock_input.csv"

def natkey(p: Path):
    m = re.match(r"^([A-Za-z]+)(\d+)", p.stem)
    return (m.group(1), int(m.group(2))) if m else (p.stem, 0)

proteins = sorted(PROT_DIR.glob("*.pdb"), key=natkey)
ligands = sorted(
    LIG_DIR.glob("*.sdf"),
    key=lambda p: int(p.stem.split("_", 1)[0]),
)

print(f"{len(proteins)} proteins, {len(ligands)} ligands -> {len(proteins)*len(ligands)} pairs")

with OUT_CSV.open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["complex_name", "protein_path", "ligand_description", "protein_sequence"])
    for prot in proteins:
        pid = prot.stem
        for lig in ligands:
            lid = lig.stem
            name = f"{pid}__{lid}"
            w.writerow([name, str(prot), str(lig), ""])

print(f"Wrote {OUT_CSV}")
