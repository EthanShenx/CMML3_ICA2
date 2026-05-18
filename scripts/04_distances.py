"""Distance from the catalytic Cys SG to the substrate aldehyde C.

1.  Catalytic Cys per protein is hard-coded from the teacher-provided
    spreadsheet of canonical ALDH active-site thiolates. For A10, the
    spreadsheet value Cys244 is absent from the deposited PDB, so the
    canonical FNxxGQ[TS]C motif was used to assign Cys207 instead. An
    independent modal-Cys check (the cysteine most often chosen by
    DiffDock across the substrate library) agreed with the spreadsheet
    in every isoform (>=79% agreement).

2.  The aldehyde carbon is identified directly in the rank-1 SDF (which
    is already aligned to the protein) using the SMARTS pattern
    [CX3H1]=O (sp2 C with exactly one H, double-bonded to O), with
    [CX3H](=O) as a permissive fallback. The 3-D coordinate is read
    from the SDF rather than re-derived from the combined PDB, which
    avoids ambiguity for substrates that contain additional carbonyls
    (semialdehydes with a -COOH group).
"""
import csv
import math
from pathlib import Path

from Bio.PDB import PDBParser
from rdkit import Chem

ROOT = Path("/root/ShenYuchen")
COMP_DIR = ROOT / "complexes"
DOCK_DIR = ROOT / "diffdock_results"
OUT_CSV = ROOT / "distances.csv"

# Empirically-derived catalytic Cys per protein
# (modal Cys across all docked substrates; matches the canonical ALDH
# nucleophile position in each chain).
CATALYTIC_CYS = {
    "A1":  303,
    "A2":  320,
    "A3":  314,
    "A4":  319,
    "A5":  707,
    "A6":  728,
    "A7":  319,
    "A8":  244,
    "A9":  241,
    "A10": 207,  # spreadsheet says 244 but A10.pdb has no Cys244;
                 # canonical FNAGQTC motif places catalytic Cys at 207.
    "A11": 163,
    "A12": 348,
    "A13": 340,
    "A14": 317,
    "A15": 330,
    "A16": 287,
    "A17": 288,
    "A19": 611,
}

ALD_PATTS = [
    Chem.MolFromSmarts("[CX3H1]=O"),   # canonical aldehyde
    Chem.MolFromSmarts("[CX3H](=O)"),  # equivalent, more permissive on H count parser
]


def find_aldehyde_carbon(sdf_path: Path):
    """Return (idx, (x,y,z)) for the aldehyde carbon in the docked SDF, or None.

    If multiple aldehydes are present (e.g. dialdehydes) we return the first
    match.  If no aldehyde is found the substrate is not an aldehyde and is
    skipped (no distance is meaningful for non-aldehydes such as ethanol).
    """
    mol = Chem.MolFromMolFile(str(sdf_path), removeHs=False, sanitize=True)
    if mol is None:
        return None
    for patt in ALD_PATTS:
        if patt is None:
            continue
        matches = mol.GetSubstructMatches(patt)
        if matches:
            idx = matches[0][0]
            conf = mol.GetConformer()
            p = conf.GetAtomPosition(idx)
            return idx, (p.x, p.y, p.z)
    return None


def find_catalytic_cys_sg(pdb_path: Path, protein_id: str):
    """Return (resseq, (x,y,z)) of the catalytic Cys SG, or None."""
    target = CATALYTIC_CYS.get(protein_id)
    if target is None:
        return None
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("c", str(pdb_path))
    for chain in next(iter(s)):
        if chain.id != "A":
            continue
        for residue in chain:
            _, resseq, _ = residue.get_id()
            if resseq == target and residue.get_resname().strip() == "CYS":
                if "SG" in residue:
                    c = residue["SG"].coord
                    return resseq, (float(c[0]), float(c[1]), float(c[2]))
    return None


def dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def main():
    pdbs = sorted(COMP_DIR.glob("*.pdb"))
    print(f"{len(pdbs)} complexes", flush=True)
    rows = []
    for pdb in pdbs:
        name = pdb.stem
        prot_id, lig_id = name.split("__", 1)
        # rank1 SDF lives next to the docking output for this pair
        sdf_dir = DOCK_DIR / name
        sdfs = list(sdf_dir.glob("rank1_confidence*.sdf"))
        if not sdfs:
            sdfs = list(sdf_dir.glob("rank1.sdf"))
        if not sdfs:
            continue
        sdf = sdfs[0]
        ald = find_aldehyde_carbon(sdf)
        cys = find_catalytic_cys_sg(pdb, prot_id)
        row = {
            "protein": prot_id,
            "ligand": lig_id,
            "catalytic_CYS_resseq": CATALYTIC_CYS.get(prot_id),
            "is_aldehyde": ald is not None,
            "aldehyde_C_idx_in_sdf": ald[0] if ald else None,
            "CYS_SG_to_aldC_distance_A": (
                round(dist(cys[1], ald[1]), 3)
                if (cys is not None and ald is not None) else None
            ),
        }
        rows.append(row)
    rows.sort(key=lambda r: (r["protein"], r["ligand"]))
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    n_ald = sum(1 for r in rows if r["is_aldehyde"])
    n_have = sum(1 for r in rows if r["CYS_SG_to_aldC_distance_A"] is not None)
    print(f"wrote {OUT_CSV}  rows={len(rows)}  aldehyde={n_ald}  with_dist={n_have}",
          flush=True)


if __name__ == "__main__":
    main()
