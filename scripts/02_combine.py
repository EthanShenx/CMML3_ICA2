"""Combine the rank1 DiffDock pose with its protein into a single PDB.

Pure-Python implementation (no PyMOL needed). The output PDB has:
  * the protein atoms on chain A
  * the rank-1 ligand atoms as HETATM, chain L, residue UNK, resseq 1
in the format expected by PRODIGY-LIG (`-c A,L  -y L:UNK`).

We do *not* depend on PyMOL because the open-source build is not installable
through the SJTU conda mirror in our environment.
"""
import re
import sys
from pathlib import Path

from rdkit import Chem

ROOT = Path("/root/ShenYuchen")
PROT_DIR = ROOT / "Data" / "proteinPDB"
DOCK_DIR = ROOT / "diffdock_results"
OUT_DIR = ROOT / "complexes"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def read_protein_chainA(pdb_path: Path):
    """Return ATOM/HETATM lines of the protein, all on chain A.

    Drops ANISOU lines (PRODIGY-LIG ignores them anyway and they bloat the file).
    """
    out = []
    with pdb_path.open() as fh:
        for line in fh:
            tag = line[:6]
            if tag == "ANISOU":
                continue
            if tag in ("ATOM  ", "HETATM"):
                # rewrite chain ID at column 22 (index 21) to 'A'
                line = line[:21] + "A" + line[22:]
                out.append(line.rstrip("\n"))
    return out


def ligand_pdb_lines(sdf_path: Path, atom_offset: int):
    """Return PDB HETATM lines for the rank1 SDF, chain L, resname UNK."""
    mol = Chem.MolFromMolFile(str(sdf_path), removeHs=False)
    if mol is None:
        raise RuntimeError(f"could not parse {sdf_path}")
    conf = mol.GetConformer()
    out = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        elem = atom.GetSymbol()
        # PDB uses right-justified 2-char element symbol in cols 77-78
        atom_name = f"{elem:<2s}{idx+1:<2d}"[:4]
        line = (
            "HETATM"
            f"{atom_offset + idx + 1:>5d} "
            f"{atom_name:<4s} UNK L"
            f"{1:>4d}    "
            f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
            f"  1.00  0.00          "
            f"{elem:>2s}"
        )
        out.append(line)
    return out


def combine(pair_name: str):
    pair_dir = DOCK_DIR / pair_name
    out_pdb = OUT_DIR / f"{pair_name}.pdb"
    if out_pdb.exists():
        return "skip"
    prot_id = pair_name.split("__", 1)[0]
    prot_path = PROT_DIR / f"{prot_id}.pdb"
    rank1 = sorted(pair_dir.glob("rank1*.sdf"))
    if not rank1:
        return "no rank1"
    prot_lines = read_protein_chainA(prot_path)
    # last serial number on chain A (column 7-11)
    last_serial = int(prot_lines[-1][6:11])
    lig_lines = ligand_pdb_lines(rank1[0], atom_offset=last_serial)
    with out_pdb.open("w") as fh:
        fh.write("\n".join(prot_lines))
        fh.write("\nTER\n")
        fh.write("\n".join(lig_lines))
        fh.write("\nEND\n")
    return "ok"


def main():
    pairs = sorted(d for d in DOCK_DIR.iterdir() if d.is_dir())
    print(f"{len(pairs)} pair directories", flush=True)
    counts = {"ok": 0, "skip": 0, "no rank1": 0, "fail": 0}
    for d in pairs:
        try:
            counts[combine(d.name)] += 1
        except Exception as e:
            counts["fail"] += 1
            print(f"FAIL {d.name}: {e}", flush=True)
        if (counts["ok"] + counts["skip"]) % 100 == 0:
            print(f"  progress: {counts}", flush=True)
    print(f"done: {counts}")


if __name__ == "__main__":
    main()
