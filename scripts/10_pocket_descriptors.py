"""Compute simple active-site pocket descriptors for A11 and A12.

Used to back the A11-versus-A12 mechanism sentence in the main text with
one quantitative geometric descriptor each (size, hydrophobicity, depth).
"""
from pathlib import Path
import math
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parent.parent
PDB_DIR = ROOT.parent / "Data" / "proteinPDB"

CATALYTIC_CYS = {"A11": 163, "A12": 348}

HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY"}
POLAR_NEUTRAL = {"SER", "THR", "ASN", "GLN", "TYR", "CYS", "HIS"}
CHARGED = {"ASP", "GLU", "LYS", "ARG"}

POCKET_R = 12.0  # Angstrom shell around the catalytic Cys CA, covering the substrate channel

def descriptors(iso):
    pdb = PDB_DIR / f"{iso}.pdb"
    structure = PDBParser(QUIET=True).get_structure(iso, str(pdb))
    chain = next(structure[0].get_chains())
    cys = None
    for res in chain:
        if res.get_resname() == "CYS" and res.id[1] == CATALYTIC_CYS[iso]:
            cys = res
            break
    if cys is None or "SG" not in cys:
        raise RuntimeError(f"{iso}: catalytic Cys{CATALYTIC_CYS[iso]} SG not found")
    sg = cys["SG"].coord
    ca_ref = cys["CA"].coord

    pocket = []
    for res in chain:
        if res.id[0] != " " or "CA" not in res:
            continue
        d = math.dist(res["CA"].coord, ca_ref)
        if d <= POCKET_R:
            pocket.append((res.get_resname(), res.id[1], d))

    nh = sum(1 for r, *_ in pocket if r in HYDROPHOBIC)
    npo = sum(1 for r, *_ in pocket if r in POLAR_NEUTRAL)
    nch = sum(1 for r, *_ in pocket if r in CHARGED)
    total = len(pocket)

    # Pocket "depth" approximated as the mean Cys-SG to nearest-protein-edge
    # distance, computed as the gyration radius of the pocket Cα cloud.
    if total >= 2:
        xs = [res_ca for res in chain if res.id[0] == " " and "CA" in res
              and math.dist(res["CA"].coord, ca_ref) <= POCKET_R
              for res_ca in [res["CA"].coord]]
        cx = sum(p[0] for p in xs) / total
        cy = sum(p[1] for p in xs) / total
        cz = sum(p[2] for p in xs) / total
        rg = math.sqrt(sum((p[0] - cx) ** 2 + (p[1] - cy) ** 2 + (p[2] - cz) ** 2
                           for p in xs) / total)
    else:
        rg = float("nan")

    return {
        "isoform": iso,
        "cys": CATALYTIC_CYS[iso],
        "pocket_residues": total,
        "hydrophobic": nh,
        "polar_neutral": npo,
        "charged": nch,
        "hydrophobic_pct": 100.0 * nh / total if total else float("nan"),
        "Rg_A": rg,
    }


if __name__ == "__main__":
    for iso in ("A11", "A12"):
        d = descriptors(iso)
        print(f"{d['isoform']}  Cys{d['cys']}  pocket residues within "
              f"{POCKET_R:.0f} A of CA: {d['pocket_residues']}  "
              f"hydrophobic={d['hydrophobic']} ({d['hydrophobic_pct']:.0f}%)  "
              f"polar={d['polar_neutral']}  charged={d['charged']}  "
              f"Rg={d['Rg_A']:.2f} A")
