"""Batch combine the rank1 DiffDock pose with its protein PDB.

For each pair: load the original protein, load the rank1 ligand SDF,
rename protein chain to A and ligand chain to L, residue UNK, then
write a combined PDB. Run inside a `pymol -cq` interpreter.
"""
import os
import re
import sys
import glob
from pathlib import Path

from pymol import cmd

ROOT     = Path("/mnt/ShenYuchen")
PROT_DIR = ROOT / "Data" / "proteinPDB"
DOCK_DIR = ROOT / "diffdock_results"
OUT_DIR  = ROOT / "complexes"
OUT_DIR.mkdir(parents=True, exist_ok=True)

pair_dirs = sorted(d for d in DOCK_DIR.iterdir() if d.is_dir())
print(f"{len(pair_dirs)} pair directories found", flush=True)

ok = fail = 0
for d in pair_dirs:
    name = d.name  # e.g. A1__1_formaldehyde
    out = OUT_DIR / f"{name}.pdb"
    if out.exists():
        ok += 1
        continue
    try:
        prot_id = name.split("__", 1)[0]
        prot_pdb = PROT_DIR / f"{prot_id}.pdb"
        rank1 = sorted(d.glob("rank1*.sdf"))
        if not rank1:
            raise FileNotFoundError("no rank1 sdf")
        cmd.reinitialize()
        cmd.load(str(prot_pdb), "prot")
        cmd.load(str(rank1[0]), "lig")
        # Force protein chain A
        cmd.alter("prot", "chain='A'")
        # Force ligand chain L and resname UNK (PRODIGY-LIG needs a 3-letter resname)
        cmd.alter("lig", "chain='L'")
        cmd.alter("lig", "resn='UNK'")
        cmd.alter("lig", "resi='1'")
        cmd.alter("lig", "type='HETATM'")
        cmd.sort()
        cmd.save(str(out), "prot or lig")
        ok += 1
    except Exception as e:
        fail += 1
        print(f"FAIL {name}: {e}", flush=True)

print(f"done ok={ok} fail={fail}", flush=True)
