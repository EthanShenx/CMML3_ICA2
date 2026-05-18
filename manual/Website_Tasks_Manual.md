# Website Tasks — Step-by-step Manual (Task 1 of the ICA)

This document covers everything that must be done by *you* on a web browser
(no SSH involved). It corresponds to Section 1 of the project sheet:
"*Use DiffDock webpage user interface to practice using some examples,
observing different angles of this analysis*."

The downstream batch processing on the GPU server (Section 2) is fully
automated by the scripts in `ShenYuchen/scripts/` and is documented in
`ShenYuchen/README_server.md`.

> Suggested example pair for the walkthrough below (used as the "Figure 1"
> case study in the report):
>
> **Protein:** `Data/proteinPDB/A2.pdb` (mitochondrial ALDH2)
> **Ligand:**  `Data/substrate/2_acetaldehyde.sdf` (its physiological substrate)
>
> Feel free to repeat with any other pair.

---

## 1. DiffDock — NVIDIA build playground

URL: <https://build.nvidia.com/mit/diffdock>

1. Open the page. Sign in with the NVIDIA account if prompted (free).
2. In the **Molecule** box click `Upload New File` and pick
   `Data/substrate/2_acetaldehyde.sdf`.
3. In the **Target Protein** box click `Upload New File` and pick
   `Data/proteinPDB/A2.pdb`.
4. Set the parameters:
   * **Generated Poses:** 10 (more poses = more diverse rankings; 10 is the default for the local CLI we use later).
   * **Diffusion Steps:** 18.
   * **Diffusion Time Divisions:** 20.
5. Click **Run**. Wait ~1–3 min.
6. When the right pane finishes loading, tick **View All Poses**, scroll
   through ranks 1–10 and note that:
   * each pose has a *confidence score* (higher = better),
   * pose 1 should sit deep in the substrate-binding cleft of the protein.
7. Click **Download** (top-right of the output pane). Save the zip,
   unzip into a folder e.g. `webpage_test/A2_acetaldehyde/`.
   The folder will contain `target_protein.pdb`, `rank{1..10}_*.sdf`,
   `rank{1..10}_*.pdb` and `pose_confidence.txt`.
8. Take a screenshot of the 3D viewer with rank-1 selected — this can be
   used as a supplementary figure.

## 2. PyMOL — combining protein + rank-1 pose

(Local PyMOL on your laptop is fine; the GUI is easier than the terminal version
for this manual step.)

1. Launch PyMOL.
2. `File → Open…` → `target_protein.pdb`.
3. `File → Open…` → `rank1_confidence-*.sdf`.
4. In the command bar at the top, run, in order:
   ```
   alter target_protein, chain='A'
   alter rank1_confidence_*, chain='L'
   alter rank1_confidence_*, resn='UNK'
   alter rank1_confidence_*, resi='1'
   alter rank1_confidence_*, type='HETATM'
   sort
   ```
   (the wildcard `rank1_confidence_*` will match the actual loaded object.)
5. `File → Export Molecule…` → make sure both objects are selected → save as
   `A2_acetaldehyde_complex.pdb`.
6. Optional rendering for figures:
   ```
   bg_color white
   hide everything
   show cartoon, chain A
   show sticks, chain L
   show spheres, name SG and resn CYS within 6 of chain L
   color cyan, chain A
   color yellow, chain L
   set ray_shadow, 0
   ray 1200, 900
   png A2_acetaldehyde_complex.png, dpi=300
   ```

## 3. Prodigy-LIG — docking free energy on the web

URL: <https://wenmr.science.uu.nl/prodigy/lig>

1. Switch to the **PRODIGY-LIGAND** tab.
2. Under *Structures*, click *选择文件* and upload
   `A2_acetaldehyde_complex.pdb`.
3. **Protein Chain ID:** `A`
4. **Ligand ID:** `L:UNK`  (chain `L`, three-letter residue `UNK`).
5. Leave electrostatic energy blank.
6. **Job ID:** anything memorable, e.g. `A2_acetaldehyde_v1`.
7. Click **Submit Prodigy-Ligand**.
8. The page auto-refreshes; after the queue clears, you will see a result
   table. Note the **ΔG_noelec (kcal/mol)** column — this is the value we
   compare to in the report.

## 4. (Optional) inspecting key distances visually

In PyMOL after step 2:
```
sele cyss, name SG and resn CYS within 8 of chain L
sele aldC, chain L and elem C and (chain L within 2.0 of (chain L and elem O))
distance d_sg_c, cyss, aldC, mode=2
```
The shortest dashed line is the catalytic-Cys ↔ aldehyde-carbon distance
that the batch script (`scripts/04_distances.py`) computes for every pair.

---

## What to record

* one PNG of the rank-1 pose for the chosen example pair (used as Fig. 1a),
* one screenshot of the Prodigy-LIG result for the same pair (used as
  evidence that the web tool agrees with the cluster batch in §2 of the report),
* the Prodigy ΔG_noelec value for the chosen pair (will be cited in-text).

These three artefacts come from the manual website work. Everything else in
the report comes from the batch pipeline that runs on the GPU server.
