# DiffDock Web UI Practice Notes

Date: 2026-05-10

Task: Use the DiffDock webpage user interface to practise with examples and inspect docking output from different viewing angles.

URL used: https://build.nvidia.com/mit/diffdock

## Examples Checked

1. Built-in default example:
   - Protein: `mpro_chain_OpenFold_processed_ms.pdb`
   - Ligand: `Ensitrelvir_analog_3d.sdf`
   - Settings shown by the UI: 20 generated poses, 18 diffusion steps, 20 diffusion time divisions.
   - Output inspection: the viewer opened with `View All Poses` selected and ranked poses 1-20 listed. Observed score range: rank 1 `-0.614`, rank 10 `-1.885`, rank 20 `-3.136`.

2. Built-in serotonin receptor example:
   - Protein: `htr2a_hunam_ESMFold_P28223_processed_ms.pdb`
   - Ligand: `Lisuride_analog_3d.sdf`
   - Settings used: 20 generated poses, 18 diffusion steps, 20 diffusion time divisions.
   - Output inspection: the submitted run completed and returned ranked poses 1-20. Observed score range: rank 1 `-0.495`, rank 10 `-1.271`, rank 20 `-3.692`.

The example selector also listed a cereblon/lenalidomide analog case (`crbn_human_ESMFold_Q96SW2_processed_ms.pdb` with `Lenalidomide_analog_3d.sdf`), confirming that the web UI supports multiple predefined protein-ligand examples.

## Viewer Observations

- `View All Poses` overlays all predicted ligand poses, which is useful for judging pose spread and pocket-level uncertainty.
- Selecting a single rank isolates that predicted ligand pose in the protein viewer.
- Rank 1 is the highest-confidence pose in the returned list. Later ranks have lower confidence scores and are useful for checking whether the ligand remains in the same pocket or moves to alternative local orientations.
- Rotating the 3D viewer changes interpretation materially: in the default orientation the protein fold and ligand cluster are visible, while a rotated view makes it easier to judge whether the ligand sits inside a cleft rather than on the protein surface.

## Saved Evidence

- `diffdock_serotonin_rank1_reset_view.png`
- `diffdock_serotonin_rank1_rotated_view.png`
- `diffdock_serotonin_rank10_rotated_view.png`

