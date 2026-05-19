# Render A8 (productive) and A12 (non-productive) DiffDock rank-1 poses
# of substrate 31 (β-carotene-derived polyene aldehyde).

bg_color white
set ray_shadows, 0
set ambient, 0.4
set specular, 0.3
set ray_trace_mode, 1
set ray_trace_color, grey50
set cartoon_transparency, 0.55
set dash_radius, 0.12
set dash_gap, 0.35
set dash_color, black
set dash_width, 4
set label_size, 22
set label_color, black
set label_font_id, 7
set label_outline_color, white

python
def render(iso, cys, lig_color, outfile):
    cmd.delete("all")
    cmd.load(iso + ".pdb", "rec")
    cmd.load(iso + "_rank1.sdf", "lig")
    # Receptor
    cmd.hide("everything", "rec")
    cmd.show("cartoon", "rec")
    cmd.color("grey80", "rec and polymer")
    # Active-site Cys
    sel_cys = "rec and resi {} and resn CYS".format(cys)
    cmd.show("sticks", sel_cys)
    cmd.color("yellow", sel_cys + " and elem C")
    cmd.color("orange", sel_cys + " and elem S")
    # Ligand sticks
    cmd.show("sticks", "lig")
    cmd.color(lig_color, "lig and elem C")
    cmd.color("red", "lig and elem O")
    # Highlight atoms (small spheres, applied AFTER selection)
    cmd.select("sg_atom", sel_cys + " and name SG")
    cmd.select("ald_c_atom", "lig and rank 20")
    cmd.show("spheres", "sg_atom")
    cmd.show("spheres", "ald_c_atom")
    cmd.set("sphere_scale", 0.45, "sg_atom")
    cmd.set("sphere_scale", 0.45, "ald_c_atom")
    cmd.color("orange", "sg_atom")
    cmd.color("firebrick", "ald_c_atom")
    # Distance dashes + label
    d = cmd.get_distance("sg_atom", "ald_c_atom")
    cmd.distance("dist_obj", "sg_atom", "ald_c_atom")
    cmd.hide("labels", "dist_obj")
    # Orient
    cmd.zoom("sg_atom or ald_c_atom or " + sel_cys, buffer=10)
    if iso == "A8":
        cmd.turn("y", -25); cmd.turn("x", 5)
    else:
        cmd.turn("y", 15); cmd.turn("x", -10)
    cmd.ray(1800, 1400)
    cmd.png(outfile, dpi=300)
    print("DIST {} {:.2f}".format(iso, d))
    return d

dA8  = render("A8",  244, "salmon", "panel_A8_productive.png")
dA12 = render("A12", 348, "slate",  "panel_A12_nonproductive.png")
print("FINAL A8={:.2f} A   A12={:.2f} A".format(dA8, dA12))
python end
