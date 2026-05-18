# Fig 1 panel (a) — pipeline workflow schematic.
# Reproduces the four FancyBboxPatch nodes + arrows from 07_make_figures.py.
suppressPackageStartupMessages({
  library(ggplot2)
})
source("_palettes.R")

steps <- data.frame(
  x      = c(1.0, 4.0, 7.0, 9.5),
  y      = c(3.2, 3.2, 3.2, 3.2),
  fill   = c("#A5D6A7", "#90CAF9", "#FFCC80", "#F48FB1"),
  label  = c("Protein PDB\n+ Ligand SDF\n(420 pairs)",
             "DiffDock-L 1.1\n10 poses/pair\n2× A16 GPU",
             "PRODIGY-LIG\nΔG scoring",
             "Cys-SG ↔\naldehyde-C\ndistance")
)
arrows <- data.frame(
  x    = steps$x[-nrow(steps)] + 0.9,
  xend = steps$x[-1]            - 0.9,
  y    = 3.28, yend = 3.28
)
p <- ggplot() +
  geom_rect(data = steps,
            aes(xmin = x - 0.9, xmax = x + 0.9,
                ymin = y - 0.7, ymax = y + 0.85, fill = I(fill)),
            colour = "#555", linewidth = 0.4) +
  geom_text(data = steps, aes(x = x, y = y + 0.08, label = label),
            fontface = "plain", size = 2.7, lineheight = 0.9) +
  geom_segment(data = arrows,
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               colour = "#333", linewidth = 0.7) +
  annotate("label", x = 5.25, y = 0.6,
           label = "master_results.csv\n(ΔG + confidence + distance)",
           size = 2.6, fontface = "italic", colour = "#555",
           fill = "#FFF9C4", label.size = 0.3, label.r = unit(0.15, "lines")) +
  annotate("text", x = 5.25, y = 0.15,
           label = "→ affinity map + geometry filter",
           size = 2.5, colour = "#666") +
  annotate("text", x = -0.4, y = 4.3, label = "(a)",
           fontface = "plain", size = 4, hjust = 0) +
  coord_cartesian(xlim = c(-0.4, 10.6), ylim = c(0, 4.5), expand = FALSE) +
  theme_void()
save_panel(p, "01_fig1a_workflow", width = 6.5, height = 3)
