# Fig 2 panel (c1) — structural pose A8 × β-carotene ald., productive.
suppressPackageStartupMessages({
  library(ggplot2)
  library(png)
  library(grid)
})
source("_palettes.R")

img <- readPNG(file.path(FIG_OUT, "..", "pose_panel",
                         "panel_A8_productive.png"))
g <- rasterGrob(img, interpolate = TRUE)

p <- ggplot() +
  annotation_custom(g, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  annotate("label", x = 0.03, y = 0.97,
           label = "Cys-SG ↔ aldehyde-C: 2.27 Å",
           hjust = 0, vjust = 1, size = 2.6, label.size = 0.3,
           fill = "white") +
  annotate("label", x = 0.97, y = 0.97,
           label = "productive",
           hjust = 1, vjust = 1, size = 2.6, label.size = 0.3,
           fill = "#2E7D32", colour = "white", fontface = "plain") +
  annotate("text", x = -0.02, y = 1.06, label = "(c1)",
           hjust = 0, vjust = 1, size = 4, fontface = "plain") +
  annotate("text", x = 0.5, y = 1.06,
           label = "A8 × β-carotene ald.", size = 3.2, fontface = "plain") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.08), expand = FALSE) +
  theme_void()
save_panel(p, "12_fig2c1_pose_A8", width = 6, height = 4)
