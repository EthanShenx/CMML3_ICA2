# Fig 2 panel (c2) — structural pose A12 × β-carotene ald., non-productive.
suppressPackageStartupMessages({
  library(ggplot2)
  library(png)
  library(grid)
})
source("_palettes.R")

img <- readPNG(file.path(FIG_OUT, "..", "pose_panel",
                         "panel_A12_nonproductive.png"))
g <- rasterGrob(img, interpolate = TRUE)

p <- ggplot() +
  annotation_custom(g, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  annotate("label", x = 0.03, y = 0.97,
           label = "Cys-SG ↔ aldehyde-C: 5.49 Å",
           hjust = 0, vjust = 1, size = 2.6, label.size = 0.3,
           fill = "white") +
  annotate("label", x = 0.97, y = 0.97,
           label = "non-productive",
           hjust = 1, vjust = 1, size = 2.6, label.size = 0.3,
           fill = "#C62828", colour = "white", fontface = "plain") +
  annotate("text", x = -0.02, y = 1.06, label = "(c2)",
           hjust = 0, vjust = 1, size = 4, fontface = "plain") +
  annotate("text", x = 0.5, y = 1.06,
           label = "A12 × β-carotene ald.", size = 3.2, fontface = "plain") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.08), expand = FALSE) +
  theme_void()
save_panel(p, "13_fig2c2_pose_A12", width = 6, height = 4)
