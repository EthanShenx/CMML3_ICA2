# Fig 1 panel (h) — inter-isoform Pearson correlation heatmap.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
source("_palettes.R")

m <- read.csv(file.path(DATA2PLOT, "fig1h_corr.csv"), check.names = FALSE)
names(m)[1] <- "row"
long <- m %>% pivot_longer(-row, names_to = "col", values_to = "r")
long$row <- factor(long$row, levels = PROTEINS)
long$col <- factor(long$col, levels = PROTEINS)

p <- ggplot(long, aes(x = col, y = row, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.2) +
  scale_fill_distiller(palette = "RdYlBu", direction = -1,
                       limits = c(0, 1), name = "Pearson r") +
  coord_fixed() +
  labs(x = NULL, y = NULL, title = "Inter-isoform\ncorrelation",
       subtitle = "(h)") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9),
        axis.text     = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"))
save_panel(p, "08_fig1h_corr_heatmap", width = 4.5, height = 4)
