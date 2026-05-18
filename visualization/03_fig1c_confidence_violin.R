# Fig 1 panel (c) — DiffDock rank-1 confidence violin per isoform (horizontal).
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("_palettes.R")

d <- read.csv(file.path(DATA2PLOT, "fig1c_confidence.csv"))
d$protein <- factor(d$protein, levels = PROTEINS)

med_order <- d %>%
  group_by(protein) %>%
  summarise(med = median(rank1_confidence, na.rm = TRUE)) %>%
  arrange(med) %>%
  pull(protein) %>%
  as.character()
d$protein <- factor(d$protein, levels = med_order)

p <- ggplot(d, aes(y = protein, x = rank1_confidence, fill = protein)) +
  geom_violin(scale = "width", trim = TRUE, colour = "grey40",
              linewidth = 0.3, alpha = 1) +
  geom_boxplot(width = 0.12, outlier.size = 0.4, fill = "white",
               colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = PROT_COLORS, guide = "none") +
  labs(x = "DiffDock rank-1 confidence", y = NULL,
       title = "Confidence per isoform",
       subtitle = "(c)") +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5))
save_panel(p, "03_fig1c_confidence_violin", width = 5.5, height = 4)
