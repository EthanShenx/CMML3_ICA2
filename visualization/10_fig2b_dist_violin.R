# Fig 2 panel (b) — per-isoform Cys-SG distance violin.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("_palettes.R")

d <- read.csv(file.path(DATA2PLOT, "fig2b_dist_violin.csv"))
med_order <- d %>%
  group_by(protein) %>%
  summarise(med = median(distance)) %>%
  arrange(med) %>%
  pull(protein) %>% as.character()
d$protein <- factor(d$protein, levels = med_order)

p <- ggplot(d, aes(y = protein, x = distance, fill = protein)) +
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = 0, xmax = 4,
           fill = "#A8D5A2", alpha = 0.12) +
  geom_violin(scale = "width", trim = TRUE, colour = "grey40",
              linewidth = 0.3, alpha = 0.85) +
  geom_boxplot(width = 0.12, outlier.size = 0.4, fill = "white",
               colour = "black", linewidth = 0.3) +
  geom_vline(xintercept = 4.0, linetype = "dashed", colour = "#D62728",
             linewidth = 0.5) +
  scale_fill_manual(values = PROT_COLORS, guide = "none") +
  labs(y = "ALDH isoform", x = "Distance (Å)",
       title = "Catalytic-Cys distance\nper isoform",
       subtitle = "(b)") +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9))
save_panel(p, "10_fig2b_dist_violin", width = 5, height = 4)
