# Fig 1 panel (e) — per-isoform ΔG violin (aldehydes only).
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("_palettes.R")

d <- read.csv(file.path(DATA2PLOT, "fig1e_dg_violin.csv"))
d <- d[!is.na(d$DG_prediction_kcalmol), ]

med_order <- d %>%
  group_by(protein) %>%
  summarise(med = median(DG_prediction_kcalmol)) %>%
  arrange(med) %>%
  pull(protein) %>% as.character()
d$protein <- factor(d$protein, levels = med_order)

p <- ggplot(d, aes(y = protein, x = DG_prediction_kcalmol, fill = protein)) +
  geom_violin(scale = "width", trim = TRUE, colour = "grey40",
              linewidth = 0.3, alpha = 0.85) +
  geom_boxplot(width = 0.12, outlier.size = 0.4, fill = "white",
               colour = "black", linewidth = 0.3) +
  geom_vline(xintercept = -7, linetype = "dashed", colour = "grey50",
             linewidth = 0.4) +
  scale_fill_manual(values = PROT_COLORS, guide = "none") +
  labs(x = expression(Delta * G~"(kcal mol"^{-1}*")"), y = NULL,
       title = "Affinity distributions\nper isoform",
       subtitle = "(e)") +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9))
save_panel(p, "05_fig1e_dg_violin", width = 5, height = 4)
