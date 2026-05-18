# SI Fig S1 panel (c) — confidence-affinity correlation per substrate class.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("_palettes.R")

d  <- read.csv(file.path(DATA2PLOT, "figS1c_conf_dg.csv"),
               fileEncoding = "UTF-8", check.names = FALSE)
sp <- read.csv(file.path(DATA2PLOT, "figS1c_spearman.csv"),
               fileEncoding = "UTF-8", check.names = FALSE)
d$category <- factor(d$category, levels = CAT_ORDER)

# Legend label per category: "{cat} (ρ=.xx)"
lab <- setNames(rep(NA_character_, length(CAT_ORDER)), CAT_ORDER)
for (i in seq_len(nrow(sp))) {
  lab[sp$category[i]] <- sprintf("%s (ρ=%.2f)", sp$category[i], sp$rho[i])
}
lab[is.na(lab)] <- names(lab[is.na(lab)])

p <- ggplot(d, aes(x = confidence, y = DG, colour = category)) +
  geom_point(size = 1.3, alpha = 0.7) +
  scale_colour_manual(values = CAT_COLORS, labels = lab, name = NULL,
                      drop = FALSE) +
  labs(x = "DiffDock rank-1 confidence",
       y = expression(Delta * G~"(kcal mol"^{-1}*")"),
       title = "Confidence–affinity correlation\nby substrate class",
       subtitle = "(c)") +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9),
        legend.position = c(0.98, 0.02),
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = alpha("white", 0.8),
                                          colour = NA),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 5.5))
save_panel(p, "23_figS1c_conf_dg_class", width = 5.5, height = 4)
