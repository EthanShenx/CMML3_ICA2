# Fig 2 panel (c) — ΔG vs Cys-SG distance scatter, coloured by DiffDock confidence.
suppressPackageStartupMessages({
  library(ggplot2)
})
source("_palettes.R")

d   <- read.csv(file.path(DATA2PLOT, "fig2c_scatter.csv"))
sp  <- read.csv(file.path(DATA2PLOT, "fig2c_spearman.csv"))

p <- ggplot(d, aes(x = distance, y = DG, colour = confidence)) +
  geom_point(size = 1.3, alpha = 0.75) +
  geom_vline(xintercept = 4.0, linetype = "dashed", colour = "grey50",
             linewidth = 0.4) +
  scale_colour_viridis_c(name = "DiffDock confidence") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Spearman ρ = %.2f\n%s",
                            sp$rho[1], format_p(sp$p[1])),
           hjust = 1.05, vjust = 1.05, size = 2.4,
           label.size = 0.3, fill = "white") +
  labs(x = "Cys-SG distance (Å)",
       y = expression(Delta * G~"(kcal mol"^{-1}*")"),
       title = "Affinity vs catalytic geometry",
       subtitle = "(c)") +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 7))
save_panel(p, "11_fig2c_dg_vs_dist", width = 5, height = 4)
