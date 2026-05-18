# Fig 2 panel (a) — Cys-SG ↔ aldehyde-C distance histogram with 4 Å threshold.
suppressPackageStartupMessages({
  library(ggplot2)
})
source("_palettes.R")

d <- read.csv(file.path(DATA2PLOT, "fig2a_distance.csv"))
n_prod <- sum(d$distance <= 4.0)
n_tot  <- nrow(d)

p <- ggplot(d, aes(x = distance)) +
  annotate("rect", xmin = 0, xmax = 4, ymin = 0, ymax = Inf,
           fill = "#A8D5A2", alpha = 0.15) +
  geom_histogram(bins = 40, fill = "#5C85D6", colour = "white",
                 linewidth = 0.2, alpha = 0.85) +
  geom_vline(xintercept = 4.0, colour = "#D62728", linetype = "dashed",
             linewidth = 0.7) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("%d/%d\nproductive\n(≤ 4 Å)", n_prod, n_tot),
           hjust = 1.1, vjust = 1.1, size = 2.5,
           colour = "#D62728", fill = "white",
           label.size = 0.3) +
  labs(x = "Cys-SG ↔ aldehyde-C distance (Å)",
       y = "Number of complexes",
       title = "Distance distribution\n(all aldehyde pairs)",
       subtitle = "(a)") +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9))
save_panel(p, "09_fig2a_distance_hist", width = 5, height = 4)
