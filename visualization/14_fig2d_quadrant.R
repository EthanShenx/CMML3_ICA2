# Fig 2 panel (d) — two-metric quadrant scatter (distance × ΔG thresholds).
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("_palettes.R")

d <- read.csv(file.path(DATA2PLOT, "fig2d_quadrant.csv"))
d <- d[!is.na(d$distance), ]

QUAD_COLS <- c(
  "Productive high-affinity"     = "#FFD700",
  "High-affinity non-productive" = "#E07B39",
  "Productive low-affinity"      = "#2196F3",
  "Neither"                      = "#BDBDBD"
)
counts <- d %>%
  group_by(quad) %>% summarise(n = n()) %>%
  mutate(label = sprintf("%s (n=%d)", quad, n))
d$quad <- factor(d$quad, levels = names(QUAD_COLS))
counts$quad <- factor(counts$quad, levels = names(QUAD_COLS))

lab_map <- setNames(counts$label[match(names(QUAD_COLS), counts$quad)],
                    names(QUAD_COLS))

p <- ggplot(d, aes(x = distance, y = DG, colour = quad)) +
  geom_vline(xintercept = 4.0, linetype = "dashed", colour = "grey50",
             linewidth = 0.4) +
  geom_hline(yintercept = -7.0, linetype = "dashed", colour = "grey50",
             linewidth = 0.4) +
  geom_point(size = 1.2, alpha = 0.72, stroke = 0) +
  scale_colour_manual(values = QUAD_COLS, labels = lab_map, name = NULL,
                      drop = FALSE) +
  labs(x = "Distance (Å)",
       y = expression(Delta * G~"(kcal mol"^{-1}*")"),
       title = "Two-metric quadrant\nanalysis",
       subtitle = "(d)") +
  guides(colour = guide_legend(override.aes = list(size = 2.5))) +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9),
        legend.position = c(0.98, 0.02),
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = alpha("white", 0.8),
                                          colour = NA),
        legend.text = element_text(size = 5.5),
        legend.key.size = unit(0.3, "cm"))
save_panel(p, "14_fig2d_quadrant", width = 5, height = 4)
