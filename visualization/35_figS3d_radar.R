# SI Fig S3 panel (d) — radar chart of isoform selectivity per substrate class.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("_palettes.R")

d <- read.csv(file.path(DATA2PLOT, "figS3d_radar.csv"),
              fileEncoding = "UTF-8", check.names = FALSE)
# Match Python: cats_r = CAT_ORDER[:-1]  (six classes, drop "Other")
cats_r <- CAT_ORDER[CAT_ORDER != "Other"]
d <- d[d$category %in% cats_r, ]
d$category <- factor(d$category, levels = cats_r)
d$protein  <- factor(d$protein,  levels = PROTEINS)
d$angle    <- as.numeric(d$category)

# Close each polygon by repeating the first category at the end
close_df <- d %>%
  group_by(protein) %>%
  filter(angle == min(angle)) %>%
  mutate(angle = length(cats_r) + 1)
plot_df <- bind_rows(d, close_df)

p <- ggplot(plot_df, aes(x = angle, y = frac,
                         colour = protein, group = protein)) +
  geom_polygon(aes(fill = protein), alpha = 0.08, linewidth = 0.6) +
  geom_path(linewidth = 0.6, alpha = 0.85) +
  coord_polar() +
  scale_x_continuous(breaks = seq_along(cats_r), labels = cats_r,
                     limits = c(0.5, length(cats_r) + 0.5)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0.25, 0.5, 0.75, 1.0),
                     labels = c("25%", "50%", "75%", "100%")) +
  scale_colour_manual(values = PROT_COLORS, name = NULL) +
  scale_fill_manual(values = PROT_COLORS, guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Selectivity profile\n(fraction with ΔG ≤ −7 kcal mol⁻¹)",
       subtitle = "(d)") +
  theme_minimal(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9),
        axis.text.x   = element_text(size = 7),
        axis.text.y   = element_text(size = 6),
        legend.key.size = unit(0.35, "cm"),
        legend.text   = element_text(size = 7))
save_panel(p, "35_figS3d_radar", width = 6, height = 5)
