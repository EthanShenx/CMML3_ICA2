# Fig 1 panel (g) — isoform × substrate binding dot matrix (vertical layout).
# Columns = substrates with ≥ 1 isoform achieving ΔG ≤ −6, sorted by
# promiscuity (n_prots) descending then by substrate ID. Rows = the seven
# isoforms. A filled coloured dot marks each (protein, ligand) pair that
# meets the threshold; non-binding cells are rendered as small grey rings.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("_palettes.R")

m <- read.csv(file.path(DATA2PLOT, "fig1g_matrix.csv"),
              fileEncoding = "UTF-8", check.names = FALSE)
m$protein  <- factor(m$protein, levels = rev(PROTEINS))   # A8 on top
m$category <- factor(m$category, levels = CAT_ORDER)

# Keep substrates that bind at least one isoform; sort by promiscuity desc
keep <- m %>%
  group_by(ligand) %>%
  summarise(n_prots = max(n_prots), lig_num = lig_num[1], .groups = "drop") %>%
  filter(n_prots >= 1) %>%
  arrange(desc(n_prots), lig_num)

m <- m[m$ligand %in% keep$ligand, ]

# Build x-axis labels: short_name only (the leading ID is redundant)
lab <- m %>%
  group_by(ligand, lig_num, short_name) %>%
  summarise(.groups = "drop") %>%
  mutate(label = short_name)
keep$label <- lab$label[match(keep$ligand, lab$ligand)]
m$ligand   <- factor(m$ligand, levels = keep$ligand)

strong <- m[m$is_strong == 1, ]
weak   <- m[m$is_strong == 0, ]

p <- ggplot(m, aes(x = ligand, y = protein)) +
  geom_point(data = weak,   shape = 21, fill = NA, colour = "grey80",
             size = 1.6, stroke = 0.4) +
  geom_point(data = strong, aes(fill = category),
             shape = 21, colour = "black", size = 3, stroke = 0.3) +
  scale_fill_manual(values = CAT_COLORS, name = NULL, drop = FALSE) +
  scale_x_discrete(breaks = keep$ligand, labels = keep$label) +
  labs(x = NULL, y = NULL,
       title = "Pan-ALDH binding at ΔG ≤ −6 kcal mol⁻¹",
       subtitle = "(g)") +
  guides(fill = guide_legend(override.aes = list(size = 3),
                              nrow = 1)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.x = element_line(colour = "grey92",
                                            linewidth = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.subtitle = element_text(face = "plain", hjust = -0.02),
        plot.title    = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.text   = element_text(size = 7))

# Adapt width to column count so substrate labels stay readable
w <- max(7, 0.28 * nrow(keep) + 2)
save_panel(p, "07_fig1g_pan_aldh_matrix", width = w, height = 5.5)
