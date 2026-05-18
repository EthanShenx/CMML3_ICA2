# Fig 1 panel (d) — ALDH × substrate ΔG heatmap with Ward dendrogram on top.
# Dendrogram segments are exported from Python (scipy) and replayed here, so
# the tree visualised matches the clustering that ordered the heatmap columns.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
source("_palettes.R")

long <- read.csv(file.path(DATA2PLOT, "fig1d_heatmap_long.csv"),
                 check.names = FALSE, fileEncoding = "UTF-8")
top1 <- read.csv(file.path(DATA2PLOT, "fig1d_top1.csv"),
                 check.names = FALSE, fileEncoding = "UTF-8")
seg  <- read.csv(file.path(DATA2PLOT, "fig1d_dendro_segments.csv"))

long$protein   <- factor(long$protein, levels = rev(PROTEINS))
long$col_index <- as.integer(long$col_index)
long$lig_num   <- factor(long$lig_num,
                         levels = unique(long$lig_num[order(long$col_index)]))

xlabels <- long %>%
  group_by(lig_num) %>%
  summarise(category = category[1]) %>%
  arrange(match(lig_num, levels(long$lig_num)))
xlabels$col <- CAT_COLORS[xlabels$category]

top1$protein   <- factor(top1$protein, levels = rev(PROTEINS))
top1$col_index <- as.integer(top1$col_index)

n_cols <- length(levels(long$lig_num))
x_lim  <- c(-0.5, n_cols - 0.5)

# ── Dendrogram subplot (top) ──
p_dendro <- ggplot(seg) +
  geom_segment(aes(x = x, xend = xend, y = y, yend = yend),
               colour = "black", linewidth = 0.3) +
  scale_x_continuous(limits = x_lim, expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 0, 2))

# ── Heatmap subplot (bottom) ──
p_heat <- ggplot(long, aes(x = col_index, y = protein, fill = DG)) +
  geom_tile(colour = NA) +
  geom_text(data = top1, aes(x = col_index, y = protein),
            label = "★", colour = "white", size = 3, fontface = "plain",
            inherit.aes = FALSE) +
  scale_fill_distiller(palette = "RdYlBu", direction = -1, na.value = "grey90",
                       limits = c(-10.5, -5.0),
                       name = expression(Delta * G[noelec]~"(kcal mol"^{-1}*")")) +
  scale_x_continuous(breaks = seq_len(n_cols) - 1L,
                     labels = levels(long$lig_num),
                     limits = x_lim, expand = c(0, 0)) +
  labs(x = "Substrate (number, coloured by class; columns hierarchically clustered)",
       y = "ALDH isoform",
       subtitle = "(d)") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                    size = 5.5,
                                    colour = xlabels$col),
        axis.text.y = element_text(size = 8),
        plot.subtitle = element_text(face = "plain", hjust = -0.02),
        legend.position = "right",
        legend.key.width  = unit(0.4, "cm"),
        legend.key.height = unit(1.2, "cm"),
        plot.margin = margin(0, 2, 2, 2))

# ── Stack: dendrogram on top, heatmap on bottom ──
p <- (p_dendro / p_heat) +
  plot_layout(heights = c(1.2, 4))

save_panel(p, "04_fig1d_affinity_heatmap", width = 14, height = 5)
