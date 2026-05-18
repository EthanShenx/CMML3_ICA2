# SI Fig S3 panel (a) — Tanimoto similarity dendrogram (Morgan FP, r=2).
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggdendro)
})
source("_palettes.R")

dist_mat <- as.matrix(read.csv(file.path(DATA2PLOT, "figS3a_distmat.csv"),
                                row.names = 1, check.names = FALSE))
meta <- read.csv(file.path(DATA2PLOT, "figS3a_meta.csv"),
                  fileEncoding = "UTF-8", check.names = FALSE)
d <- as.dist(dist_mat)
hc <- hclust(d, method = "ward.D2")
dendro_data <- dendro_data(hc, type = "rectangle")

lbl <- label(dendro_data)
lbl$category <- meta$category[match(as.character(lbl$label), meta$label)]

p <- ggplot() +
  geom_segment(data = segment(dendro_data),
               aes(x = x, y = y, xend = xend, yend = yend),
               colour = "grey50", linewidth = 0.3) +
  geom_text(data = lbl,
            aes(x = x, y = -0.02, label = label, colour = category),
            angle = 90, hjust = 1, vjust = 0.5, size = 1.7) +
  scale_colour_manual(values = CAT_COLORS, drop = FALSE, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.18, 0.05))) +
  labs(x = "Substrate number (colour = class)",
       y = "Ward linkage distance",
       title = "Tanimoto similarity dendrogram\n(Morgan FP, r=2)",
       subtitle = "(a)") +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9),
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank())
save_panel(p, "32_figS3a_dendrogram", width = 7, height = 4.5)
