# SI Fig S3 panel (c) — median Cys-SG distance by substrate class for each isoform.
suppressPackageStartupMessages({
  library(ggplot2)
})
source("_palettes.R")

d <- read.csv(file.path(DATA2PLOT, "figS3c_med_dist.csv"),
              fileEncoding = "UTF-8", check.names = FALSE)
d$protein  <- factor(d$protein, levels = PROTEINS)
d$category <- factor(d$category, levels = CAT_ORDER)

p <- ggplot(d, aes(x = category, y = median_distance,
                   colour = protein, group = protein)) +
  geom_line(linewidth = 0.5, alpha = 0.85) +
  geom_point(size = 1.8, alpha = 0.9) +
  geom_hline(yintercept = 4.0, linetype = "dashed", colour = "#D62728",
             linewidth = 0.4) +
  scale_colour_manual(values = PROT_COLORS, name = NULL) +
  labs(x = NULL, y = "Median Cys-SG distance (Å)",
       title = "Catalytic geometry by substrate class\nand isoform",
       subtitle = "(c)") +
  guides(colour = guide_legend(ncol = 2)) +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9),
        axis.text.x   = element_text(angle = 30, hjust = 1, size = 7),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 0.8),
                                          colour = NA),
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 6.5))
save_panel(p, "34_figS3c_median_dist", width = 6, height = 4)
