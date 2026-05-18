# Fig 1 panel (f) — top-5 substrates per protein, lollipop chart.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("_palettes.R")

d <- read.csv(file.path(DATA2PLOT, "fig1f_top5.csv"),
              fileEncoding = "UTF-8", check.names = FALSE)

d$protein <- factor(d$protein, levels = PROTEINS)
d <- d %>% arrange(protein, DG_prediction_kcalmol)
# Build sequential y positions with a gap between proteins (matches Python)
d$y <- NA_real_
y <- 0
for (i in seq_len(nrow(d))) {
  d$y[i] <- y
  y <- y + 1
  if (i < nrow(d) && d$protein[i + 1] != d$protein[i]) y <- y + 0.5
}

p <- ggplot(d, aes(x = DG_prediction_kcalmol, y = y, colour = protein)) +
  geom_segment(aes(xend = -5.1, yend = y), linewidth = 0.6, alpha = 0.7) +
  geom_point(size = 2.2) +
  geom_vline(xintercept = -7, linetype = "dashed", colour = "grey50",
             linewidth = 0.3) +
  scale_y_continuous(breaks = d$y, labels = d$short_name) +
  scale_colour_manual(values = PROT_COLORS, name = NULL) +
  labs(x = expression(Delta * G~"(kcal mol"^{-1}*")"), y = NULL,
       title = "Top-5 substrates per isoform",
       subtitle = "(f)") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3))) +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 5.5),
        plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom",
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 6))
save_panel(p, "06_fig1f_top5_lollipop", width = 6.5, height = 6)
