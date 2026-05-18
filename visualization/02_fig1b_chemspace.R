# Fig 1 panel (b) — chemical space scatter (MW vs logP, coloured by class).
suppressPackageStartupMessages({
  library(ggplot2)
})
source("_palettes.R")

d <- read.csv(file.path(DATA2PLOT, "fig1b_chemspace.csv"),
              check.names = FALSE, fileEncoding = "UTF-8")
d$category <- factor(d$category, levels = CAT_ORDER)

p <- ggplot(d, aes(x = logP, y = MW, fill = category)) +
  geom_point(shape = 21, size = 2.4, alpha = 0.85, colour = "white",
             stroke = 0.3) +
  scale_fill_manual(values = CAT_COLORS, drop = FALSE, name = NULL) +
  labs(x = "Calculated log P", y = "Exact MW (Da)",
       title = "Chemical space of 60 substrates",
       subtitle = "(b)") +
  theme_classic(base_size = 9) +
  theme(plot.title    = element_text(size = 9, hjust = 0.5),
        plot.subtitle = element_text(face = "plain", hjust = -0.05),
        legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7),
                                          colour = NA),
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 6))
save_panel(p, "02_fig1b_chemspace", width = 5.5, height = 4)
