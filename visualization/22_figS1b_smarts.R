# SI Fig S1 panel (b) — SMARTS aldehyde identification per isoform.
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
})
source("_palettes.R")

d <- read.csv(file.path(DATA2PLOT, "figS1b_smarts.csv"))
d$protein <- factor(d$protein, levels = PROTEINS)
long <- pivot_longer(d, c(n_aldehyde, n_nonaldehyde),
                     names_to = "type", values_to = "n")
long$type <- factor(long$type,
                     levels = c("n_aldehyde", "n_nonaldehyde"),
                     labels = c("Aldehyde (SMARTS ✓)", "Non-aldehyde"))

p <- ggplot(long, aes(x = protein, y = n, fill = type)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("Aldehyde (SMARTS ✓)" = "#5C85D6",
                                "Non-aldehyde"        = "#FF8C00"),
                    name = NULL) +
  labs(x = NULL, y = "Number of substrates",
       title = "SMARTS aldehyde identification\nper isoform",
       subtitle = "(b)") +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9),
        legend.position = "top",
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7))
save_panel(p, "22_figS1b_smarts", width = 6, height = 4)
