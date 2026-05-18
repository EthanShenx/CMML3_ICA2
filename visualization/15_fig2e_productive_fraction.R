# Fig 2 panel (e) — productive fraction per isoform with Wilson 95% CI bars.
suppressPackageStartupMessages({
  library(ggplot2)
})
source("_palettes.R")

d   <- read.csv(file.path(DATA2PLOT, "fig2e_productive.csv"))
ov  <- read.csv(file.path(DATA2PLOT, "fig2e_overall.csv"))
d$protein <- factor(d$protein, levels = PROTEINS)

p <- ggplot(d, aes(x = protein, y = frac, fill = protein)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = frac - (frac - ci_lo),
                    ymax = frac + (ci_hi - frac)),
                width = 0.18, linewidth = 0.6) +
  geom_hline(yintercept = ov$overall_frac[1], linetype = "dashed",
             colour = "black", alpha = 0.6, linewidth = 0.4) +
  scale_fill_manual(values = PROT_COLORS, guide = "none") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  annotate("text", x = 0.6, y = ov$overall_frac[1] + 0.04,
           label = sprintf("Overall (%.0f%%)", ov$overall_frac[1] * 100),
           hjust = 0, size = 2.4) +
  labs(x = NULL, y = "Productive fraction",
       title = "Productive poses\nper isoform (95% CI)",
       subtitle = "(e)") +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9))
save_panel(p, "15_fig2e_productive_fraction", width = 5, height = 4)
