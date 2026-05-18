# SI Fig S1 panel (d) — ΔG distribution by substrate class with Kruskal–Wallis.
suppressPackageStartupMessages({
  library(ggplot2)
})
source("_palettes.R")

d  <- read.csv(file.path(DATA2PLOT, "figS1d_dg_class.csv"),
               fileEncoding = "UTF-8", check.names = FALSE)
kw <- read.csv(file.path(DATA2PLOT, "figS1d_kruskal.csv"))
d$category <- factor(d$category, levels = CAT_ORDER)

p <- ggplot(d, aes(x = category, y = DG, fill = category)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.3) +
  scale_fill_manual(values = CAT_COLORS, guide = "none", drop = FALSE) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Kruskal–Wallis\n%s", format_p(kw$p[1])),
           hjust = 1.05, vjust = 1.1, size = 2.4, fill = "white") +
  labs(x = NULL, y = expression(Delta * G~"(kcal mol"^{-1}*")"),
       title = "Binding affinity by substrate class",
       subtitle = "(d)") +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9),
        axis.text.x   = element_text(angle = 30, hjust = 1, size = 7))
save_panel(p, "24_figS1d_dg_boxplot", width = 5.5, height = 4)
