# Fig 2 panel (f) — productive fraction per substrate class.
suppressPackageStartupMessages({
  library(ggplot2)
})
source("_palettes.R")

d  <- read.csv(file.path(DATA2PLOT, "fig2f_class_prod.csv"),
               fileEncoding = "UTF-8", check.names = FALSE)
ov <- read.csv(file.path(DATA2PLOT, "fig2f_overall.csv"))
# drop "Other" to match Python (cats_sorted = CAT_ORDER[:-1])
d <- d[d$category != "Other", ]
# sort by descending frac as in Python
d <- d[order(-d$frac), ]
d$category <- factor(d$category, levels = d$category)

p <- ggplot(d, aes(y = category, x = frac, fill = category)) +
  geom_col(height = 0.6, width = 0.6) +
  geom_vline(xintercept = ov$overall_frac[1], linetype = "dashed",
             colour = "black", alpha = 0.6, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.0f%%", frac * 100)),
            hjust = -0.1, size = 2.5) +
  scale_fill_manual(values = CAT_COLORS, guide = "none") +
  scale_x_continuous(limits = c(0, 1.05), labels = scales::percent) +
  labs(x = "Productive fraction (≤ 4 Å)", y = NULL,
       title = "Productive poses per\nsubstrate class",
       subtitle = "(f)") +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9))
save_panel(p, "16_fig2f_class_productive", width = 5.5, height = 4)
