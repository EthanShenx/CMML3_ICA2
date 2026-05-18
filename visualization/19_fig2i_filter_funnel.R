# Fig 2 panel (i) — three-gate filter funnel: affinity × geometry × LE.
suppressPackageStartupMessages({
  library(ggplot2)
})
source("_palettes.R")

d <- read.csv(file.path(DATA2PLOT, "fig2i_funnel.csv"),
              fileEncoding = "UTF-8", check.names = FALSE)
# Preserve the stage order (0..3 from Python)
d <- d[order(d$order), ]
d$stage <- factor(d$stage, levels = d$stage)

# Compute retention labels
prev <- c(NA, head(d$n, -1))
retain <- ifelse(is.na(prev), "",
                  sprintf("  (%d%% retained)", round(100 * d$n / prev)))
d$annot <- sprintf("n = %d%s", d$n, retain)

cols_funnel <- c("#9E9E9E", "#5C85D6", "#FFB300", "#2E7D32")

p <- ggplot(d, aes(y = stage, x = n)) +
  geom_col(fill = cols_funnel, colour = "black", linewidth = 0.4,
           width = 0.6) +
  geom_text(aes(label = annot), hjust = -0.05, size = 2.6) +
  scale_y_discrete(limits = rev(levels(d$stage))) +
  scale_x_continuous(limits = c(0, max(d$n) * 1.35)) +
  labs(x = "Pairs retained", y = NULL,
       title = "Three-gate filter funnel: affinity × geometry × ligand efficiency",
       subtitle = "(i)") +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.02),
        plot.title    = element_text(size = 9),
        axis.text.y   = element_text(size = 8))
save_panel(p, "19_fig2i_filter_funnel", width = 7, height = 4)
