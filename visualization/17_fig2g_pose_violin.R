# Fig 2 panel (g) — ΔG violin: productive vs non-productive poses with MWU p.
suppressPackageStartupMessages({
  library(ggplot2)
})
source("_palettes.R")

d   <- read.csv(file.path(DATA2PLOT, "fig2g_pose_violin.csv"),
                fileEncoding = "UTF-8", check.names = FALSE)
mwu <- read.csv(file.path(DATA2PLOT, "fig2g_mwu.csv"))

POSE_COLS <- c("Productive (≤ 4 Å)"      = "#5C85D6",
               "Non-productive (> 4 Å)"  = "#D62728")
d$pose <- factor(d$pose, levels = names(POSE_COLS))

p <- ggplot(d, aes(x = pose, y = DG, fill = pose)) +
  geom_violin(scale = "width", trim = TRUE, colour = "grey40",
              linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.4, fill = "white",
               colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = POSE_COLS, guide = "none") +
  annotate("label", x = 1.5, y = min(d$DG) + 0.2,
           label = sprintf("Mann–Whitney\n%s", format_p(mwu$p[1])),
           hjust = 0.5, vjust = 0, size = 2.3,
           label.size = 0.3, fill = "white") +
  labs(x = NULL, y = expression(Delta * G~"(kcal mol"^{-1}*")"),
       title = "ΔG: productive vs\nnon-productive poses",
       subtitle = "(g)") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(size = 7),
        plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9))
save_panel(p, "17_fig2g_pose_violin", width = 5, height = 4)
