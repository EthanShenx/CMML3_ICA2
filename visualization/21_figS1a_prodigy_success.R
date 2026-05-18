# SI Fig S1 panel (a) — PRODIGY-LIG pipeline success rates per isoform.
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
})
source("_palettes.R")

d <- read.csv(file.path(DATA2PLOT, "figS1a_prodigy.csv"))
d$protein <- factor(d$protein, levels = PROTEINS)
long <- pivot_longer(d, c(success_pct, fail_pct),
                     names_to = "outcome", values_to = "pct")
long$outcome <- factor(long$outcome,
                        levels = c("success_pct", "fail_pct"),
                        labels = c("Success (rc=0)", "Failure (rc≠0)"))

p <- ggplot(long, aes(x = protein, y = pct, fill = outcome)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("Success (rc=0)" = "#4CAF50",
                                "Failure (rc≠0)" = "#F44336"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 110)) +
  labs(x = NULL, y = "PRODIGY-LIG run outcome (%)",
       title = "PRODIGY-LIG pipeline success rates",
       subtitle = "(a)") +
  theme_classic(base_size = 9) +
  theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
        plot.title    = element_text(hjust = 0.5, size = 9),
        legend.position = c(0.98, 0.05),
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = alpha("white", 0.7),
                                          colour = NA),
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7))
save_panel(p, "21_figS1a_prodigy_success", width = 6, height = 4)
