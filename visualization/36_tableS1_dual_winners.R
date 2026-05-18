# Supplementary Table S1 — top-20 dual-winner pairs.
# Renders a graphical table image AND emits report/table_S1_dual_winners.tex
# matching the LaTeX produced by the original Python script.
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})
source("_palettes.R")

dw <- read.csv(file.path(DATA2PLOT, "tableS1_dual_winners.csv"),
                fileEncoding = "UTF-8", check.names = FALSE)
dw$category <- factor(dw$category, levels = CAT_ORDER)

# ── Graphical table (PNG image)
tbl_df <- data.frame(
  Protein   = dw$protein,
  Substrate = dw$substrate,
  "ΔG"      = sprintf("%.2f", dw$DG),
  "Dist."   = sprintf("%.2f", dw$distance),
  Confidence = ifelse(is.na(dw$confidence), "—",
                       sprintf("%.2f", dw$confidence)),
  Class     = as.character(dw$category),
  check.names = FALSE
)
tt <- ttheme_minimal(
  core    = list(fg_params = list(fontsize = 8,
                                   fontface = ifelse(dw$near, "bold", "plain"))),
  colhead = list(fg_params = list(fontsize = 8.5, fontface = "plain"))
)
grb <- tableGrob(tbl_df, rows = NULL, theme = tt)
png(file.path(FIG_OUT, "36_tableS1_dual_winners.png"),
    width = 9, height = 7, units = "in", res = 300, bg = "white")
grid.draw(grb)
invisible(dev.off())
message("✓ 36_tableS1_dual_winners.png")

# ── LaTeX (identical scheme to 07_make_figures.py make_table_S1)
tex_safe <- function(s) {
  s <- as.character(s)
  s <- gsub("α", "$\\alpha$", s, fixed = TRUE)
  s <- gsub("β", "$\\beta$",  s, fixed = TRUE)
  s <- gsub("→", "$\\rightarrow$", s, fixed = TRUE)
  s <- gsub("↔", "$\\leftrightarrow$", s, fixed = TRUE)
  s <- gsub("★", "$\\star$", s, fixed = TRUE)
  s <- gsub("≤", "$\\leq$",  s, fixed = TRUE)
  s <- gsub("≥", "$\\geq$",  s, fixed = TRUE)
  s
}
report_dir <- normalizePath(file.path(dirname(FIG_OUT), "..", "report"),
                             mustWork = FALSE)
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

lines <- c(
  "\\begin{table*}[t]",
  "\\centering",
  "\\small",
  "\\caption{\\textbf{Table S1 --- Top-20 dual-winner protein--substrate pairs.} ",
  "Dual-winner: PRODIGY-LIG $\\Delta G_{\\text{noelec}} \\leq -7.0$\\,kcal\\,mol$^{-1}$ ",
  "\\textit{and} catalytic Cys-SG $\\leq 4$\\,\\AA. ",
  "\\textbf{Bold} rows satisfy the near-attack criterion ($\\leq 2.5$\\,\\AA). ",
  "Substrate class colour coding matches Figs.~1--2.}",
  "\\label{tab:S1}",
  "\\begin{tabular}{@{}llrrrl@{}}",
  "\\toprule",
  "\\textbf{Protein} & \\textbf{Substrate} & $\\boldsymbol{\\Delta G}$ & \\textbf{Dist.} & \\textbf{Confidence} & \\textbf{Class} \\\\",
  " & & \\textbf{(kcal mol$^{-1}$)} & \\textbf{(\\AA)} & & \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(dw))) {
  row <- dw[i, ]
  conf_str <- if (is.na(row$confidence)) "---" else sprintf("%.2f", row$confidence)
  cells <- c(tex_safe(row$protein), tex_safe(row$substrate),
             sprintf("%.2f", row$DG), sprintf("%.2f", row$distance),
             conf_str, tex_safe(as.character(row$category)))
  if (isTRUE(row$near)) cells <- paste0("\\textbf{", cells, "}")
  lines <- c(lines, paste0("  ", paste(cells, collapse = " & "), " \\\\"))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table*}")
out_tex <- file.path(report_dir, "table_S1_dual_winners.tex")
writeLines(lines, out_tex)
message("✓ ", out_tex)
