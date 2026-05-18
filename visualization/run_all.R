# Run every numbered panel script from the visualization/ directory.
# Must be run with cwd = ShenYuchen/visualization (or sourced from R there).
scripts <- sort(list.files(".", pattern = "^[0-9]{2}_.*\\.R$"))
for (s in scripts) {
  cat("\n──", s, "──\n")
  source(s, local = new.env())
}
cat("\nAll panels rendered to ../figures/panels_R/\n")
