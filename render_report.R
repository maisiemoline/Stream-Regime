# ==============================================================================
# Script: render_report.R
# Description: Helper script to knit and render stream_regime_classification.Rmd
# ==============================================================================

# Ensure Pandoc from RStudio or system is detected
pandoc_dirs <- c(
  "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64",
  "/Applications/RStudio.app/Contents/Resources/app/resources/pandoc",
  "/usr/local/bin",
  "/opt/homebrew/bin"
)

for (p in pandoc_dirs) {
  if (file.exists(file.path(p, "pandoc"))) {
    Sys.setenv(RSTUDIO_PANDOC = p)
    rmarkdown::find_pandoc(dir = p)
    break
  }
}

message("Rendering stream_regime_classification.Rmd to HTML...")
rmarkdown::render(
  input = "stream_regime_classification.Rmd",
  output_file = "stream_regime_classification.html",
  quiet = FALSE
)

message("Report successfully rendered: stream_regime_classification.html")
