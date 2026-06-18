project_root <- normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
if (basename(getwd()) != "paper") {
  project_root <- normalizePath(getwd(), mustWork = TRUE)
}

local_library <- file.path(project_root, ".Rlib")
.libPaths(c(local_library, .libPaths()))

if (!requireNamespace("rjtools", quietly = TRUE)) {
  stop(
    "Install rjtools first with:\n",
    "install.packages('rjtools')",
    call. = FALSE
  )
}

if (!rmarkdown::pandoc_available()) {
  rstudio_pandoc <- file.path(
    "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools",
    if (grepl("arm64|aarch64", R.version$platform)) "aarch64" else "x86_64"
  )
  if (file.exists(file.path(rstudio_pandoc, "pandoc"))) {
    Sys.setenv(RSTUDIO_PANDOC = rstudio_pandoc)
  }
}

article <- file.path(project_root, "paper", "tensormodels.Rmd")

rmarkdown::render(
  article,
  output_format = "rjtools::rjournal_web_article",
  clean = TRUE
)

rmarkdown::render(
  article,
  output_format = "rjtools::rjournal_pdf_article",
  clean = TRUE
)
