#!/usr/bin/env Rscript

tmp_lib <- tempdir(check = TRUE)
dir.create(tmp_lib, recursive = TRUE, showWarnings = FALSE)
withr::with_libpaths(tmp_lib, {
  devtools::install_local(
    quick = TRUE,
    force = TRUE,
    upgrade_dependencies = FALSE
  )
  library(fctBio)
  capture.output(
    lintr::lint_package(),
    file = "lintr_report.txt"
  )
})
