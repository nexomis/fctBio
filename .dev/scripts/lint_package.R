#!/usr/bin/env Rscript

tmp_lib <- ".tmp/lib-lint"
dir.create(tmp_lib, recursive = TRUE, showWarnings = FALSE)
withr::with_libpaths(tmp_lib, {
  devtools::install_local(
    build_vignettes = FALSE,
    force = FALSE,
    upgrade = "never"
  )
  library(fctBio)
  lints <- lintr::lint_package()
  lintr::sarif_output(lints, filename = "lintr_results.sarif")
})
