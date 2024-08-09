path <- getwd()
setwd(path)

temp_lib_path <- tempfile(pattern="Rlib")
dir.create(temp_lib_path)

.libPaths(c(temp_lib_path, .libPaths()))

message("### Install depedencies ###")

devtools::install_local(path = path,
  force = TRUE,
  build = TRUE,
  upgrade = "never",
  build_vignettes = FALSE)

devtools::build_vignettes()
dir.create("inst/doc")
file.copy(dir("doc", full.names=TRUE), "inst/doc", overwrite=TRUE)
