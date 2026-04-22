conda_bin <- "/home/trhova/miniconda3/bin"
conda_lib <- "/home/trhova/miniconda3/lib"
conda_pkgconfig <- file.path(conda_lib, "pkgconfig")

if (dir.exists(conda_bin)) {
  Sys.setenv(PATH = paste(conda_bin, Sys.getenv("PATH"), sep = .Platform$path.sep))
}

if (dir.exists(conda_lib)) {
  Sys.setenv(
    LD_LIBRARY_PATH = paste(conda_lib, Sys.getenv("LD_LIBRARY_PATH"), sep = .Platform$path.sep)
  )
}

if (dir.exists(conda_pkgconfig)) {
  Sys.setenv(
    PKG_CONFIG_PATH = paste(conda_pkgconfig, Sys.getenv("PKG_CONFIG_PATH"), sep = .Platform$path.sep)
  )
}

external_libs <- unique(c(Sys.getenv("R_LIBS_USER"), .Library.site, .Library))
external_libs <- external_libs[nzchar(external_libs) & dir.exists(external_libs)]

if (length(external_libs) > 0) {
  Sys.setenv(
    RENV_CONFIG_EXTERNAL_LIBRARIES = paste(external_libs, collapse = .Platform$path.sep)
  )
}

if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}
