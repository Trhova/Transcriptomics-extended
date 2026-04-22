Sys.setenv(
  PATH = paste("/usr/bin", Sys.getenv("PATH"), sep = .Platform$path.sep),
  PKG_CONFIG = "/usr/bin/pkg-config"
)
Sys.unsetenv(c("PKG_CONFIG_LIBDIR", "ICU_CFLAGS", "ICU_LIBS"))

conda_lib <- "/home/trhova/miniconda3/lib"
conda_pkgconfig <- file.path(conda_lib, "pkgconfig")
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

dirs <- c(
  "data",
  "data/signatures",
  "results",
  "results/figures",
  "results/tables",
  "results/objects",
  "scripts"
)

invisible(lapply(dirs, function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}))

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

manifest <- utils::read.delim(
  "data/package_manifest.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

if (!file.exists("renv/activate.R")) {
  renv::init(bare = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

installed <- rownames(installed.packages())
missing_cran <- manifest$package[manifest$source == "cran" & !manifest$package %in% installed]
missing_bioc <- manifest$package[manifest$source == "bioconductor" & !manifest$package %in% installed]

if (length(missing_cran) > 0) {
  install.packages(
    missing_cran,
    repos = "https://cloud.r-project.org",
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

if (length(missing_bioc) > 0) {
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE,
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

patchelf <- Sys.which("patchelf")
magick_lib <- tryCatch(find.package("magick"), error = function(...) "")
magick_so <- file.path(magick_lib, "libs", "magick.so")

if (nzchar(patchelf) && nzchar(magick_lib) && file.exists(magick_so) && dir.exists(conda_lib)) {
  system2(patchelf, c("--set-rpath", conda_lib, magick_so))
}

utils::write.table(
  manifest,
  file = "data/package_manifest.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
renv::snapshot(prompt = FALSE, force = TRUE)

message("Setup complete. Packages recorded in renv.lock and data/package_manifest.tsv.")
