suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(forcats)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

dir_create_if_missing <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

initialize_repo_dirs <- function() {
  dirs <- c(
    "data",
    "data/signatures",
    "results",
    "results/figures",
    "results/tables",
    "results/objects",
    "scripts"
  )
  invisible(lapply(dirs, dir_create_if_missing))
}

read_package_manifest <- function(path = "data/package_manifest.tsv") {
  readr::read_tsv(path, show_col_types = FALSE)
}

check_required_packages <- function(extra = character()) {
  manifest <- read_package_manifest()
  required <- unique(c(manifest$package, extra))
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing, collapse = ", "),
      ". Run `Rscript scripts/00_setup.R` first."
    )
  }
}

load_packages <- function(packages) {
  invisible(lapply(packages, function(pkg) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }))
}

analysis_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, color = "#13293D"),
      plot.subtitle = element_text(color = "#406E8E"),
      axis.title = element_text(face = "bold", color = "#13293D"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
}

save_plot <- function(plot, filename, width = 8, height = 5, dpi = 300) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

write_csv_table <- function(data, path) {
  readr::write_csv(data, path)
  invisible(data)
}

read_airway_dataset <- function() {
  load_packages(c("airway", "SummarizedExperiment"))
  data("airway", package = "airway")
  airway
}

build_sample_metadata <- function(se) {
  meta <- as.data.frame(SummarizedExperiment::colData(se)) %>%
    tibble::rownames_to_column("sample_id")

  if ("Run" %in% colnames(meta) && !"run" %in% colnames(meta)) {
    meta <- dplyr::rename(meta, run = Run)
  }

  if (!"run" %in% colnames(meta)) {
    meta$run <- NA_character_
  }

  meta <- meta %>%
    dplyr::mutate(
      sample_id = colnames(se),
      dex = factor(dex, levels = c("untrt", "trt")),
      cell = factor(cell),
      albut = factor(albut),
      run = factor(run)
    ) %>%
    dplyr::select(sample_id, cell, dex, albut, run)

  meta
}

map_gene_symbols <- function(gene_ids) {
  load_packages(c("AnnotationDbi", "org.Hs.eg.db"))

  clean_ids <- sub("\\..*$", "", gene_ids)
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = clean_ids,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = clean_ids,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )

  tibble(
    gene_id = gene_ids,
    ensembl_id = clean_ids,
    gene_symbol = unname(symbols),
    entrez_id = unname(entrez)
  )
}

annotate_result_frame <- function(df, gene_ids = rownames(df)) {
  annotations <- map_gene_symbols(gene_ids)
  tibble::as_tibble(df, rownames = "gene_id") %>%
    dplyr::left_join(annotations, by = "gene_id") %>%
    dplyr::relocate(ensembl_id, gene_symbol, entrez_id, .after = gene_id)
}

collapse_matrix_by_symbol <- function(expr_matrix, gene_ids = rownames(expr_matrix)) {
  annotations <- map_gene_symbols(gene_ids)
  keep <- !is.na(annotations$gene_symbol) & annotations$gene_symbol != ""

  collapsed <- tibble::as_tibble(expr_matrix[keep, , drop = FALSE], rownames = "gene_id") %>%
    dplyr::left_join(annotations[keep, c("gene_id", "gene_symbol")], by = "gene_id") %>%
    dplyr::select(-gene_id) %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), mean), .groups = "drop")

  out <- as.matrix(collapsed[, -1, drop = FALSE])
  rownames(out) <- collapsed$gene_symbol
  out
}

read_signature_panel <- function(path = "data/signatures/sample_states.tsv") {
  panel <- readr::read_tsv(path, show_col_types = FALSE)
  split(panel$gene, panel$signature)
}

paired_effect_table <- function(score_matrix, metadata, feature_name) {
  meta <- metadata %>%
    dplyr::select(sample_id, cell, dex) %>%
    dplyr::arrange(cell, dex)

  if (all(meta$sample_id %in% colnames(score_matrix))) {
    score_df <- tibble::as_tibble(t(score_matrix), rownames = "sample_id")
  } else if (all(meta$sample_id %in% rownames(score_matrix))) {
    score_df <- tibble::as_tibble(score_matrix, rownames = "sample_id")
  } else {
    stop("Could not align score matrix to sample metadata.")
  }

  score_df <- score_df %>%
    dplyr::left_join(meta, by = "sample_id")

  features <- setdiff(colnames(score_df), c("sample_id", "cell", "dex"))

  rows <- lapply(features, function(feature) {
    wide <- score_df %>%
      dplyr::select(cell, dex, dplyr::all_of(feature)) %>%
      dplyr::rename(score = dplyr::all_of(feature)) %>%
      tidyr::pivot_wider(names_from = dex, values_from = score)

    if (!all(c("untrt", "trt") %in% colnames(wide))) {
      stop("Expected paired untreated/treated columns for feature: ", feature)
    }

    delta <- wide$trt - wide$untrt
    test <- stats::t.test(wide$trt, wide$untrt, paired = TRUE)

    tibble(
      !!feature_name := feature,
      mean_untrt = mean(wide$untrt),
      mean_trt = mean(wide$trt),
      mean_delta = mean(delta),
      p_value = unname(test$p.value)
    )
  })

  dplyr::bind_rows(rows) %>%
    dplyr::mutate(padj = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(dplyr::desc(abs(mean_delta)))
}
