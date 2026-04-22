source("scripts/utils.R")

check_required_packages(c("GSVA", "pheatmap"))
load_packages(c("GSVA", "pheatmap"))
initialize_repo_dirs()

vsd <- readRDS("results/objects/vsd.rds")
metadata <- readr::read_csv("results/tables/sample_metadata.csv", show_col_types = FALSE)

vst_mat <- SummarizedExperiment::assay(vsd)
vst_symbol_mat <- collapse_matrix_by_symbol(vst_mat)
signature_list <- read_signature_panel()

if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
  param <- GSVA::ssgseaParam(
    exprData = vst_symbol_mat,
    geneSets = signature_list,
    normalize = TRUE
  )
  ssgsea_scores <- GSVA::gsva(param, verbose = FALSE)
} else {
  ssgsea_scores <- GSVA::gsva(
    vst_symbol_mat,
    signature_list,
    method = "ssgsea",
    kcdf = "Gaussian",
    abs.ranking = TRUE,
    verbose = FALSE
  )
}

score_tbl <- as_tibble(t(ssgsea_scores), rownames = "sample_id") %>%
  dplyr::left_join(metadata, by = "sample_id")

write_csv_table(score_tbl, "results/tables/sample_level_scores.csv")

annotation_col <- as.data.frame(metadata[, c("cell", "dex")])
rownames(annotation_col) <- metadata$sample_id

png("results/figures/ssgsea_heatmap.png", width = 2400, height = 1800, res = 300)
pheatmap::pheatmap(
  ssgsea_scores,
  scale = "row",
  color = colorRampPalette(c("#406E8E", "white", "#C1666B"))(100),
  annotation_col = annotation_col,
  fontsize_row = 11,
  fontsize_col = 10,
  main = "Sample-level pathway and state activity by ssGSEA"
)
dev.off()
