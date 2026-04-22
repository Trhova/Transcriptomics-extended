required_paths <- c(
  "results/tables/sample_metadata.csv",
  "results/tables/deseq2_results_summary.csv",
  "results/tables/top_de_genes.csv",
  "results/tables/top_ora_terms_up.csv",
  "results/tables/top_ora_terms_down.csv",
  "results/tables/top_gsea_terms.csv",
  "results/tables/sample_level_scores.csv",
  "results/tables/top_progeny_pathways.csv",
  "results/tables/top_viper_tfs.csv",
  "results/figures/sample_metadata_overview.png",
  "results/figures/library_sizes.png",
  "results/figures/pca_vst.png",
  "results/figures/sample_distance_heatmap.png",
  "results/figures/ma_plot.png",
  "results/figures/volcano_plot.png",
  "results/figures/ora_summary.png",
  "results/figures/gsea_summary.png",
  "results/figures/gsea_top_positive.png",
  "results/figures/gsea_top_negative.png",
  "results/figures/ssgsea_heatmap.png",
  "results/figures/progeny_activity.png",
  "results/figures/viper_tf_activity.png"
)

missing_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_paths) > 0) {
  stop(
    "Cannot render README. Missing outputs:\n",
    paste0("- ", missing_paths, collapse = "\n")
  )
}

if (rmarkdown::pandoc_available()) {
  rmarkdown::render(
    input = "README.Rmd",
    output_file = "README.md",
    quiet = TRUE
  )
} else {
  knitr::knit("README.Rmd", output = "README.md", quiet = TRUE)
}

message("Rendered README.md")
