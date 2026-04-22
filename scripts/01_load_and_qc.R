source("scripts/utils.R")

check_required_packages(c("airway", "DESeq2", "pheatmap"))
load_packages(c("DESeq2", "pheatmap", "ggrepel", "scales"))
initialize_repo_dirs()

se <- read_airway_dataset()
metadata <- build_sample_metadata(se)
counts_mat <- SummarizedExperiment::assay(se)

dds <- DESeq2::DESeqDataSet(
  se,
  design = ~ cell + dex
)
dds$dex <- factor(dds$dex, levels = c("untrt", "trt"))
dds <- DESeq2::estimateSizeFactors(dds)
vsd <- DESeq2::vst(dds, blind = FALSE)

write_csv_table(metadata, "results/tables/sample_metadata.csv")
saveRDS(dds, "results/objects/dds_input.rds")
saveRDS(vsd, "results/objects/vsd.rds")

metadata_overview <- metadata %>%
  tidyr::pivot_longer(-sample_id, names_to = "field", values_to = "value") %>%
  dplyr::mutate(
    field = factor(field, levels = c("cell", "dex", "albut", "run")),
    sample_id = forcats::fct_rev(sample_id)
  ) %>%
  ggplot(aes(field, sample_id, fill = value)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = value), size = 3.2) +
  scale_fill_manual(
    values = stats::setNames(
      scales::hue_pal(l = 75, c = 90)(length(unique(c(
        as.character(metadata$cell),
        as.character(metadata$dex),
        as.character(metadata$albut),
        as.character(metadata$run)
      )))),
      unique(c(as.character(metadata$cell), as.character(metadata$dex), as.character(metadata$albut), as.character(metadata$run)))
    )
  ) +
  labs(
    title = "Airway sample metadata overview",
    subtitle = "Eight samples arranged as four donor-matched treated/untreated pairs",
    x = NULL,
    y = NULL
  ) +
  guides(fill = "none") +
  analysis_theme()

library_sizes <- tibble(
  sample_id = colnames(dds),
  library_size = colSums(counts_mat)
) %>%
  dplyr::left_join(metadata, by = "sample_id") %>%
  ggplot(aes(forcats::fct_reorder(sample_id, library_size), library_size, fill = dex)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = scales::comma(library_size)), hjust = -0.05, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("untrt" = "#7C98B3", "trt" = "#C1666B")) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Library sizes by sample",
    subtitle = "Depth is fairly balanced across the paired design",
    x = NULL,
    y = "Total mapped counts",
    fill = "Treatment"
  ) +
  analysis_theme()

pca_data <- DESeq2::plotPCA(vsd, intgroup = c("dex", "cell"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = dex, shape = cell, label = name)) +
  geom_point(size = 3.8) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c("untrt" = "#7C98B3", "trt" = "#C1666B")) +
  labs(
    title = "PCA on VST-transformed counts",
    subtitle = "Samples separate by both donor and dexamethasone exposure",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)"),
    color = "Treatment",
    shape = "Donor"
  ) +
  analysis_theme()

sample_dist <- dist(t(SummarizedExperiment::assay(vsd)))
sample_dist_mat <- as.matrix(sample_dist)
rownames(sample_dist_mat) <- metadata$sample_id
colnames(sample_dist_mat) <- metadata$sample_id
annotation_df <- as.data.frame(metadata[, c("cell", "dex")])
rownames(annotation_df) <- metadata$sample_id

png("results/figures/sample_distance_heatmap.png", width = 2200, height = 1800, res = 300)
pheatmap::pheatmap(
  sample_dist_mat,
  clustering_distance_rows = sample_dist,
  clustering_distance_cols = sample_dist,
  annotation_col = annotation_df,
  annotation_row = annotation_df,
  main = "Sample-to-sample distances on VST data",
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 9
)
dev.off()

save_plot(metadata_overview, "results/figures/sample_metadata_overview.png", width = 8.5, height = 4.5)
save_plot(library_sizes, "results/figures/library_sizes.png", width = 8.5, height = 4.8)
save_plot(pca_plot, "results/figures/pca_vst.png", width = 7.2, height = 5.6)
