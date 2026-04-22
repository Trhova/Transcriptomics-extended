source("scripts/utils.R")

check_required_packages(c("DESeq2", "apeglm"))
load_packages(c("DESeq2", "apeglm", "ggrepel", "scales"))
initialize_repo_dirs()

dds <- readRDS("results/objects/dds_input.rds")
dds <- DESeq2::DESeq(dds)

res <- DESeq2::results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05)
res_shrunk <- DESeq2::lfcShrink(
  dds,
  coef = "dex_trt_vs_untrt",
  type = "apeglm"
)

res_tbl <- annotate_result_frame(as.data.frame(res))
res_shrunk_tbl <- annotate_result_frame(as.data.frame(res_shrunk))

summary_tbl <- tibble(
  metric = c(
    "Genes tested",
    "Significant genes (padj < 0.05)",
    "Significant genes with |log2FC| >= 1",
    "Up-regulated genes",
    "Down-regulated genes"
  ),
  value = c(
    sum(!is.na(res_tbl$padj)),
    sum(res_tbl$padj < 0.05, na.rm = TRUE),
    sum(res_shrunk_tbl$padj < 0.05 & abs(res_shrunk_tbl$log2FoldChange) >= 1, na.rm = TRUE),
    sum(res_shrunk_tbl$padj < 0.05 & res_shrunk_tbl$log2FoldChange >= 1, na.rm = TRUE),
    sum(res_shrunk_tbl$padj < 0.05 & res_shrunk_tbl$log2FoldChange <= -1, na.rm = TRUE)
  )
)

top_de <- res_shrunk_tbl %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::select(gene_id, gene_symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
  dplyr::slice_head(n = 15)

write_csv_table(res_shrunk_tbl, "results/tables/deseq2_results_full.csv.gz")
write_csv_table(summary_tbl, "results/tables/deseq2_results_summary.csv")
write_csv_table(top_de, "results/tables/top_de_genes.csv")
saveRDS(dds, "results/objects/dds_fit.rds")
saveRDS(res_shrunk, "results/objects/res_shrunk.rds")

png("results/figures/ma_plot.png", width = 2200, height = 1700, res = 300)
DESeq2::plotMA(res_shrunk, ylim = c(-6, 6), main = "MA plot: dexamethasone-treated vs untreated")
dev.off()

volcano_tbl <- res_shrunk_tbl %>%
  dplyr::mutate(
    padj_plot = dplyr::if_else(is.na(padj), NA_real_, pmax(padj, .Machine$double.xmin)),
    neg_log10_padj = -log10(padj_plot),
    status = dplyr::case_when(
      padj < 0.05 & log2FoldChange >= 1 ~ "Up",
      padj < 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

labels <- volcano_tbl %>%
  dplyr::filter(status != "Not significant") %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 12)

volcano_plot <- ggplot(volcano_tbl, aes(log2FoldChange, neg_log10_padj, color = status)) +
  geom_point(alpha = 0.75, size = 1.6) +
  ggrepel::geom_text_repel(
    data = labels,
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#6C757D") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#6C757D") +
  scale_color_manual(
    values = c("Down" = "#406E8E", "Not significant" = "grey75", "Up" = "#C1666B")
  ) +
  labs(
    title = "Volcano plot of dexamethasone response",
    subtitle = "Shrunken log2 fold changes emphasize stable effect sizes",
    x = "Shrunken log2 fold change",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  analysis_theme()

save_plot(volcano_plot, "results/figures/volcano_plot.png", width = 7.5, height = 5.8)
