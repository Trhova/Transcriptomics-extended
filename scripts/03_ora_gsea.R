source("scripts/utils.R")

check_required_packages(c("fgsea", "msigdbr"))
load_packages(c("fgsea", "msigdbr", "scales"))
initialize_repo_dirs()

res_shrunk_tbl <- readr::read_csv("results/tables/deseq2_results_full.csv.gz", show_col_types = FALSE)

deg_threshold <- res_shrunk_tbl %>%
  dplyr::filter(!is.na(padj), !is.na(log2FoldChange), !is.na(gene_symbol), gene_symbol != "")

up_genes <- deg_threshold %>%
  dplyr::filter(padj < 0.05, log2FoldChange >= 1) %>%
  dplyr::pull(gene_symbol) %>%
  unique()

down_genes <- deg_threshold %>%
  dplyr::filter(padj < 0.05, log2FoldChange <= -1) %>%
  dplyr::pull(gene_symbol) %>%
  unique()

background <- deg_threshold$gene_symbol %>%
  unique()

go_bp_sets <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:BP"
) %>%
  dplyr::select(gs_name, gs_description, gene_symbol) %>%
  dplyr::distinct()

run_ora <- function(genes_of_interest, universe_genes, pathway_df, min_size = 10, max_size = 500) {
  terms <- split(pathway_df$gene_symbol, pathway_df$gs_name)
  descriptions <- pathway_df %>%
    dplyr::distinct(gs_name, gs_description)

  universe_genes <- unique(universe_genes)
  genes_of_interest <- unique(genes_of_interest)
  n_interest <- length(genes_of_interest)
  n_universe <- length(universe_genes)

  rows <- lapply(names(terms), function(term) {
    set_genes <- intersect(unique(terms[[term]]), universe_genes)
    set_size <- length(set_genes)

    if (set_size < min_size || set_size > max_size) {
      return(NULL)
    }

    overlap_genes <- intersect(genes_of_interest, set_genes)
    overlap <- length(overlap_genes)

    if (overlap == 0) {
      return(NULL)
    }

    tibble::tibble(
      ID = term,
      Description = descriptions$gs_description[descriptions$gs_name == term],
      GeneRatio = paste0(overlap, "/", n_interest),
      BgRatio = paste0(set_size, "/", n_universe),
      pvalue = stats::phyper(overlap - 1, set_size, n_universe - set_size, n_interest, lower.tail = FALSE),
      Count = overlap,
      geneID = paste(overlap_genes, collapse = "/")
    )
  })

  dplyr::bind_rows(rows) %>%
    dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
    dplyr::arrange(p.adjust, dplyr::desc(Count))
}

ora_up_tbl <- run_ora(up_genes, background, go_bp_sets) %>%
  dplyr::slice_head(n = 15)

ora_down_tbl <- run_ora(down_genes, background, go_bp_sets) %>%
  dplyr::slice_head(n = 15)

write_csv_table(ora_up_tbl, "results/tables/top_ora_terms_up.csv")
write_csv_table(ora_down_tbl, "results/tables/top_ora_terms_down.csv")

ora_plot_tbl <- dplyr::bind_rows(
  ora_up_tbl %>% dplyr::mutate(direction = "Up"),
  ora_down_tbl %>% dplyr::mutate(direction = "Down")
) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_head(n = 8) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(term = forcats::fct_reorder(Description, -log10(p.adjust)))

ora_plot <- ggplot(ora_plot_tbl, aes(-log10(p.adjust), term, color = direction, size = Count)) +
  geom_point(alpha = 0.85) +
  facet_wrap(~ direction, scales = "free_y") +
  scale_color_manual(values = c("Down" = "#406E8E", "Up" = "#C1666B")) +
  labs(
    title = "ORA highlights coherent biological themes among DE genes",
    subtitle = "Separate over-representation analyses for up- and down-regulated genes",
    x = "-log10 adjusted p-value",
    y = NULL,
    color = NULL
  ) +
  analysis_theme()

save_plot(ora_plot, "results/figures/ora_summary.png", width = 9, height = 5.6)

hallmark_sets <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

rank_tbl <- res_shrunk_tbl %>%
  dplyr::filter(!is.na(gene_symbol), gene_symbol != "", !is.na(log2FoldChange)) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

rank_vector <- rank_tbl$log2FoldChange
names(rank_vector) <- rank_tbl$gene_symbol
rank_vector <- sort(rank_vector, decreasing = TRUE)

pathway_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
gsea_res <- fgsea::fgsea(pathways = pathway_list, stats = rank_vector, minSize = 10, maxSize = 500)

gsea_tbl <- as_tibble(gsea_res) %>%
  dplyr::arrange(padj, dplyr::desc(abs(NES))) %>%
  dplyr::mutate(direction = dplyr::if_else(NES > 0, "Positive", "Negative"))

write_csv_table(gsea_tbl, "results/tables/top_gsea_terms.csv")

gsea_plot_tbl <- dplyr::bind_rows(
  gsea_tbl %>% dplyr::filter(NES > 0) %>% dplyr::slice_head(n = 8),
  gsea_tbl %>% dplyr::filter(NES < 0) %>% dplyr::slice_head(n = 8)
) %>%
  dplyr::mutate(pathway = forcats::fct_reorder(pathway, NES))

gsea_plot <- ggplot(gsea_plot_tbl, aes(NES, pathway, fill = direction)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = c("Negative" = "#406E8E", "Positive" = "#C1666B")) +
  labs(
    title = "GSEA tracks pathway shifts across the full ranked gene list",
    subtitle = "Positive NES indicates enrichment toward dexamethasone-up genes",
    x = "Normalized enrichment score",
    y = NULL,
    fill = NULL
  ) +
  analysis_theme()

save_plot(gsea_plot, "results/figures/gsea_summary.png", width = 8.5, height = 5.8)

best_positive <- gsea_tbl %>%
  dplyr::filter(NES > 0) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::pull(pathway)

best_negative <- gsea_tbl %>%
  dplyr::filter(NES < 0) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::pull(pathway)

png("results/figures/gsea_top_positive.png", width = 2200, height = 1700, res = 300)
plotEnrichment(pathway_list[[best_positive]], rank_vector) +
  labs(title = paste("Top positive GSEA term:", best_positive))
dev.off()

png("results/figures/gsea_top_negative.png", width = 2200, height = 1700, res = 300)
plotEnrichment(pathway_list[[best_negative]], rank_vector) +
  labs(title = paste("Top negative GSEA term:", best_negative))
dev.off()
