source("scripts/utils.R")

check_required_packages(c("fgsea", "msigdbr"))
load_packages(c("fgsea", "msigdbr", "scales"))
initialize_repo_dirs()

format_term_label <- function(value, prefix = NULL) {
  label <- value

  if (!is.null(prefix)) {
    label <- sub(prefix, "", label)
  }

  label <- gsub("_", " ", label)
  label <- tolower(label)
  label <- tools::toTitleCase(label)
  label <- gsub("\\bTnfa\\b", "TNFA", label)
  label <- gsub("\\bNfkb\\b", "NFkB", label)
  label <- gsub("\\bIl2\\b", "IL2", label)
  label <- gsub("\\bIl6\\b", "IL6", label)
  label <- gsub("\\bJak\\b", "JAK", label)
  label <- gsub("\\bStat3\\b", "STAT3", label)
  label <- gsub("\\bStat5\\b", "STAT5", label)
  label <- gsub("\\bMtorc1\\b", "MTORC1", label)
  label <- gsub("\\bE2f\\b", "E2F", label)
  label <- gsub("\\bKras\\b", "KRAS", label)
  label <- gsub("\\bUv\\b", "UV", label)
  label <- gsub("\\bDn\\b", "Down", label)
  label <- gsub("\\bUp\\b", "Up", label)
  label
}

wrap_term_label <- function(value, width = 30) {
  vapply(value, function(item) paste(strwrap(item, width = width), collapse = "\n"), character(1))
}

reorder_within <- function(x, by, within, fun = mean, sep = "___") {
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun)
}

scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(labels = function(x) gsub(paste0(sep, ".*$"), "", x), ...)
}

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

background <- deg_threshold$gene_symbol %>% unique()

hallmark_sets <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::distinct()

run_ora <- function(genes_of_interest, universe_genes, pathway_df, min_size = 10, max_size = 500) {
  terms <- split(pathway_df$gene_symbol, pathway_df$gs_name)

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
      Term = format_term_label(term, "^HALLMARK_"),
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

ora_up_tbl <- run_ora(up_genes, background, hallmark_sets) %>%
  dplyr::slice_head(n = 15)

ora_down_tbl <- run_ora(down_genes, background, hallmark_sets) %>%
  dplyr::slice_head(n = 15)

write_csv_table(ora_up_tbl, "results/tables/top_ora_terms_up.csv")
write_csv_table(ora_down_tbl, "results/tables/top_ora_terms_down.csv")

ora_plot_tbl <- dplyr::bind_rows(
  ora_up_tbl %>% dplyr::mutate(direction = "Up"),
  ora_down_tbl %>% dplyr::mutate(direction = "Down")
) %>%
  dplyr::mutate(direction = factor(direction, levels = c("Up", "Down"))) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_head(n = 6) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    significance = -log10(p.adjust),
    term_wrapped = wrap_term_label(Term, width = 28),
    term_ordered = reorder_within(term_wrapped, significance, direction)
  )

ora_plot <- ggplot(ora_plot_tbl, aes(significance, term_ordered, fill = direction)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(
    aes(label = Count),
    hjust = -0.2,
    size = 3.6,
    fontface = "bold",
    color = "#13293D"
  ) +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_reordered() +
  scale_fill_manual(values = c("Down" = "#406E8E", "Up" = "#C1666B")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    title = "ORA: top enriched Hallmark themes",
    subtitle = "Thresholded DE genes, analyzed separately for up- and down-regulated sets",
    x = "-log10 adjusted p-value",
    y = NULL
  ) +
  analysis_theme() +
  theme(
    axis.text.y = element_text(size = 11, lineheight = 0.95),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 13),
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 12)
  )

save_plot(ora_plot, "results/figures/ora_summary.png", width = 8.8, height = 7.6)

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
  dplyr::mutate(
    direction = dplyr::if_else(NES > 0, "Positive", "Negative"),
    Pathway = format_term_label(pathway, "^HALLMARK_")
  )

write_csv_table(gsea_tbl, "results/tables/top_gsea_terms.csv")

gsea_plot_tbl <- dplyr::bind_rows(
  gsea_tbl %>% dplyr::filter(NES > 0) %>% dplyr::slice_head(n = 6),
  gsea_tbl %>% dplyr::filter(NES < 0) %>% dplyr::slice_head(n = 6)
) %>%
  dplyr::mutate(pathway = forcats::fct_reorder(wrap_term_label(Pathway, width = 28), NES))

gsea_plot <- ggplot(gsea_plot_tbl, aes(NES, pathway, fill = direction)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = c("Negative" = "#406E8E", "Positive" = "#C1666B")) +
  labs(
    title = "GSEA: top enriched Hallmark programs",
    subtitle = "Ranked shrunken log2 fold changes; positive NES indicates dexamethasone-up programs",
    x = "Normalized enrichment score",
    y = NULL,
    fill = NULL
  ) +
  analysis_theme() +
  theme(
    axis.text.y = element_text(size = 11, lineheight = 0.95),
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 12)
  )

save_plot(gsea_plot, "results/figures/gsea_summary.png", width = 8.8, height = 6.4)

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
  labs(title = paste("Top positive GSEA term:", format_term_label(best_positive, "^HALLMARK_")))
dev.off()

png("results/figures/gsea_top_negative.png", width = 2200, height = 1700, res = 300)
plotEnrichment(pathway_list[[best_negative]], rank_vector) +
  labs(title = paste("Top negative GSEA term:", format_term_label(best_negative, "^HALLMARK_")))
dev.off()
