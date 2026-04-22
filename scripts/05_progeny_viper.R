source("scripts/utils.R")

check_required_packages(c("progeny", "dorothea", "viper"))
load_packages(c("progeny", "dorothea", "viper", "ggrepel"))
initialize_repo_dirs()

vsd <- readRDS("results/objects/vsd.rds")
metadata <- readr::read_csv("results/tables/sample_metadata.csv", show_col_types = FALSE)

vst_mat <- SummarizedExperiment::assay(vsd)
vst_symbol_mat <- collapse_matrix_by_symbol(vst_mat)

progeny_scores <- progeny::progeny(
  vst_symbol_mat,
  scale = TRUE,
  organism = "Human",
  top = 500
)

data(dorothea_hs, package = "dorothea")
regulon <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))

tf_scores <- dorothea::run_viper(
  vst_symbol_mat,
  regulon,
  options = list(method = "scale", minsize = 4, eset.filter = FALSE, verbose = FALSE)
)

progeny_tbl <- paired_effect_table(progeny_scores, metadata, "pathway")
viper_tbl <- paired_effect_table(tf_scores, metadata, "tf")

write_csv_table(progeny_tbl, "results/tables/top_progeny_pathways.csv")
write_csv_table(viper_tbl %>% dplyr::slice_head(n = 25), "results/tables/top_viper_tfs.csv")

progeny_plot <- progeny_tbl %>%
  dplyr::mutate(pathway = forcats::fct_reorder(pathway, mean_delta)) %>%
  ggplot(aes(mean_delta, pathway, fill = mean_delta > 0)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("TRUE" = "#C1666B", "FALSE" = "#406E8E")) +
  labs(
    title = "PROGENy pathway activity shifts with dexamethasone treatment",
    subtitle = "Bars show the paired treated-minus-untreated change across donors",
    x = "Mean paired activity change",
    y = NULL,
    fill = NULL
  ) +
  analysis_theme() +
  guides(fill = "none")

top_tf_plot_tbl <- viper_tbl %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::mutate(tf = forcats::fct_reorder(tf, mean_delta))

viper_plot <- ggplot(top_tf_plot_tbl, aes(mean_delta, tf, fill = mean_delta > 0)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("TRUE" = "#C1666B", "FALSE" = "#406E8E")) +
  labs(
    title = "DoRothEA + VIPER suggests upstream transcription factor activity changes",
    subtitle = "Top transcription factors ranked by absolute paired activity shift",
    x = "Mean paired activity change",
    y = NULL,
    fill = NULL
  ) +
  analysis_theme() +
  guides(fill = "none")

save_plot(progeny_plot, "results/figures/progeny_activity.png", width = 8.2, height = 5.4)
save_plot(viper_plot, "results/figures/viper_tf_activity.png", width = 8.2, height = 6.2)
