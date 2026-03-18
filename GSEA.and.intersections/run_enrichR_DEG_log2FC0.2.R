# =============================================================================
# Enrichr + barplots in ONE script (merged run_enrichR_DEG_log2FC0.2.R + run_enrichR_genes_in_common.R)
# =============================================================================
# Set MODE to:
#   "DEG"            -> tables/DEG_log2FC0.2/genes_only/*.txt
#   "genes_in_common" -> genes_in_common/*_genes_in_common.txt (also terms + Hallmark/KEGG folders)
# Run from project root:  Rscript scripts_FGSEA/run_enrichR_all_in_one.R
# =============================================================================

MODE <- "DEG"   # "DEG" or "genes_in_common"

library(enrichR)
library(readr)
library(dplyr)
library(ggplot2)

base_dir <- normalizePath(".", winslash = "/")

if (MODE == "DEG") {
  genes_dir <- file.path(base_dir, "tables", "DEG_log2FC0.2", "genes_only")
  out_dir <- file.path(base_dir, "tables", "DEG_log2FC0.2", "enrichr_results_R")
  input_pattern <- "\\.txt$"
  use_hall_kegg_dirs <- FALSE
} else {
  genes_dir <- file.path(base_dir, "genes_in_common")
  out_dir <- file.path(base_dir, "genes_in_common", "enrichr_results_R")
  input_pattern <- "_genes_in_common\\.txt$"
  use_hall_kegg_dirs <- TRUE
}

tables_dir <- file.path(out_dir, "tables")
plots_dir <- file.path(out_dir, "plots")
terms_dir <- file.path(out_dir, "terms_AdjP_below_0.2")
dir.create(out_dir, showWarnings = FALSE)
dir.create(tables_dir, showWarnings = FALSE)
dir.create(plots_dir, showWarnings = FALSE)
if (use_hall_kegg_dirs) {
  dir.create(terms_dir, showWarnings = FALSE)
  hall_dir <- file.path(genes_dir, "enrichr_Hallmark_only")
  kegg_dir <- file.path(genes_dir, "enrichr_KEGG_only")
  dir.create(hall_dir, showWarnings = FALSE)
  dir.create(kegg_dir, showWarnings = FALSE)
}

enrichr_databases <- c(
  "KEGG_2021_Human",
  "Reactome_2022",
  "MSigDB_Hallmark_2020",
  "WikiPathway_2023_Human",
  "GO_Biological_Process_2023"
)
adj_p_cutoff <- 0.2
top_n_plot <- 25

txt_files <- list.files(genes_dir, pattern = input_pattern, full.names = TRUE)
txt_files <- sort(txt_files)
if (length(txt_files) == 0) stop("No input files in ", genes_dir)

cat("MODE:", MODE, "\n")
cat("Input:", genes_dir, "\n")
cat("Output:", out_dir, "\n")
cat("Databases:", paste(enrichr_databases, collapse = ", "), "\n\n")

for (f in txt_files) {
  if (MODE == "genes_in_common") {
    base_name <- gsub("_genes_in_common\\.txt$", "", basename(f))
  } else {
    base_name <- sub("\\.txt$", "", basename(f))
  }
  cat("Processing:", base_name, "\n")
  gene_list <- unique(trimws(readLines(f)))
  gene_list <- gene_list[nzchar(gene_list)]
  if (length(gene_list) == 0) { cat("  Skip: no genes.\n"); next }
  cat("  Genes:", length(gene_list), "\n")

  enrichr_results <- tryCatch(enrichr(gene_list, enrichr_databases),
    error = function(e) { cat("  Error:", conditionMessage(e), "\n"); NULL })
  if (is.null(enrichr_results)) next

  if (use_hall_kegg_dirs && nrow(bind_rows(lapply(names(enrichr_results), function(db) {
    df <- enrichr_results[[db]]; if (nrow(df) == 0) return(NULL); df$Gene_set <- db; df
  }), .id = NULL)) > 0) {
    all_res <- bind_rows(lapply(names(enrichr_results), function(db) {
      df <- enrichr_results[[db]]
      if (nrow(df) == 0) return(NULL)
      df$Gene_set <- db
      df
    }))
    if ("Adjusted P-value" %in% names(all_res) && !"Adjusted.P.value" %in% names(all_res))
      all_res$Adjusted.P.value <- all_res[["Adjusted P-value"]]
    below <- all_res %>% filter(Adjusted.P.value < adj_p_cutoff) %>% arrange(Gene_set, Adjusted.P.value)
    if (nrow(below) > 0) {
      write_csv(below, file.path(terms_dir, paste0(base_name, "_terms_AdjP_below_0.2.csv")))
      below_hall <- below %>% filter(grepl("Hallmark", Gene_set, ignore.case = TRUE))
      if (nrow(below_hall) > 0)
        write_csv(below_hall, file.path(hall_dir, paste0(base_name, "_terms_AdjP_below_0.2_Hallmark_only.csv")))
      below_kegg <- below %>% filter(grepl("KEGG", Gene_set, ignore.case = TRUE))
      if (nrow(below_kegg) > 0)
        write_csv(below_kegg, file.path(kegg_dir, paste0(base_name, "_terms_AdjP_below_0.2_KEGG_only.csv")))
    }
  }

  for (db in names(enrichr_results)) {
    result_df <- enrichr_results[[db]]
    if (nrow(result_df) == 0) next
    if (!"Adjusted.P.value" %in% names(result_df) && "Adjusted P-value" %in% names(result_df))
      result_df$Adjusted.P.value <- result_df[["Adjusted P-value"]]
    significant <- result_df %>%
      filter(Adjusted.P.value < adj_p_cutoff) %>%
      arrange(Adjusted.P.value) %>%
      head(top_n_plot)
    if (nrow(significant) == 0) next

    db_short <- switch(db,
      "WikiPathway_2023_Human" = "WikiPathways",
      "MSigDB_Hallmark_2020" = if (use_hall_kegg_dirs) "MSigDB_Hallmark" else "Hallmark",
      "GO_Biological_Process_2023" = if (use_hall_kegg_dirs) "GO_Biological_Process" else "GO_BP",
      "Reactome_2022" = "Reactome",
      "KEGG_2021_Human" = "KEGG", db)

    write_csv(significant, file.path(tables_dir, paste0(base_name, "_", db_short, ".csv")))
    if (use_hall_kegg_dirs && db_short == "MSigDB_Hallmark")
      write_csv(significant, file.path(hall_dir, paste0(base_name, "_", db_short, ".csv")))
    if (use_hall_kegg_dirs && db_short == "KEGG")
      write_csv(significant, file.path(kegg_dir, paste0(base_name, "_", db_short, ".csv")))

    plot_data <- significant %>%
      mutate(
        Term_short = ifelse(nchar(Term) > 60, paste0(substr(Term, 1, 57), "..."), Term),
        log_p = -log10(pmax(Adjusted.P.value, 1e-300))
      )
    p <- ggplot(plot_data, aes(x = log_p, y = reorder(Term_short, log_p))) +
      geom_col(fill = "steelblue", alpha = 0.8) +
      geom_vline(xintercept = -log10(adj_p_cutoff), linetype = "dashed", color = "gray40") +
      labs(title = paste0(base_name, " — ", db_short),
           subtitle = paste0("Adj P < ", adj_p_cutoff, ", top ", nrow(plot_data), " terms"),
           x = "-log10(Adjusted P-value)", y = "") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA))
    plot_path <- file.path(plots_dir, paste0(base_name, "_", db_short, "_barplot.png"))
    ggsave(plot_path, p, width = 10, height = max(5, nrow(plot_data) * 0.3), dpi = 150, bg = "white")
    if (use_hall_kegg_dirs && db_short == "MSigDB_Hallmark")
      file.copy(plot_path, file.path(hall_dir, paste0(base_name, "_", db_short, "_barplot.png")), overwrite = TRUE)
    if (use_hall_kegg_dirs && db_short == "KEGG")
      file.copy(plot_path, file.path(kegg_dir, paste0(base_name, "_", db_short, "_barplot.png")), overwrite = TRUE)
    cat("  ", db_short, ":", nrow(significant), "terms -> table + barplot\n", sep = "")
  }
  cat("\n")
}

cat("Done. Results in:", out_dir, "\n")
cat("  tables/  plots/\n")
if (use_hall_kegg_dirs) cat("  terms_AdjP_below_0.2/  enrichr_Hallmark_only/  enrichr_KEGG_only/\n")
