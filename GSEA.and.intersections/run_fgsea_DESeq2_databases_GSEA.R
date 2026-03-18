# =============================================================================
# FGSEA using DESeq2 CSV results + GMT files in databases_GSEA
# =============================================================================
# - Reads DESeq2 CSVs (Gene_symbol, log2FoldChange, padj, ...) from project root
# - Ranks genes by log2FoldChange (or stat if you prefer)
# - Loads GMTs from ../databases_GSEA/ (Hallmark, KEGG, Reactome, WikiPathways)
# - Runs fgsea per cell type × pathway set; saves tables and plots
# =============================================================================

library(tidyverse)
library(fgsea)

# ----- Paths (run from scripts_FGSEA or project root) -----
base_dir <- if (dir.exists("databases_GSEA")) getwd() else dirname(getwd())
csv_dir <- base_dir
gmt_dir <- file.path(base_dir, "databases_GSEA")
out_dir <- file.path(base_dir, "GSEA_results_DESeq2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)

# ----- DESeq2 CSV files (cell-type names from filenames) -----
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
csv_files <- csv_files[basename(csv_files) != "summary_DEG_counts.csv"]
if (length(csv_files) == 0) {
  stop("No CSV files found in ", csv_dir)
}
cat("DESeq2 CSV files to process:", length(csv_files), "\n")
cat("GMT directory:", gmt_dir, "\n\n")

# ----- GMT files in databases_GSEA -----
gmt_sets <- list(
  Hallmark     = file.path(gmt_dir, "h.all.v2026.1.Hs.symbols.gmt"),
  KEGG         = file.path(gmt_dir, "c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt"),
  Reactome     = file.path(gmt_dir, "c2.cp.reactome.v2026.1.Hs.symbols.gmt"),
  WikiPathways = file.path(gmt_dir, "c2.cp.wikipathways.v2026.1.Hs.symbols.gmt"),
  GO_BP        = file.path(gmt_dir, "c5.go.bp.v2026.1.Hs.symbols.gmt")
)
gmt_sets <- gmt_sets[vapply(gmt_sets, file.exists, logical(1))]
if (length(gmt_sets) == 0) {
  stop("No GMT files found in ", gmt_dir)
}
cat("GMT sets to use:", paste(names(gmt_sets), collapse = ", "), "\n\n")

# ----- Process each DESeq2 CSV -----
for (csv_path in csv_files) {
  cell_type <- tools::file_path_sans_ext(basename(csv_path))
  cat("======== ", cell_type, " ========\n", sep = "")

  res <- read_csv(csv_path, col_types = cols(), show_col_types = FALSE)
  if (!"Gene_symbol" %in% names(res)) {
    cat("  Skip: no Gene_symbol column.\n\n")
    next
  }
  if (!"log2FoldChange" %in% names(res)) {
    cat("  Skip: no log2FoldChange column.\n\n")
    next
  }

  # Rank by log2FoldChange (DESeq2 style); drop NA and duplicate symbols
  res2 <- res %>%
    dplyr::select(Gene_symbol, log2FoldChange) %>%
    na.omit() %>%
    distinct() %>%
    group_by(Gene_symbol) %>%
    summarise(log2FoldChange = mean(log2FoldChange), .groups = "drop") %>%
    arrange(desc(log2FoldChange))
  ranks <- deframe(res2)
  cat("  Ranked genes:", length(ranks), "\n")

  for (set_name in names(gmt_sets)) {
    gmt_file <- gmt_sets[[set_name]]
    pathways <- gmtPathways(gmt_file)
    cat("  Running fgsea:", set_name, "... ")
    fgseaRes <- fgsea(pathways = pathways, stats = ranks, nperm = 1000)
    fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))

    out_table <- file.path(out_dir, "tables",
                           paste0(cell_type, "_", set_name, "_fgsea.csv"))
    write_csv(fgseaResTidy, out_table)

    # Barplot (top pathways by NES)
    n_plot <- min(30, nrow(fgseaResTidy))
    p <- ggplot(head(fgseaResTidy, n_plot), aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill = padj < 0.05)) +
      coord_flip() +
      labs(x = "", y = "NES", title = paste(cell_type, "-", set_name)) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8))
    ggsave(file.path(out_dir, "plots",
                     paste0(cell_type, "_", set_name, "_barplot.png")),
           p, width = 10, height = max(6, n_plot * 0.2), dpi = 150, bg = "white")
    cat("OK\n")
  }
  cat("\n")
}

cat("Done. Results in:", out_dir, "\n")
cat("  tables/: *_fgsea.csv\n")
cat("  plots/:  *_barplot.png\n")
