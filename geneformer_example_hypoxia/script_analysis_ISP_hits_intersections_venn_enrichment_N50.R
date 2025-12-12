# ============================================================================
# COMBINED SCRIPT: ISP HITS ANALYSIS + EXCITATORY/IMM COMPARISON
# ============================================================================
# This script combines:
# 1. ISP_hits_model_18x3_read_all_csv_files_perform_enrichR_clusterProfiler_N50.R
# 2. compare_excitatory_imm_acute_chronic_with_enrichment_N50.R
#
# Features:
# - Reads all cell_type_ISP_* CSV files
# - Performs intersections for Imm, Olig, Glia, and Excitatory cell types
# - Creates Venn diagrams for ALL intersections
# - Prints and saves common genes for ALL intersections
# - Runs enrichment analyses (enrichR and clusterProfiler) on all intersections
#
# Filter criteria for significant genes:
#   - Sig == 1
#   - N_Detections >= 50
#   - Shift_to_goal_end > 0 (positive)
#   - Shift_to_alt_end < 0 (negative)
# ============================================================================

library(readr)
library(dplyr)
library(ggplot2)

# Load clusterProfiler packages (install if needed)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  cat("Note: clusterProfiler package not found. Install with: BiocManager::install('clusterProfiler')\n")
  cat("clusterProfiler analyses will be skipped.\n")
  clusterprofiler_available <- FALSE
} else {
  library(clusterProfiler)
  clusterprofiler_available <- TRUE
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  cat("Note: org.Hs.eg.db package not found. Install with: BiocManager::install('org.Hs.eg.db')\n")
  orgdb_available <- FALSE
} else {
  library(org.Hs.eg.db)
  orgdb_available <- TRUE
}

if (!requireNamespace("ReactomePA", quietly = TRUE)) {
  cat("Note: ReactomePA package not found. Install with: BiocManager::install('ReactomePA')\n")
  reactomepa_available <- FALSE
} else {
  library(ReactomePA)
  reactomepa_available <- TRUE
}

if (!requireNamespace("enrichplot", quietly = TRUE)) {
  cat("Note: enrichplot package not found. Install with: BiocManager::install('enrichplot')\n")
  enrichplot_available <- FALSE
} else {
  library(enrichplot)
  enrichplot_available <- TRUE
}

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  cat("Note: msigdbr package not found. Install with: install.packages('msigdbr')\n")
  msigdbr_available <- FALSE
} else {
  library(msigdbr)
  msigdbr_available <- TRUE
}

# Load enrichR (install if needed)
if (!requireNamespace("enrichR", quietly = TRUE)) {
  cat("Note: enrichR package not found. Install with: install.packages('enrichR')\n")
  cat("enrichR analyses will be skipped.\n")
  enrichr_available <- FALSE
} else {
  library(enrichR)
  enrichr_available <- TRUE
}

# Load Venn diagram packages (install if needed)
if (!requireNamespace("ggvenn", quietly = TRUE)) {
  cat("Note: ggvenn package not found. Install with: install.packages('ggvenn')\n")
  cat("Venn diagrams will be skipped.\n")
  ggvenn_available <- FALSE
} else {
  library(ggvenn)
  ggvenn_available <- TRUE
}

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  cat("Note: UpSetR package not found. Install with: install.packages('UpSetR')\n")
  cat("UpSet plots will be skipped.\n")
  upsetr_available <- FALSE
} else {
  library(UpSetR)
  upsetr_available <- TRUE
}

# ===== SETUP =====
script_dir <- getwd()
cat("Working directory:", script_dir, "\n\n")

# ===== HELPER FUNCTION TO FILTER SIGNIFICANT GENES =====
# Filter significant genes with all criteria:
# - Sig == 1
# - N_Detections >= 50
# - Shift_to_goal_end > 0 (positive)
# - Shift_to_alt_end < 0 (negative)
filter_significant_genes <- function(df) {
  if (!"Sig" %in% colnames(df)) {
    return(df %>% filter(FALSE))  # Return empty if Sig column doesn't exist
  }
  
  # Find the alt_end column name (varies between acute and chronic)
  alt_end_col <- grep("^Shift_to_alt_end", colnames(df), value = TRUE)
  
  if (length(alt_end_col) == 0) {
    # If no alt_end column, just filter with available criteria
    df_filtered <- df %>%
      filter(Sig == 1,
             N_Detections >= 50,
             Shift_to_goal_end > 0)
  } else {
    # Filter with all criteria including alt_end
    alt_end_name <- alt_end_col[1]
    df_filtered <- df %>%
      filter(Sig == 1,
             N_Detections >= 50,
             Shift_to_goal_end > 0,
             .data[[alt_end_name]] < 0)
  }
  
  return(df_filtered)
}

# ===== FUNCTION TO EXTRACT SIGNIFICANT GENES FROM DATAFRAME =====
extract_significant_genes <- function(df) {
  # First filter for significant genes with all criteria
  df_sig <- filter_significant_genes(df)
  
  # Then extract gene names
  if ("Gene_name" %in% colnames(df_sig)) {
    genes <- unique(df_sig$Gene_name[!is.na(df_sig$Gene_name) & df_sig$Gene_name != ""])
    return(genes)
  } else {
    return(character(0))
  }
}

# ===== FUNCTION TO CREATE VENN DIAGRAMS =====
# Function to create Venn diagrams from significant genes in multiple variables
create_venn_diagram <- function(var_names, title = "Venn Diagram", output_file = NULL, use_sig_only = TRUE) {
  # Check if ggvenn is available
  if (!ggvenn_available) {
    cat("Warning: ggvenn package not available. Skipping Venn diagram.\n")
    return(NULL)
  }
  
  # Extract significant genes from each variable
  gene_lists <- list()
  valid_vars <- c()
  
  for (var_name in var_names) {
    if (exists(var_name, envir = .GlobalEnv)) {
      df <- get(var_name, envir = .GlobalEnv)
      
      if (use_sig_only && "Sig" %in% colnames(df) && "N_Detections" %in% colnames(df) && 
          "Shift_to_goal_end" %in% colnames(df)) {
        df <- filter_significant_genes(df)
      }
      
      if ("Gene_name" %in% colnames(df)) {
        genes <- unique(df$Gene_name[!is.na(df$Gene_name) & df$Gene_name != ""])
        gene_lists[[var_name]] <- genes
        valid_vars <- c(valid_vars, var_name)
      }
    }
  }
  
  if (length(gene_lists) < 2) {
    cat("Warning: Need at least 2 valid variables to create Venn diagram.\n")
    return(NULL)
  }
  
  # For 2-3 sets, use ggvenn
  if (length(gene_lists) <= 3) {
    tryCatch({
      p <- ggvenn(gene_lists, 
                  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")[seq_along(gene_lists)],
                  stroke_size = 0.5,
                  set_name_size = 4,
                  text_size = 3,
                  show_percentage = TRUE) +
        labs(title = title) +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      
      # Save plot if output file is specified
      if (!is.null(output_file)) {
        ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
        cat("Venn diagram saved to:", output_file, "\n")
      }
      
      return(p)
    }, error = function(e) {
      cat("Error creating Venn diagram:", conditionMessage(e), "\n")
      return(NULL)
    })
  } else {
    # For more than 3 sets, suggest UpSet plot
    cat("Note: More than 3 sets detected. Consider using UpSet plot instead.\n")
    if (upsetr_available) {
      tryCatch({
        # Convert to binary matrix for UpSetR
        all_genes <- unique(unlist(gene_lists))
        binary_matrix <- data.frame(
          Gene = all_genes,
          stringsAsFactors = FALSE
        )
        
        for (var_name in names(gene_lists)) {
          binary_matrix[[var_name]] <- as.integer(all_genes %in% gene_lists[[var_name]])
        }
        
        # Create UpSet plot
        if (!is.null(output_file)) {
          png(output_file, width = 1200, height = 800)
          print(upset(binary_matrix, sets = names(gene_lists), order.by = "freq"))
          dev.off()
          cat("UpSet plot saved to:", output_file, "\n")
        }
      }, error = function(e) {
        cat("Error creating UpSet plot:", conditionMessage(e), "\n")
        return(NULL)
      })
    }
  }
}

# ===== FUNCTION TO CREATE VENN DIAGRAM FROM GENE LISTS (not variables) =====
create_venn_diagram_from_lists <- function(gene_lists, title = "Venn Diagram", output_file = NULL) {
  # Check if ggvenn is available
  if (!ggvenn_available) {
    cat("Warning: ggvenn package not available. Skipping Venn diagram.\n")
    return(NULL)
  }
  
  if (length(gene_lists) < 2) {
    cat("Warning: Need at least 2 gene lists to create Venn diagram.\n")
    return(NULL)
  }
  
  # For 2-3 sets, use ggvenn
  if (length(gene_lists) <= 3) {
    tryCatch({
      p <- ggvenn(gene_lists, 
                  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")[seq_along(gene_lists)],
                  stroke_size = 0.5,
                  set_name_size = 4,
                  text_size = 3,
                  show_percentage = TRUE) +
        labs(title = title) +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      
      # Save plot if output file is specified
      if (!is.null(output_file)) {
        ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
        cat("Venn diagram saved to:", output_file, "\n")
      }
      
      return(p)
    }, error = function(e) {
      cat("Error creating Venn diagram:", conditionMessage(e), "\n")
      return(NULL)
    })
  } else {
    # For more than 3 sets, suggest UpSet plot
    cat("Note: More than 3 sets detected. Consider using UpSet plot instead.\n")
  }
}

# ===== HELPER FUNCTION TO CALCULATE JACCARD SIMILARITY =====
calculate_jaccard <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union_size <- length(union(set1, set2))
  if (union_size == 0) return(0)
  return(intersection / union_size)
}

# ===== FUNCTION TO RUN enrichR ANALYSIS =====
run_enrichr_analysis <- function(gene_list, output_name, results_base_dir = NULL, adj_p_cutoff = 0.2, top_n = 50) {
  # Check if enrichR is available
  if (!enrichr_available || length(gene_list) == 0) {
    return(NULL)
  }
  
  cat("  Running enrichR for:", output_name, "(", length(gene_list), "genes)...\n")
  
  # Clean gene list
  gene_list <- as.character(gene_list)
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  gene_list <- unique(gene_list)
  
  if (length(gene_list) == 0) {
    cat("Warning: Empty gene list provided for", output_name, "\n")
    return(NULL)
  }
  
  # Set results directory
  if (is.null(results_base_dir)) {
    results_base_dir <- script_dir
  }
  results_dir <- file.path(results_base_dir, paste0(output_name, "_enrichR_results"))
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(results_dir, "plots"), showWarnings = FALSE)
  dir.create(file.path(results_dir, "tables"), showWarnings = FALSE)
  
  # Define enrichR databases
  enrichr_databases <- c(
    "GO_Biological_Process_2021", "GO_Molecular_Function_2021", 
    "GO_Cellular_Component_2021", "KEGG_2021_Human", 
    "Reactome_2022", "WikiPathways_2021_Human"
  )
  
  # Run enrichR analysis
  tryCatch({
    enrichr_results <- enrichr(gene_list, enrichr_databases)
    
    # Save gene list
    write_csv(data.frame(Gene_name = gene_list), 
              file.path(results_dir, "tables", paste0(output_name, "_gene_list.csv")))
    
    # Process each database
    all_significant_results <- data.frame()
    summary_stats <- data.frame()
    
    for (database in names(enrichr_results)) {
      result_df <- enrichr_results[[database]]
      
      if (nrow(result_df) > 0) {
        # Filter significant results
        significant <- result_df %>%
          filter(Adjusted.P.value < adj_p_cutoff) %>%
          arrange(Adjusted.P.value) %>%
          head(top_n)
        
        if (nrow(significant) > 0) {
          # Add database info
          significant$Database <- database
          
          # Save individual database results
          write_csv(significant, 
                    file.path(results_dir, "tables", paste0(output_name, "_", database, ".csv")))
          
          # Add to combined results (top 5 per database)
          top5 <- head(significant, 5)
          all_significant_results <- rbind(all_significant_results, top5)
          
          # Create plot
          plot_data <- significant %>%
            head(15) %>%
            mutate(
              Term_short = ifelse(nchar(Term) > 50, paste0(substr(Term, 1, 47), "..."), Term),
              log_p = -log10(Adjusted.P.value)
            )
          
          # Make plot
          p <- ggplot(plot_data, aes(x = log_p, y = reorder(Term_short, log_p))) +
            geom_col(fill = "steelblue", alpha = 0.7) +
            labs(
              title = paste("Enriched Terms -", database),
              subtitle = paste(output_name, "(", length(gene_list), "genes)"),
              x = "-log10(Adjusted P-value)",
              y = ""
            ) +
            theme_minimal() +
            theme(
              axis.text.y = element_text(size = 9),
              plot.background = element_rect(fill = "white", color = NA),
              panel.background = element_rect(fill = "white", color = NA),
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5)
            )
          
          # Save plot
          ggsave(file.path(results_dir, "plots", paste0(output_name, "_", database, ".png")), 
                 p, width = 10, height = 6, dpi = 300, bg = "white")
          
          # Add to summary
          summary_stats <- rbind(summary_stats, data.frame(
            Database = database,
            Significant_Terms = nrow(significant),
            Top_Term = significant$Term[1],
            Top_P_Value = significant$Adjusted.P.value[1],
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Save combined results
    if (nrow(all_significant_results) > 0) {
      write_csv(all_significant_results, 
                file.path(results_dir, "tables", paste0(output_name, "_all_top_results.csv")))
    }
    
    # Save summary
    if (nrow(summary_stats) > 0) {
      write_csv(summary_stats, 
                file.path(results_dir, "tables", paste0(output_name, "_summary.csv")))
      cat("    Found", nrow(summary_stats), "databases with significant results\n")
    } else {
      cat("    No significant results found\n")
    }
    
    cat("    enrichR results saved to:", results_dir, "\n")
    return(list(results = enrichr_results, summary = summary_stats, directory = results_dir))
    
  }, error = function(e) {
    cat("    Error in enrichR:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# ===== FUNCTION TO RUN clusterProfiler ANALYSIS =====
run_clusterprofiler_analysis <- function(gene_list, output_name, results_base_dir = NULL, 
                                         pvalue_cutoff = 0.2, qvalue_cutoff = 1.0, 
                                         min_gssize = 10, max_gssize = 500) {
  # Check if required packages are available
  if (!clusterprofiler_available || !orgdb_available || length(gene_list) == 0) {
    return(NULL)
  }
  
  cat("  Running clusterProfiler for:", output_name, "(", length(gene_list), "genes)...\n")
  
  # Clean gene list
  gene_list <- as.character(gene_list)
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  gene_list <- unique(gene_list)
  
  if (length(gene_list) == 0) {
    cat("Warning: Empty gene list provided for", output_name, "\n")
    return(NULL)
  }
  
  # Set results directory
  if (is.null(results_base_dir)) {
    results_base_dir <- script_dir
  }
  results_dir <- file.path(results_base_dir, paste0(output_name, "_clusterProfiler_results"))
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(results_dir, "plots"), showWarnings = FALSE)
  dir.create(file.path(results_dir, "tables"), showWarnings = FALSE)
  
  # Convert gene names to Entrez IDs
  tryCatch({
    gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    if (nrow(gene_ids) == 0) {
      cat("    Warning: Could not convert any genes to Entrez IDs\n")
      return(NULL)
    }
    
    cat("    Converted", nrow(gene_ids), "genes to Entrez IDs\n")
    
    all_results <- list()
    
    # ===== GO Over-Representation Analysis =====
    cat("    Running GO ORA...\n")
    tryCatch({
      ego <- enrichGO(gene = gene_ids$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      keyType = "ENTREZID",
                      pAdjustMethod = "BH",
                      pvalueCutoff = pvalue_cutoff,
                      qvalueCutoff = qvalue_cutoff,
                      minGSSize = min_gssize,
                      maxGSSize = max_gssize,
                      readable = TRUE)
      
      if (!is.null(ego) && nrow(ego@result) > 0) {
        write.table(ego@result,
                    file = file.path(results_dir, "tables", paste0(output_name, "_GO_ORA_Results.txt")),
                    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        
        png(file.path(results_dir, "plots", paste0(output_name, "_GO_ORA.png")),
            width = 1000, height = 800)
        print(dotplot(ego, showCategory = 20))
        dev.off()
        
        all_results$GO_ORA <- ego
        cat("      Found", nrow(ego@result), "enriched GO terms\n")
      }
    }, error = function(e) {
      cat("      Error in GO ORA:", conditionMessage(e), "\n")
    })
    
    # ===== KEGG Over-Representation Analysis =====
    cat("    Running KEGG ORA...\n")
    tryCatch({
      kegg_enrich <- enrichKEGG(gene = gene_ids$ENTREZID,
                                organism = "hsa",
                                pAdjustMethod = "BH",
                                pvalueCutoff = pvalue_cutoff,
                                minGSSize = min_gssize,
                                maxGSSize = max_gssize)
      
      if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
        kegg_result <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
        write.table(kegg_result@result,
                    file = file.path(results_dir, "tables", paste0(output_name, "_KEGG_ORA_Results.txt")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        
        png(file.path(results_dir, "plots", paste0(output_name, "_KEGG_ORA.png")),
            width = 600, height = 600)
        print(dotplot(kegg_enrich, showCategory = 20))
        dev.off()
        
        all_results$KEGG_ORA <- kegg_enrich
        cat("      Found", nrow(kegg_enrich@result), "enriched KEGG pathways\n")
      }
    }, error = function(e) {
      cat("      Error in KEGG ORA:", conditionMessage(e), "\n")
    })
    
    # ===== Reactome Over-Representation Analysis =====
    if (reactomepa_available) {
      cat("    Running Reactome ORA...\n")
      tryCatch({
        reactome_ora <- enrichPathway(gene = gene_ids$ENTREZID,
                                      organism = "human",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = pvalue_cutoff,
                                      minGSSize = min_gssize,
                                      maxGSSize = max_gssize)
        
        if (!is.null(reactome_ora) && nrow(reactome_ora@result) > 0) {
          reactome_result <- setReadable(reactome_ora, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
          write.table(reactome_result@result,
                      file = file.path(results_dir, "tables", paste0(output_name, "_Reactome_ORA_Results.txt")),
                      sep = "\t", row.names = FALSE, quote = FALSE)
          
          png(file.path(results_dir, "plots", paste0(output_name, "_Reactome_ORA.png")),
              width = 1000, height = 800)
          print(dotplot(reactome_ora, showCategory = 20))
          dev.off()
          
          all_results$Reactome_ORA <- reactome_ora
          cat("      Found", nrow(reactome_ora@result), "enriched Reactome pathways\n")
        }
      }, error = function(e) {
        cat("      Error in Reactome ORA:", conditionMessage(e), "\n")
      })
    }
    
    # ===== WikiPathways Over-Representation Analysis =====
    cat("    Running WikiPathways ORA...\n")
    tryCatch({
      wikipathways_enrich <- enrichWP(gene = gene_ids$ENTREZID,
                                       organism = "Homo sapiens",
                                       pvalueCutoff = pvalue_cutoff,
                                       minGSSize = min_gssize,
                                       maxGSSize = max_gssize)
      
      if (!is.null(wikipathways_enrich) && nrow(wikipathways_enrich@result) > 0) {
        wikipathways_result <- setReadable(wikipathways_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
        write.table(wikipathways_result@result,
                    file = file.path(results_dir, "tables", paste0(output_name, "_WikiPathways_ORA_Results.txt")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        
        png(file.path(results_dir, "plots", paste0(output_name, "_WikiPathways_ORA.png")),
            width = 1000, height = 800)
        print(dotplot(wikipathways_enrich, showCategory = 20))
        dev.off()
        
        all_results$WikiPathways_ORA <- wikipathways_enrich
        cat("      Found", nrow(wikipathways_enrich@result), "enriched WikiPathways\n")
      }
    }, error = function(e) {
      cat("      Error in WikiPathways ORA:", conditionMessage(e), "\n")
    })
    
    # ===== MSigDB Hallmark Over-Representation Analysis =====
    if (msigdbr_available) {
      cat("    Running MSigDB Hallmark ORA...\n")
      tryCatch({
        hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
        hallmark_t2g <- hallmark_sets %>%
          dplyr::select(gs_name, entrez_gene)
        
        hallmark_enrich <- enricher(gene = gene_ids$ENTREZID,
                                    TERM2GENE = hallmark_t2g,
                                    pvalueCutoff = pvalue_cutoff,
                                    pAdjustMethod = "BH",
                                    minGSSize = min_gssize,
                                    maxGSSize = max_gssize)
        
        if (!is.null(hallmark_enrich) && nrow(hallmark_enrich@result) > 0) {
          hallmark_result <- setReadable(hallmark_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
          write.table(hallmark_result@result,
                      file = file.path(results_dir, "tables", paste0(output_name, "_MSigDB_Hallmark_ORA_Results.txt")),
                      sep = "\t", row.names = FALSE, quote = FALSE)
          
          png(file.path(results_dir, "plots", paste0(output_name, "_MSigDB_Hallmark_ORA.png")),
              width = 1000, height = 800)
          print(dotplot(hallmark_enrich, showCategory = 20))
          dev.off()
          
          all_results$MSigDB_Hallmark_ORA <- hallmark_enrich
          cat("      Found", nrow(hallmark_enrich@result), "enriched Hallmark gene sets\n")
        }
      }, error = function(e) {
        cat("      Error in MSigDB Hallmark ORA:", conditionMessage(e), "\n")
      })
    }
    
    cat("    Results saved to:", results_dir, "\n")
    return(list(results = all_results, directory = results_dir))
    
  }, error = function(e) {
    cat("    Error running clusterProfiler:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# ===== READ ALL CSV FILES =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("READING ALL CSV FILES\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")

# Find all CSV files matching the pattern
csv_files <- list.files(path = script_dir, pattern = "^cell_type_ISP_.*\\.csv$", full.names = TRUE)
csv_files <- csv_files[order(basename(csv_files))]

cat("Found", length(csv_files), "CSV files\n\n")

# Create vectors to track variable names and corresponding files
created_variables <- c()
successful_files <- c()
all_data <- list()  # Store all dataframes by filename
all_genes <- list()  # Store all gene lists by filename

# Read all CSV files
for (file in csv_files) {
  file_name <- basename(file)
  cat("Reading:", file_name, "...")
  
  # Read the CSV file
  tryCatch({
    df <- read_csv(file, col_types = cols(), show_col_types = FALSE)
    
    # Create a valid R variable name from the file name
    var_name <- gsub("\\.csv$", "", file_name)
    var_name <- gsub("^cell_type_ISP_", "", var_name)
    var_name <- gsub("_GPU_ONLY_genes$", "", var_name)
    var_name <- gsub("-", "_", var_name)
    var_name <- gsub("[^A-Za-z0-9_]", "_", var_name)
    if (grepl("^[0-9]", var_name)) {
      var_name <- paste0("X", var_name)
    }
    
    # Assign the data frame to a variable in the global environment
    assign(var_name, df, envir = .GlobalEnv)
    created_variables <- c(created_variables, var_name)
    successful_files <- c(successful_files, file_name)
    all_data[[file_name]] <- df
    
    # Extract significant genes
    genes <- extract_significant_genes(df)
    all_genes[[file_name]] <- genes
    
    cat(" OK -> variable:", var_name, "(", nrow(df), "rows,", length(genes), "significant genes)\n")
    
  }, error = function(e) {
    cat(" ERROR:", conditionMessage(e), "\n")
  })
}

cat("\nSuccessfully read", length(created_variables), "files\n\n")

# ===== DEFINE FILES TO COMPARE (for excitatory/imm comparison) =====
acute_files <- c(
  "cell_type_ISP_excitatory_acute_GPU_ONLY_genes.csv",
  "cell_type_ISP_Imm-CGE_acute_GPU_ONLY_genes.csv",
  "cell_type_ISP_Imm-LGE_acute_GPU_ONLY_genes.csv",
  "cell_type_ISP_Imm-MGE_acute_GPU_ONLY_genes.csv"
)

chronic_files <- c(
  "cell_type_ISP_excitatory_chronic_GPU_ONLY_genes.csv",
  "cell_type_ISP_Imm-CGE_chronic_GPU_ONLY_genes.csv",
  "cell_type_ISP_Imm-LGE_chronic_GPU_ONLY_genes.csv",
  "cell_type_ISP_Imm-MGE_chronic_GPU_ONLY_genes.csv"
)

# Create output directory
output_dir <- file.path(script_dir, "combined_analysis_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("Results will be saved to:", output_dir, "\n\n")

# ===== INTERSECTION ANALYSIS: Imm Acute (CGE, LGE, MGE) =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION ANALYSIS: Imm Acute (CGE, LGE, MGE)\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")

imm_vars_acute <- c("Imm_CGE_acute", "Imm_LGE_acute", "Imm_MGE_acute")
imm_sig_genes_acute <- list()

for (var_name in imm_vars_acute) {
  if (exists(var_name, envir = .GlobalEnv)) {
    df <- get(var_name, envir = .GlobalEnv)
    genes <- extract_significant_genes(df)
    imm_sig_genes_acute[[var_name]] <- genes
    cat(var_name, ":", length(genes), "significant genes\n")
  }
}

if (length(imm_sig_genes_acute) >= 2) {
  # Find intersections
  all_genes_list <- unlist(imm_sig_genes_acute)
  unique_genes <- unique(all_genes_list)
  
  # Pairwise intersections
  if (length(imm_sig_genes_acute) >= 2) {
    var_names <- names(imm_sig_genes_acute)
    n_vars <- length(var_names)
    
    # Calculate all pairwise intersections
    pairwise_intersections <- list()
    for (i in 1:(n_vars - 1)) {
      for (j in (i + 1):n_vars) {
        var1 <- var_names[i]
        var2 <- var_names[j]
        common <- intersect(imm_sig_genes_acute[[var1]], imm_sig_genes_acute[[var2]])
        pair_name <- paste0(gsub("Imm_", "", var1), "_", gsub("Imm_|_acute", "", var2))
        pairwise_intersections[[pair_name]] <- common
        cat("  ", var1, "∩", var2, ":", length(common), "genes\n")
        
        # Print common genes
        if (length(common) > 0) {
          cat("    Common genes:", paste(head(common, 10), collapse = ", "))
          if (length(common) > 10) cat(" ... (", length(common), "total)")
          cat("\n")
          
          # Save to CSV
          output_file <- file.path(script_dir, paste0("imm_acute_", toupper(gsub("Imm_|_acute", "", var1)), "_", toupper(gsub("Imm_|_acute", "", var2)), "_intersection.csv"))
          write_csv(data.frame(Gene_name = common, stringsAsFactors = FALSE), output_file)
          cat("    Saved to:", basename(output_file), "\n")
        }
      }
    }
    
    # Three-way intersection (if we have 3 variables)
    if (n_vars == 3) {
      common_all_three <- Reduce(intersect, imm_sig_genes_acute)
      cat("  Common to all three:", length(common_all_three), "genes\n")
      
      if (length(common_all_three) > 0) {
        cat("    Common genes:", paste(head(common_all_three, 10), collapse = ", "))
        if (length(common_all_three) > 10) cat(" ... (", length(common_all_three), "total)")
        cat("\n")
        
        # Save to CSV
        output_file <- file.path(script_dir, "imm_acute_common_all_three_genes.csv")
        write_csv(data.frame(Gene_name = common_all_three, stringsAsFactors = FALSE), output_file)
        cat("    Saved to:", basename(output_file), "\n")
        
        assign("imm_common_all_three", common_all_three, envir = .GlobalEnv)
      }
      
      # Create Venn diagram
      cat("\n  Creating Venn diagram...\n")
      venn_file <- file.path(script_dir, "imm_acute_venn_diagram.png")
      create_venn_diagram(imm_vars_acute, 
                         title = "Imm Cell Types - Acute (Significant Genes)",
                         output_file = venn_file,
                         use_sig_only = TRUE)
    }
  }
  cat("\n")
} else {
  cat("Warning: Need at least 2 Imm acute variables for intersection analysis\n\n")
}

# ===== INTERSECTION ANALYSIS: Imm Chronic (CGE, LGE, MGE) =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION ANALYSIS: Imm Chronic (CGE, LGE, MGE)\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")

imm_vars_chronic <- c("Imm_CGE_chronic", "Imm_LGE_chronic", "Imm_MGE_chronic")
imm_sig_genes_chronic <- list()

for (var_name in imm_vars_chronic) {
  if (exists(var_name, envir = .GlobalEnv)) {
    df <- get(var_name, envir = .GlobalEnv)
    genes <- extract_significant_genes(df)
    imm_sig_genes_chronic[[var_name]] <- genes
    cat(var_name, ":", length(genes), "significant genes\n")
  }
}

if (length(imm_sig_genes_chronic) >= 2) {
  # Find intersections (same logic as acute)
  var_names <- names(imm_sig_genes_chronic)
  n_vars <- length(var_names)
  
  # Calculate all pairwise intersections
  for (i in 1:(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      var1 <- var_names[i]
      var2 <- var_names[j]
      common <- intersect(imm_sig_genes_chronic[[var1]], imm_sig_genes_chronic[[var2]])
      cat("  ", var1, "∩", var2, ":", length(common), "genes\n")
      
      if (length(common) > 0) {
        cat("    Common genes:", paste(head(common, 10), collapse = ", "))
        if (length(common) > 10) cat(" ... (", length(common), "total)")
        cat("\n")
        
        # Save to CSV
        output_file <- file.path(script_dir, paste0("imm_chronic_", toupper(gsub("Imm_|_chronic", "", var1)), "_", toupper(gsub("Imm_|_chronic", "", var2)), "_intersection.csv"))
        write_csv(data.frame(Gene_name = common, stringsAsFactors = FALSE), output_file)
        cat("    Saved to:", basename(output_file), "\n")
      }
    }
  }
  
  # Three-way intersection (if we have 3 variables)
  if (n_vars == 3) {
    common_all_three <- Reduce(intersect, imm_sig_genes_chronic)
    cat("  Common to all three:", length(common_all_three), "genes\n")
    
    if (length(common_all_three) > 0) {
      cat("    Common genes:", paste(head(common_all_three, 10), collapse = ", "))
      if (length(common_all_three) > 10) cat(" ... (", length(common_all_three), "total)")
      cat("\n")
      
      # Save to CSV
      output_file <- file.path(script_dir, "imm_chronic_common_all_three_genes.csv")
      write_csv(data.frame(Gene_name = common_all_three, stringsAsFactors = FALSE), output_file)
      cat("    Saved to:", basename(output_file), "\n")
      
      assign("imm_common_all_three_chronic", common_all_three, envir = .GlobalEnv)
    }
    
    # Create Venn diagram
    cat("\n  Creating Venn diagram...\n")
    venn_file <- file.path(script_dir, "imm_chronic_venn_diagram.png")
    create_venn_diagram(imm_vars_chronic, 
                       title = "Imm Cell Types - Chronic (Significant Genes)",
                       output_file = venn_file,
                       use_sig_only = TRUE)
  }
  cat("\n")
} else {
  cat("Warning: Need at least 2 Imm chronic variables for intersection analysis\n\n")
}

# ===== INTERSECTION ANALYSIS: Olig Acute (OPC, Premyel, Myel) =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION ANALYSIS: Olig Acute (OPC, Premyelinating, Myelinating)\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")

olig_vars_acute <- c("OPC_acute", "Premyelinating_Olig_acute", "Myelinating_Olig_acute")
olig_sig_genes_acute <- list()

for (var_name in olig_vars_acute) {
  if (exists(var_name, envir = .GlobalEnv)) {
    df <- get(var_name, envir = .GlobalEnv)
    genes <- extract_significant_genes(df)
    olig_sig_genes_acute[[var_name]] <- genes
    cat(var_name, ":", length(genes), "significant genes\n")
  }
}

if (length(olig_sig_genes_acute) >= 2) {
  var_names <- names(olig_sig_genes_acute)
  n_vars <- length(var_names)
  
  # Pairwise intersections
  for (i in 1:(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      var1 <- var_names[i]
      var2 <- var_names[j]
      common <- intersect(olig_sig_genes_acute[[var1]], olig_sig_genes_acute[[var2]])
      cat("  ", var1, "∩", var2, ":", length(common), "genes\n")
      
      if (length(common) > 0) {
        cat("    Common genes:", paste(head(common, 10), collapse = ", "))
        if (length(common) > 10) cat(" ... (", length(common), "total)")
        cat("\n")
        
        # Save to CSV
        pair_name <- paste0(tolower(gsub("_acute|Olig_", "", var1)), "_", tolower(gsub("_acute|Olig_", "", var2)))
        output_file <- file.path(script_dir, paste0("olig_acute_", toupper(substr(pair_name, 1, 1)), substr(pair_name, 2, nchar(pair_name)), "_intersection.csv"))
        write_csv(data.frame(Gene_name = common, stringsAsFactors = FALSE), output_file)
        cat("    Saved to:", basename(output_file), "\n")
      }
    }
  }
  
  # Three-way intersection
  if (n_vars == 3) {
    common_all_three <- Reduce(intersect, olig_sig_genes_acute)
    cat("  Common to all three:", length(common_all_three), "genes\n")
    
    if (length(common_all_three) > 0) {
      cat("    Common genes:", paste(head(common_all_three, 10), collapse = ", "))
      if (length(common_all_three) > 10) cat(" ... (", length(common_all_three), "total)")
      cat("\n")
      
      output_file <- file.path(script_dir, "olig_acute_common_all_three_genes.csv")
      write_csv(data.frame(Gene_name = common_all_three, stringsAsFactors = FALSE), output_file)
      cat("    Saved to:", basename(output_file), "\n")
    }
    
    # Create Venn diagram
    cat("\n  Creating Venn diagram...\n")
    venn_file <- file.path(script_dir, "olig_acute_venn_diagram.png")
    create_venn_diagram(olig_vars_acute, 
                       title = "Olig Cell Types - Acute (Significant Genes)",
                       output_file = venn_file,
                       use_sig_only = TRUE)
  }
  cat("\n")
}

# ===== INTERSECTION ANALYSIS: Olig Chronic (OPC, Premyel, Myel) =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION ANALYSIS: Olig Chronic (OPC, Premyelinating, Myelinating)\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")

olig_vars_chronic <- c("OPC_chronic", "Premyelinating_Olig_chronic", "Myelinating_Olig_chronic")
olig_sig_genes_chronic <- list()

for (var_name in olig_vars_chronic) {
  if (exists(var_name, envir = .GlobalEnv)) {
    df <- get(var_name, envir = .GlobalEnv)
    genes <- extract_significant_genes(df)
    olig_sig_genes_chronic[[var_name]] <- genes
    cat(var_name, ":", length(genes), "significant genes\n")
  }
}

if (length(olig_sig_genes_chronic) >= 2) {
  var_names <- names(olig_sig_genes_chronic)
  n_vars <- length(var_names)
  
  # Pairwise intersections
  for (i in 1:(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      var1 <- var_names[i]
      var2 <- var_names[j]
      common <- intersect(olig_sig_genes_chronic[[var1]], olig_sig_genes_chronic[[var2]])
      cat("  ", var1, "∩", var2, ":", length(common), "genes\n")
      
      if (length(common) > 0) {
        cat("    Common genes:", paste(head(common, 10), collapse = ", "))
        if (length(common) > 10) cat(" ... (", length(common), "total)")
        cat("\n")
        
        pair_name <- paste0(tolower(gsub("_chronic|Olig_", "", var1)), "_", tolower(gsub("_chronic|Olig_", "", var2)))
        output_file <- file.path(script_dir, paste0("olig_chronic_", toupper(substr(pair_name, 1, 1)), substr(pair_name, 2, nchar(pair_name)), "_intersection.csv"))
        write_csv(data.frame(Gene_name = common, stringsAsFactors = FALSE), output_file)
        cat("    Saved to:", basename(output_file), "\n")
      }
    }
  }
  
  # Three-way intersection
  if (n_vars == 3) {
    common_all_three <- Reduce(intersect, olig_sig_genes_chronic)
    cat("  Common to all three:", length(common_all_three), "genes\n")
    
    if (length(common_all_three) > 0) {
      cat("    Common genes:", paste(head(common_all_three, 10), collapse = ", "))
      if (length(common_all_three) > 10) cat(" ... (", length(common_all_three), "total)")
      cat("\n")
      
      output_file <- file.path(script_dir, "olig_chronic_common_all_three_genes.csv")
      write_csv(data.frame(Gene_name = common_all_three, stringsAsFactors = FALSE), output_file)
      cat("    Saved to:", basename(output_file), "\n")
    }
    
    # Create Venn diagram
    cat("\n  Creating Venn diagram...\n")
    venn_file <- file.path(script_dir, "olig_chronic_venn_diagram.png")
    create_venn_diagram(olig_vars_chronic, 
                       title = "Olig Cell Types - Chronic (Significant Genes)",
                       output_file = venn_file,
                       use_sig_only = TRUE)
  }
  cat("\n")
}

# ===== INTERSECTION ANALYSIS: Glia (Astroglia and Microglia) - Acute =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION ANALYSIS: Glia - Astroglia and Microglia (Acute)\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")

glia_vars_acute <- c("Astroglia_acute", "Microglia_acute")
glia_sig_genes_acute <- list()

for (var_name in glia_vars_acute) {
  if (exists(var_name, envir = .GlobalEnv)) {
    df <- get(var_name, envir = .GlobalEnv)
    genes <- extract_significant_genes(df)
    glia_sig_genes_acute[[var_name]] <- genes
    cat(var_name, ":", length(genes), "significant genes\n")
  }
}

if (length(glia_sig_genes_acute) == 2) {
  var_names <- names(glia_sig_genes_acute)
  var1 <- var_names[1]
  var2 <- var_names[2]
  
  genes1 <- glia_sig_genes_acute[[var1]]
  genes2 <- glia_sig_genes_acute[[var2]]
  
  common <- intersect(genes1, genes2)
  only_var1 <- setdiff(genes1, genes2)
  only_var2 <- setdiff(genes2, genes1)
  
  cat("  ", var1, "∩", var2, ":", length(common), "genes\n")
  
  if (length(common) > 0) {
    cat("    Common genes:", paste(head(common, 10), collapse = ", "))
    if (length(common) > 10) cat(" ... (", length(common), "total)")
    cat("\n")
    
    output_file <- file.path(script_dir, "glia_acute_common_both_genes.csv")
    write_csv(data.frame(Gene_name = common, stringsAsFactors = FALSE), output_file)
    cat("    Saved to:", basename(output_file), "\n")
    
    assign("glia_common_both_acute", common, envir = .GlobalEnv)
  }
  
  if (length(only_var1) > 0) {
    output_file <- file.path(script_dir, paste0("glia_acute_only_", gsub("_acute", "", var1), ".csv"))
    write_csv(data.frame(Gene_name = only_var1, stringsAsFactors = FALSE), output_file)
    cat("    Only", var1, ":", length(only_var1), "genes\n")
  }
  
  if (length(only_var2) > 0) {
    output_file <- file.path(script_dir, paste0("glia_acute_only_", gsub("_acute", "", var2), ".csv"))
    write_csv(data.frame(Gene_name = only_var2, stringsAsFactors = FALSE), output_file)
    cat("    Only", var2, ":", length(only_var2), "genes\n")
  }
  
  # Create Venn diagram
  cat("\n  Creating Venn diagram...\n")
  venn_file <- file.path(script_dir, "glia_acute_venn_diagram.png")
  create_venn_diagram(glia_vars_acute, 
                     title = "Glia Cell Types - Acute (Astroglia and Microglia)",
                     output_file = venn_file,
                     use_sig_only = TRUE)
  cat("\n")
} else {
  cat("Warning: Need both Astroglia_acute and Microglia_acute variables for intersection analysis\n\n")
}

# ===== INTERSECTION ANALYSIS: Glia (Astroglia and Microglia) - Chronic =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION ANALYSIS: Glia - Astroglia and Microglia (Chronic)\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")

glia_vars_chronic <- c("Astroglia_chronic", "Microglia_chronic")
glia_sig_genes_chronic <- list()

for (var_name in glia_vars_chronic) {
  if (exists(var_name, envir = .GlobalEnv)) {
    df <- get(var_name, envir = .GlobalEnv)
    genes <- extract_significant_genes(df)
    glia_sig_genes_chronic[[var_name]] <- genes
    cat(var_name, ":", length(genes), "significant genes\n")
  }
}

if (length(glia_sig_genes_chronic) == 2) {
  var_names <- names(glia_sig_genes_chronic)
  var1 <- var_names[1]
  var2 <- var_names[2]
  
  genes1 <- glia_sig_genes_chronic[[var1]]
  genes2 <- glia_sig_genes_chronic[[var2]]
  
  common <- intersect(genes1, genes2)
  only_var1 <- setdiff(genes1, genes2)
  only_var2 <- setdiff(genes2, genes1)
  
  cat("  ", var1, "∩", var2, ":", length(common), "genes\n")
  
  if (length(common) > 0) {
    cat("    Common genes:", paste(head(common, 10), collapse = ", "))
    if (length(common) > 10) cat(" ... (", length(common), "total)")
    cat("\n")
    
    output_file <- file.path(script_dir, "glia_chronic_common_both_genes.csv")
    write_csv(data.frame(Gene_name = common, stringsAsFactors = FALSE), output_file)
    cat("    Saved to:", basename(output_file), "\n")
    
    assign("glia_common_both_chronic", common, envir = .GlobalEnv)
  }
  
  if (length(only_var1) > 0) {
    output_file <- file.path(script_dir, paste0("glia_chronic_only_", gsub("_chronic", "", var1), ".csv"))
    write_csv(data.frame(Gene_name = only_var1, stringsAsFactors = FALSE), output_file)
    cat("    Only", var1, ":", length(only_var1), "genes\n")
  }
  
  if (length(only_var2) > 0) {
    output_file <- file.path(script_dir, paste0("glia_chronic_only_", gsub("_chronic", "", var2), ".csv"))
    write_csv(data.frame(Gene_name = only_var2, stringsAsFactors = FALSE), output_file)
    cat("    Only", var2, ":", length(only_var2), "genes\n")
  }
  
  # Create Venn diagram
  cat("\n  Creating Venn diagram...\n")
  venn_file <- file.path(script_dir, "glia_chronic_venn_diagram.png")
  create_venn_diagram(glia_vars_chronic, 
                     title = "Glia Cell Types - Chronic (Astroglia and Microglia)",
                     output_file = venn_file,
                     use_sig_only = TRUE)
  cat("\n")
} else {
  cat("Warning: Need both Astroglia_chronic and Microglia_chronic variables for intersection analysis\n\n")
}

# ===== INTERSECTION ANALYSIS: Excitatory with Imm (Acute) =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION ANALYSIS: Excitatory with Imm Cell Types (Acute)\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")

if (exists("excitatory_acute", envir = .GlobalEnv)) {
  df_excitatory_acute <- get("excitatory_acute", envir = .GlobalEnv)
  genes_excitatory_acute <- extract_significant_genes(df_excitatory_acute)
  cat("Excitatory_acute:", length(genes_excitatory_acute), "significant genes\n\n")
  
  imm_types_acute <- c("Imm_CGE_acute", "Imm_LGE_acute", "Imm_MGE_acute")
  excitatory_imm_intersections_acute <- list()
  
  for (imm_var in imm_types_acute) {
    if (exists(imm_var, envir = .GlobalEnv)) {
      df_imm <- get(imm_var, envir = .GlobalEnv)
      genes_imm <- extract_significant_genes(df_imm)
      common <- intersect(genes_excitatory_acute, genes_imm)
      
      imm_type_short <- gsub("Imm_|_acute", "", imm_var)
      cat("Excitatory_acute ∩", imm_var, ":", length(common), "genes\n")
      
      if (length(common) > 0) {
        cat("  Common genes:", paste(head(common, 10), collapse = ", "))
        if (length(common) > 10) cat(" ... (", length(common), "total)")
        cat("\n")
        
        output_file <- file.path(script_dir, paste0("excitatory_", tolower(imm_type_short), "_intersection_acute.csv"))
        write_csv(data.frame(Gene_name = common, stringsAsFactors = FALSE), output_file)
        cat("  Saved to:", basename(output_file), "\n")
        
        excitatory_imm_intersections_acute[[imm_type_short]] <- common
      }
      cat("\n")
    }
  }
  
  # Create Venn diagram for excitatory with all Imm types
  if (length(excitatory_imm_intersections_acute) >= 2) {
    cat("Creating Venn diagram for Excitatory_acute with Imm types...\n")
    venn_lists <- list(Excitatory = genes_excitatory_acute)
    for (name in names(excitatory_imm_intersections_acute)) {
      venn_lists[[paste0("Imm_", toupper(name))]] <- excitatory_imm_intersections_acute[[name]]
    }
    
    if (length(venn_lists) <= 3) {
      venn_file <- file.path(script_dir, "excitatory_imm_acute_venn_diagram.png")
      create_venn_diagram_from_lists(venn_lists, 
                                     title = "Excitatory with Imm Cell Types - Acute",
                                     output_file = venn_file)
    }
  }
} else {
  cat("Warning: excitatory_acute variable not found\n\n")
}

# ===== INTERSECTION ANALYSIS: Excitatory with Imm (Chronic) =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION ANALYSIS: Excitatory with Imm Cell Types (Chronic)\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")

if (exists("excitatory_chronic", envir = .GlobalEnv)) {
  df_excitatory_chronic <- get("excitatory_chronic", envir = .GlobalEnv)
  genes_excitatory_chronic <- extract_significant_genes(df_excitatory_chronic)
  cat("Excitatory_chronic:", length(genes_excitatory_chronic), "significant genes\n\n")
  
  imm_types_chronic <- c("Imm_CGE_chronic", "Imm_LGE_chronic", "Imm_MGE_chronic")
  excitatory_imm_intersections_chronic <- list()
  
  for (imm_var in imm_types_chronic) {
    if (exists(imm_var, envir = .GlobalEnv)) {
      df_imm <- get(imm_var, envir = .GlobalEnv)
      genes_imm <- extract_significant_genes(df_imm)
      common <- intersect(genes_excitatory_chronic, genes_imm)
      
      imm_type_short <- gsub("Imm_|_chronic", "", imm_var)
      cat("Excitatory_chronic ∩", imm_var, ":", length(common), "genes\n")
      
      if (length(common) > 0) {
        cat("  Common genes:", paste(head(common, 10), collapse = ", "))
        if (length(common) > 10) cat(" ... (", length(common), "total)")
        cat("\n")
        
        output_file <- file.path(script_dir, paste0("excitatory_", tolower(imm_type_short), "_intersection_chronic.csv"))
        write_csv(data.frame(Gene_name = common, stringsAsFactors = FALSE), output_file)
        cat("  Saved to:", basename(output_file), "\n")
        
        excitatory_imm_intersections_chronic[[imm_type_short]] <- common
      }
      cat("\n")
    }
  }
  
  # Create Venn diagram
  if (length(excitatory_imm_intersections_chronic) >= 2) {
    cat("Creating Venn diagram for Excitatory_chronic with Imm types...\n")
    venn_lists <- list(Excitatory = genes_excitatory_chronic)
    for (name in names(excitatory_imm_intersections_chronic)) {
      venn_lists[[paste0("Imm_", toupper(name))]] <- excitatory_imm_intersections_chronic[[name]]
    }
    
    if (length(venn_lists) <= 3) {
      venn_file <- file.path(script_dir, "excitatory_imm_chronic_venn_diagram.png")
      create_venn_diagram_from_lists(venn_lists, 
                                     title = "Excitatory with Imm Cell Types - Chronic",
                                     output_file = venn_file)
    }
  }
} else {
  cat("Warning: excitatory_chronic variable not found\n\n")
}

# ===== INTERSECTION ANALYSIS: Excitatory vs All Inhibitory Neurons (Acute) =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION ANALYSIS: Excitatory vs All Inhibitory Neurons (Acute)\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")

if (exists("excitatory_acute", envir = .GlobalEnv)) {
  df_excitatory_acute <- get("excitatory_acute", envir = .GlobalEnv)
  genes_excitatory_acute <- extract_significant_genes(df_excitatory_acute)
  cat("Excitatory_acute:", length(genes_excitatory_acute), "significant genes\n")
  
  # Collect all inhibitory neuron genes (CGE, LGE, MGE combined)
  imm_types_acute <- c("Imm_CGE_acute", "Imm_LGE_acute", "Imm_MGE_acute")
  all_inhibitory_genes_acute <- c()
  
  for (imm_var in imm_types_acute) {
    if (exists(imm_var, envir = .GlobalEnv)) {
      df_imm <- get(imm_var, envir = .GlobalEnv)
      genes_imm <- extract_significant_genes(df_imm)
      all_inhibitory_genes_acute <- c(all_inhibitory_genes_acute, genes_imm)
      cat(imm_var, ":", length(genes_imm), "significant genes\n")
    }
  }
  
  # Get unique inhibitory genes (union of all inhibitory types)
  all_inhibitory_genes_acute <- unique(all_inhibitory_genes_acute)
  cat("All inhibitory neurons (CGE+LGE+MGE union):", length(all_inhibitory_genes_acute), "unique genes\n\n")
  
  # Find intersection between excitatory and all inhibitory
  common_excitatory_inhibitory_acute <- intersect(genes_excitatory_acute, all_inhibitory_genes_acute)
  only_excitatory_acute <- setdiff(genes_excitatory_acute, all_inhibitory_genes_acute)
  only_inhibitory_acute <- setdiff(all_inhibitory_genes_acute, genes_excitatory_acute)
  
  cat("Excitatory_acute ∩ All_Inhibitory_acute:", length(common_excitatory_inhibitory_acute), "genes\n")
  cat("Only in Excitatory_acute:", length(only_excitatory_acute), "genes\n")
  cat("Only in All_Inhibitory_acute:", length(only_inhibitory_acute), "genes\n")
  
  if (length(common_excitatory_inhibitory_acute) > 0) {
    cat("\nCommon genes:", paste(head(common_excitatory_inhibitory_acute, 20), collapse = ", "))
    if (length(common_excitatory_inhibitory_acute) > 20) cat(" ... (", length(common_excitatory_inhibitory_acute), "total)")
    cat("\n")
    
    # Save common genes
    output_file <- file.path(script_dir, "excitatory_inhibitory_common_acute.csv")
    write_csv(data.frame(Gene_name = common_excitatory_inhibitory_acute, stringsAsFactors = FALSE), output_file)
    cat("Saved common genes to:", basename(output_file), "\n")
    
    # Create Venn diagram
    cat("\nCreating Venn diagram for Excitatory vs All Inhibitory (Acute)...\n")
    venn_lists <- list(
      Excitatory = genes_excitatory_acute,
      All_Inhibitory = all_inhibitory_genes_acute
    )
    venn_file <- file.path(script_dir, "excitatory_vs_inhibitory_acute_venn_diagram.png")
    create_venn_diagram_from_lists(venn_lists, 
                                   title = "Excitatory vs All Inhibitory Neurons - Acute",
                                   output_file = venn_file)
    
    # Run enrichment analysis ONLY on common genes
    cat("\nRunning enrichment analysis on common genes (Excitatory ∩ All Inhibitory - Acute)...\n")
    output_name <- "excitatory_inhibitory_common_acute"
    
    if (enrichr_available) {
      run_enrichr_analysis(common_excitatory_inhibitory_acute, output_name, results_base_dir = script_dir)
    }
    
    if (clusterprofiler_available && orgdb_available) {
      run_clusterprofiler_analysis(common_excitatory_inhibitory_acute, output_name, results_base_dir = script_dir)
    }
  } else {
    cat("\nNo common genes found between Excitatory and All Inhibitory neurons (Acute)\n")
  }
  cat("\n")
} else {
  cat("Warning: excitatory_acute variable not found\n\n")
}

# ===== INTERSECTION ANALYSIS: Excitatory vs All Inhibitory Neurons (Chronic) =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION ANALYSIS: Excitatory vs All Inhibitory Neurons (Chronic)\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")

if (exists("excitatory_chronic", envir = .GlobalEnv)) {
  df_excitatory_chronic <- get("excitatory_chronic", envir = .GlobalEnv)
  genes_excitatory_chronic <- extract_significant_genes(df_excitatory_chronic)
  cat("Excitatory_chronic:", length(genes_excitatory_chronic), "significant genes\n")
  
  # Collect all inhibitory neuron genes (CGE, LGE, MGE combined)
  imm_types_chronic <- c("Imm_CGE_chronic", "Imm_LGE_chronic", "Imm_MGE_chronic")
  all_inhibitory_genes_chronic <- c()
  
  for (imm_var in imm_types_chronic) {
    if (exists(imm_var, envir = .GlobalEnv)) {
      df_imm <- get(imm_var, envir = .GlobalEnv)
      genes_imm <- extract_significant_genes(df_imm)
      all_inhibitory_genes_chronic <- c(all_inhibitory_genes_chronic, genes_imm)
      cat(imm_var, ":", length(genes_imm), "significant genes\n")
    }
  }
  
  # Get unique inhibitory genes (union of all inhibitory types)
  all_inhibitory_genes_chronic <- unique(all_inhibitory_genes_chronic)
  cat("All inhibitory neurons (CGE+LGE+MGE union):", length(all_inhibitory_genes_chronic), "unique genes\n\n")
  
  # Find intersection between excitatory and all inhibitory
  common_excitatory_inhibitory_chronic <- intersect(genes_excitatory_chronic, all_inhibitory_genes_chronic)
  only_excitatory_chronic <- setdiff(genes_excitatory_chronic, all_inhibitory_genes_chronic)
  only_inhibitory_chronic <- setdiff(all_inhibitory_genes_chronic, genes_excitatory_chronic)
  
  cat("Excitatory_chronic ∩ All_Inhibitory_chronic:", length(common_excitatory_inhibitory_chronic), "genes\n")
  cat("Only in Excitatory_chronic:", length(only_excitatory_chronic), "genes\n")
  cat("Only in All_Inhibitory_chronic:", length(only_inhibitory_chronic), "genes\n")
  
  if (length(common_excitatory_inhibitory_chronic) > 0) {
    cat("\nCommon genes:", paste(head(common_excitatory_inhibitory_chronic, 20), collapse = ", "))
    if (length(common_excitatory_inhibitory_chronic) > 20) cat(" ... (", length(common_excitatory_inhibitory_chronic), "total)")
    cat("\n")
    
    # Save common genes
    output_file <- file.path(script_dir, "excitatory_inhibitory_common_chronic.csv")
    write_csv(data.frame(Gene_name = common_excitatory_inhibitory_chronic, stringsAsFactors = FALSE), output_file)
    cat("Saved common genes to:", basename(output_file), "\n")
    
    # Create Venn diagram
    cat("\nCreating Venn diagram for Excitatory vs All Inhibitory (Chronic)...\n")
    venn_lists <- list(
      Excitatory = genes_excitatory_chronic,
      All_Inhibitory = all_inhibitory_genes_chronic
    )
    venn_file <- file.path(script_dir, "excitatory_vs_inhibitory_chronic_venn_diagram.png")
    create_venn_diagram_from_lists(venn_lists, 
                                   title = "Excitatory vs All Inhibitory Neurons - Chronic",
                                   output_file = venn_file)
    
    # Run enrichment analysis ONLY on common genes
    cat("\nRunning enrichment analysis on common genes (Excitatory ∩ All Inhibitory - Chronic)...\n")
    output_name <- "excitatory_inhibitory_common_chronic"
    
    if (enrichr_available) {
      run_enrichr_analysis(common_excitatory_inhibitory_chronic, output_name, results_base_dir = script_dir)
    }
    
    if (clusterprofiler_available && orgdb_available) {
      run_clusterprofiler_analysis(common_excitatory_inhibitory_chronic, output_name, results_base_dir = script_dir)
    }
  } else {
    cat("\nNo common genes found between Excitatory and All Inhibitory neurons (Chronic)\n")
  }
  cat("\n")
} else {
  cat("Warning: excitatory_chronic variable not found\n\n")
}

# ===== ACUTE vs CHRONIC COMPARISONS =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("ACUTE vs CHRONIC COMPARISONS (by cell type)\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")

# Get base name function
get_base_name <- function(filename) {
  base <- gsub("^cell_type_ISP_", "", filename)
  base <- gsub("_acute_GPU_ONLY_genes\\.csv$", "", base)
  base <- gsub("_chronic_GPU_ONLY_genes\\.csv$", "", base)
  return(base)
}

# Compare acute vs chronic for each cell type
for (var_name in created_variables) {
  if (grepl("_acute$", var_name)) {
    base_name <- gsub("_acute$", "", var_name)
    chronic_var <- paste0(base_name, "_chronic")
    
    if (exists(chronic_var, envir = .GlobalEnv)) {
      df_acute <- get(var_name, envir = .GlobalEnv)
      df_chronic <- get(chronic_var, envir = .GlobalEnv)
      
      genes_acute <- extract_significant_genes(df_acute)
      genes_chronic <- extract_significant_genes(df_chronic)
      
      common <- intersect(genes_acute, genes_chronic)
      acute_only <- setdiff(genes_acute, genes_chronic)
      chronic_only <- setdiff(genes_chronic, genes_acute)
      
      cat(base_name, ":\n")
      cat("  Acute:", length(genes_acute), "genes\n")
      cat("  Chronic:", length(genes_chronic), "genes\n")
      cat("  Common:", length(common), "genes\n")
      cat("  Acute only:", length(acute_only), "genes\n")
      cat("  Chronic only:", length(chronic_only), "genes\n")
      
      if (length(common) > 0) {
        cat("  Common genes:", paste(head(common, 10), collapse = ", "))
        if (length(common) > 10) cat(" ... (", length(common), "total)")
        cat("\n")
        
        output_file <- file.path(script_dir, paste0(base_name, "_common_acute_chronic.csv"))
        write_csv(data.frame(Gene_name = common, stringsAsFactors = FALSE), output_file)
        cat("  Saved to:", basename(output_file), "\n")
      }
      
      # Create Venn diagram for acute vs chronic
      if (length(genes_acute) > 0 && length(genes_chronic) > 0) {
        venn_lists <- list(Acute = genes_acute, Chronic = genes_chronic)
        venn_file <- file.path(script_dir, paste0(base_name, "_acute_vs_chronic_venn_diagram.png"))
        create_venn_diagram_from_lists(venn_lists, 
                                      title = paste0(base_name, " - Acute vs Chronic"),
                                      output_file = venn_file)
      }
      cat("\n")
    }
  }
}

# ===== RUN ENRICHMENT ANALYSES ON ALL INTERSECTIONS =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("RUNNING ENRICHMENT ANALYSES ON ALL INTERSECTIONS\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")

# Collect all intersection gene lists
intersection_files <- list.files(path = script_dir, pattern = "*intersection.*\\.csv$|*common.*\\.csv$", full.names = FALSE)
intersection_files <- intersection_files[!grepl("_enrichR_|_clusterProfiler_", intersection_files)]
# Skip excitatory_inhibitory_common files since they were already processed with enrichment
intersection_files <- intersection_files[!grepl("excitatory_inhibitory_common", intersection_files)]

cat("Found", length(intersection_files), "intersection files for enrichment analysis\n")
cat("(Note: excitatory_inhibitory_common files already processed with enrichment above)\n\n")

for (intersection_file in intersection_files) {
  file_path <- file.path(script_dir, intersection_file)
  
  if (file.exists(file_path)) {
    tryCatch({
      gene_data <- read_csv(file_path, col_types = cols(), show_col_types = FALSE)
      
      if ("Gene_name" %in% colnames(gene_data)) {
        gene_list <- gene_data$Gene_name
      } else if (ncol(gene_data) > 0) {
        gene_list <- gene_data[[1]]
      } else {
        next
      }
      
      gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
      
      if (length(gene_list) > 0) {
        output_name <- gsub("\\.csv$", "", intersection_file)
        cat("Processing:", output_name, "(", length(gene_list), "genes)\n")
        
        # Run enrichR
        if (enrichr_available) {
          run_enrichr_analysis(gene_list, output_name, results_base_dir = script_dir)
        }
        
        # Run clusterProfiler
        if (clusterprofiler_available && orgdb_available) {
          run_clusterprofiler_analysis(gene_list, output_name, results_base_dir = script_dir)
        }
        
        cat("\n")
      }
    }, error = function(e) {
      cat("Error processing", intersection_file, ":", conditionMessage(e), "\n")
    })
  }
}

# ===== SUMMARY =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("ANALYSIS COMPLETE\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")
cat("Results saved to:", script_dir, "\n")
cat("Files processed:", length(created_variables), "\n")
cat("Venn diagrams created for all intersections\n")
cat("Common genes printed and saved to CSV files\n")
cat("Enrichment analyses completed\n")
cat("\n")
