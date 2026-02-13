# ============================================================================
# SCRIPT: INTERSECT ACUTE CELL TYPES AND RUN ENRICHMENT
# ============================================================================
# This script:
# 1. Reads specified CSV files with acute cell type data
# 2. Extracts significant hits (Sig == 1) from each file
# 3. Computes intersection of all significant hits across all files
# 4. Runs enrichR and clusterProfiler enrichment analyses on common genes
# 5. Creates Venn diagram showing the intersection
# 6. Saves results to CSV files
# ============================================================================
# UPDATED: File list modified to match actual CSV files in directory
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

# ===== FUNCTION TO EXTRACT CELL TYPE NAME FROM FILENAME =====
extract_cell_type_name <- function(filename) {
  # Extract the part between "ISP_" and "_2dHIE"
  # Pattern: cell_type_ISP_<CELL_TYPE>_2dHIE_to_non_hie_jan2026_genes.csv
  if (grepl("ISP_.*_2dHIE", filename)) {
    cell_type <- gsub("^.*ISP_", "", filename)
    cell_type <- gsub("_2dHIE.*$", "", cell_type)
    return(cell_type)
  } else {
    # Fallback: return filename without extension
    return(gsub("\\.csv$", "", basename(filename)))
  }
}

# ===== HELPER FUNCTION TO FILTER SIGNIFICANT GENES =====
filter_significant_genes <- function(df) {
  if (!"Sig" %in% colnames(df)) {
    return(df %>% filter(FALSE))  # Return empty if Sig column doesn't exist
  }
  
  # Filter with only Sig == 1
  df_filtered <- df %>%
    filter(Sig == 1)
  
  return(df_filtered)
}

# ===== FUNCTION TO FILTER AND SAVE HITS WITH SIG=1 AND POSITIVE SHIFT =====
filter_and_save_hits <- function(df, original_filename, output_dir) {
  # Check required columns exist
  if (!"Sig" %in% colnames(df) || !"Shift_to_goal_end" %in% colnames(df)) {
    cat("    Warning: Required columns (Sig or Shift_to_goal_end) not found\n")
    return(NULL)
  }
  
  # Filter for Sig == 1 AND Shift_to_goal_end > 0
  hits <- df %>%
    filter(Sig == 1, Shift_to_goal_end > 0)
  
  if (nrow(hits) == 0) {
    cat("    No hits found (Sig=1 and Shift_to_goal_end > 0)\n")
    return(NULL)
  }
  
  # Create output filename - use original filename without .csv extension
  base_filename <- gsub("\\.csv$", "", basename(original_filename))
  output_file <- file.path(output_dir, paste0(base_filename, "_signif_hits_positive_shift.txt"))
  
  # Write to text file
  write.table(hits, 
              file = output_file,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)
  
  cat("    Saved", nrow(hits), "hits to:", basename(output_file), "\n")
  return(hits)
}

# ===== FUNCTION TO EXTRACT SIGNIFICANT GENES FROM DATAFRAME =====
extract_significant_genes <- function(df) {
  # First filter for significant genes
  df_sig <- filter_significant_genes(df)
  
  # Then extract gene names
  if ("Gene_name" %in% colnames(df_sig)) {
    genes <- unique(df_sig$Gene_name[!is.na(df_sig$Gene_name) & df_sig$Gene_name != ""])
    return(genes)
  } else {
    return(character(0))
  }
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
  # enrichr_databases <- c(
  #  "GO_Biological_Process_2021", "GO_Molecular_Function_2021", 
  #  "GO_Cellular_Component_2021", "KEGG_2021_Human", 
  #  "Reactome_2022", "WikiPathways_2021_Human"
  # )

  enrichr_databases <- c(
  # Updated Gene Ontology (2025 versions)
  "GO_Biological_Process_2025",
  "GO_Molecular_Function_2025", 
  "GO_Cellular_Component_2025",
  
  # Core pathway databases
  "KEGG_2021_Human",
  
  # Pathway Databases (from your images)
  "Reactome_Pathways_2024",
  "WikiPathways_2024_Human", 
  "BioPlanet_2019",
  "WikiPathways_2024_Mouse",
  
  # Specialized Collections (from your images)
  "Elsevier_Pathway_Collection",
  "MSigDB_Hallmark_2020",
  "BioCarta_2016",
  "HumanCyc_2016",
  "NCI-Nature_2016",
  "Panther_2016",
  
  # Coexpression Database (from your images)
  # "ARCHS4_Kinases_Coexp",
  
  # Phenotype and Disease Databases
  # "MGI_Mammalian_Phenotype_Level_4_2024",
  "Jensen_DISEASES_Curated_2025",
  "Jensen_DISEASES_Experimental_2025",
  "Human_Phenotype_Ontology",
  "SynGO_2024",
  # "KOMP2_Mouse_Phenotypes_2022",
  
  # Additional commonly used databases
  "ENCODE_TF_ChIP-seq_2015",
  "ChEA_2016",
  
  # Additional databases that may be available with different names
  "Reactome_2022",
  "WikiPathways_2023_Human"
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

# ===== FUNCTION TO CREATE VENN DIAGRAM FROM GENE LISTS =====
create_venn_diagram_from_lists <- function(gene_lists, title = "Venn Diagram", output_file = NULL) {
  cat("  Creating Venn diagram with", length(gene_lists), "gene lists...\n")
  
  # Check if ggvenn is available
  if (!ggvenn_available) {
    cat("Warning: ggvenn package not available. Skipping Venn diagram.\n")
    return(NULL)
  }
  
  if (length(gene_lists) < 2) {
    cat("Warning: Need at least 2 gene lists to create Venn diagram.\n")
    return(NULL)
  }
  
  # Check if any gene lists are empty
  empty_lists <- sapply(gene_lists, function(x) length(x) == 0)
  if (all(empty_lists)) {
    cat("Warning: All gene lists are empty. Cannot create Venn diagram.\n")
    return(NULL)
  }
  
  # For 2-3 sets, use ggvenn
  if (length(gene_lists) <= 3) {
    tryCatch({
      cat("  Using ggvenn for", length(gene_lists), "sets...\n")
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
        cat("  Saving Venn diagram to:", output_file, "\n")
        ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
        cat("  Venn diagram saved successfully to:", basename(output_file), "\n")
      } else {
        cat("  Warning: No output file specified for Venn diagram.\n")
      }
      
      return(p)
    }, error = function(e) {
      cat("  Error creating Venn diagram:", conditionMessage(e), "\n")
      cat("  Error details:", toString(e), "\n")
      return(NULL)
    })
  } else {
    # For more than 3 sets, use UpSet plot
    cat("  Note: More than 3 sets detected. Using UpSet plot instead.\n")
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
          cat("  Saving UpSet plot to:", output_file, "\n")
          png(output_file, width = 1200, height = 800)
          print(upset(binary_matrix, sets = names(gene_lists), order.by = "freq"))
          dev.off()
          cat("  UpSet plot saved successfully to:", basename(output_file), "\n")
        } else {
          cat("  Warning: No output file specified for UpSet plot.\n")
        }
        return(TRUE)
      }, error = function(e) {
        cat("  Error creating UpSet plot:", conditionMessage(e), "\n")
        return(NULL)
      })
    } else {
      cat("  Warning: UpSetR package not available. Cannot create plot for more than 3 sets.\n")
      return(NULL)
    }
  }
}

# ===== MAIN ANALYSIS =====
cat("=", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECT CELL TYPES AND RUN ENRICHMENT\n")
cat("=", paste(rep("=", 60), collapse=""), "\n\n")

# Define CSV files to process - UPDATED to match actual files in directory
csv_files <- c(
  "cell_type_ISP_Ast_Rgl_2dHIE_to_non_hie_jan2026_genes.csv",
  "cell_type_ISP_Imm_CGE_2dHIE_to_non_hie_jan2026_genes.csv",
  "cell_type_ISP_Imm_CGE_LGE_2dHIE_to_non_hie_jan2026_genes.csv",
  "cell_type_ISP_Imm_Exc_1_2dHIE_to_non_hie_jan2026_genes.csv",
  "cell_type_ISP_Imm_Exc_2_2dHIE_to_non_hie_jan2026_genes.csv",
  "cell_type_ISP_Imm_LGE_2dHIE_to_non_hie_jan2026_genes.csv",
  "cell_type_ISP_Imm_MGE_2dHIE_to_non_hie_jan2026_genes.csv",
  "cell_type_ISP_Imm_MGE_MAF_CRABP1_2dHIE_to_non_hie_jan2026_genes.csv",
  "cell_type_ISP_Imm_MGE_MAF_NPY_2dHIE_to_non_hie_jan2026_genes.csv",
  "cell_type_ISP_Mcg_2dHIE_to_non_hie_jan2026_genes.csv",
  "cell_type_ISP_Myel_Olig_2dHIE_to_non_hie_jan2026_genes.csv",
  "cell_type_ISP_OPC_2dHIE_to_non_hie_jan2026_genes.csv",
  "cell_type_ISP_Prem_Olig_2dHIE_to_non_hie_jan2026_genes.csv"
)

# Store significant genes from each file (using cell type names as keys)
all_gene_lists <- list()
file_summary <- data.frame(
  File = character(),
  Cell_Type = character(),
  Total_Rows = integer(),
  Significant_Genes = integer(),
  Hits_Sig1_PosShift = integer(),
  stringsAsFactors = FALSE
)

# Read and process each CSV file
cat("Reading and processing CSV files...\n\n")
for (csv_file in csv_files) {
  file_path <- file.path(script_dir, csv_file)
  file_name <- basename(csv_file)
  
  # Extract cell type name from filename
  cell_type_name <- extract_cell_type_name(file_name)
  
  cat("Processing:", file_name, "(Cell Type:", cell_type_name, ")...")
  
  if (!file.exists(file_path)) {
    cat(" ERROR: File not found\n")
    next
  }
  
  tryCatch({
    df <- read_csv(file_path, col_types = cols(), show_col_types = FALSE)
    
    # Filter and save hits with Sig=1 and positive Shift_to_goal_end
    cat("\n  Filtering hits (Sig=1, Shift_to_goal_end > 0)...")
    hits <- filter_and_save_hits(df, file_name, script_dir)
    
    # Extract gene names from hits for enrichment analysis
    hits_genes <- character(0)
    if (!is.null(hits) && nrow(hits) > 0) {
      if ("Gene_name" %in% colnames(hits)) {
        hits_genes <- unique(hits$Gene_name[!is.na(hits$Gene_name) & hits$Gene_name != ""])
      }
    }
    
    # Run enrichment analysis on hits genes if available
    if (length(hits_genes) > 0) {
      cat("\n  Running enrichment analysis on", length(hits_genes), "hits genes...")
      base_filename <- gsub("\\.csv$", "", file_name)
      output_name <- paste0(base_filename, "_signif_hits_positive_shift")
      
      # Run enrichR
      if (enrichr_available) {
        cat("\n    Running enrichR...")
        enrichr_result <- run_enrichr_analysis(hits_genes, output_name, results_base_dir = script_dir)
      }
      
      # Run clusterProfiler
      if (clusterprofiler_available && orgdb_available) {
        cat("\n    Running clusterProfiler...")
        clusterprofiler_result <- run_clusterprofiler_analysis(hits_genes, output_name, results_base_dir = script_dir)
      }
      cat("\n")
    } else {
      cat("\n  No hits genes available for enrichment analysis\n")
    }
    
    # Extract significant genes (for intersection analysis)
    genes <- extract_significant_genes(df)
    
    # Store gene list using cell type name as key
    all_gene_lists[[cell_type_name]] <- genes
    
    # Add to summary
    hits_count <- ifelse(is.null(hits), 0, nrow(hits))
    file_summary <- rbind(file_summary, data.frame(
      File = file_name,
      Cell_Type = cell_type_name,
      Total_Rows = nrow(df),
      Significant_Genes = length(genes),
      Hits_Sig1_PosShift = hits_count,
      stringsAsFactors = FALSE
    ))
    
    cat("\n  OK (", nrow(df), "rows,", length(genes), "significant genes,", hits_count, "hits with positive shift)\n")
    
  }, error = function(e) {
    cat(" ERROR:", conditionMessage(e), "\n")
  })
}

cat("\n")

# Print summary
cat("File Summary:\n")
print(file_summary)
cat("\n")

# ===== FUNCTION TO COMPUTE INTERSECTION AND RUN ANALYSIS =====
compute_intersection_and_analyze <- function(cell_type_names, intersection_name, intersection_number, output_dir) {
  cat("\n", paste(rep("=", 70), collapse=""), "\n")
  cat("INTERSECTION", intersection_number, ":", intersection_name, "\n")
  cat("Cell types:", paste(cell_type_names, collapse = ", "), "\n")
  cat(paste(rep("=", 70), collapse=""), "\n\n")
  
  # Get gene lists for specified cell types
  intersection_lists <- list()
  missing_types <- character(0)
  
  for (ct in cell_type_names) {
    if (ct %in% names(all_gene_lists)) {
      intersection_lists[[ct]] <- all_gene_lists[[ct]]
    } else {
      missing_types <- c(missing_types, ct)
    }
  }
  
  if (length(missing_types) > 0) {
    cat("Warning: Cell types not found:", paste(missing_types, collapse = ", "), "\n")
  }
  
  if (length(intersection_lists) == 0) {
    cat("No valid gene lists found for this intersection. Skipping.\n\n")
    return(NULL)
  }
  
  # Compute intersection
  if (length(intersection_lists) == 1) {
    common_genes <- intersection_lists[[1]]
  } else {
    common_genes <- intersection_lists[[1]]
    for (i in 2:length(intersection_lists)) {
      common_genes <- intersect(common_genes, intersection_lists[[i]])
    }
  }
  
  cat("Common significant genes:", length(common_genes), "\n")
  
  if (length(common_genes) > 0) {
    cat("First 20 common genes:", paste(head(common_genes, 20), collapse = ", "))
    if (length(common_genes) > 20) cat(" ... (", length(common_genes), "total)")
    cat("\n\n")
    
    # Save intersection to CSV in INTERSECTIONS folder
    safe_name <- gsub("[^A-Za-z0-9_]", "_", intersection_name)
    output_file <- file.path(output_dir, paste0("intersection_", intersection_number, "_", safe_name, "_common_genes.csv"))
    write_csv(data.frame(Gene_name = common_genes, stringsAsFactors = FALSE), output_file)
    cat("Saved common genes to:", basename(output_file), "\n")
    
    # Create Venn diagram in INTERSECTIONS folder
    if (length(intersection_lists) >= 2) {
      cat("\nCreating Venn diagram...\n")
      venn_file <- file.path(output_dir, paste0("intersection_", intersection_number, "_", safe_name, "_venn_diagram.png"))
      
      # Format cell type names for display
      venn_gene_lists <- intersection_lists
      display_names <- gsub("_", " ", names(venn_gene_lists))
      names(venn_gene_lists) <- display_names
      
      venn_result <- create_venn_diagram_from_lists(venn_gene_lists, 
                                    title = paste("Intersection", intersection_number, "-", intersection_name),
                                    output_file = venn_file)
      
      if (is.null(venn_result)) {
        cat("Warning: Venn diagram creation failed or was skipped.\n")
      }
    }
    
    # Run enrichment analyses - save to INTERSECTIONS folder
    cat("\nRunning enrichment analyses...\n")
    output_name <- paste0("intersection_", intersection_number, "_", safe_name)
    
    # Run enrichR
    if (enrichr_available) {
      enrichr_result <- run_enrichr_analysis(common_genes, output_name, results_base_dir = output_dir)
    }
    
    # Run clusterProfiler
    if (clusterprofiler_available && orgdb_available) {
      clusterprofiler_result <- run_clusterprofiler_analysis(common_genes, output_name, results_base_dir = output_dir)
    }
    
    cat("\n")
    return(common_genes)
    
  } else {
    cat("\nNo common significant genes found. Skipping enrichment analyses.\n\n")
    return(NULL)
  }
}

# ===== COMPUTE SPECIFIED INTERSECTIONS =====
if (length(all_gene_lists) > 0) {
  
  # Create INTERSECTIONS folder for all intersection results
  intersections_dir <- file.path(script_dir, "INTERSECTIONS")
  dir.create(intersections_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Created INTERSECTIONS folder:", intersections_dir, "\n\n")
  
  # Save summary
  summary_file <- file.path(script_dir, "cell_types_summary.csv")
  write_csv(file_summary, summary_file)
  cat("Saved summary to:", basename(summary_file), "\n\n")
  
  # Also save individual cell type gene lists for reference
  for (cell_type in names(all_gene_lists)) {
    if (length(all_gene_lists[[cell_type]]) > 0) {
      individual_file <- file.path(script_dir, paste0("cell_type_", cell_type, "_significant_genes.csv"))
      write_csv(data.frame(Gene_name = all_gene_lists[[cell_type]], stringsAsFactors = FALSE), individual_file)
    }
  }
  
  # INTERSECTION 1: Imm_CGE, Imm_LGE, Imm_MGE
  intersection1_genes <- compute_intersection_and_analyze(
    cell_type_names = c("Imm_CGE", "Imm_LGE", "Imm_MGE"),
    intersection_name = "Imm_CGE_Imm_LGE_Imm_MGE",
    intersection_number = 1,
    output_dir = intersections_dir
  )
  
  # INTERSECTION 2: Imm_CGE, Imm_LGE, Imm_MGE, Imm_CGE_LGE, Imm_MGE_MAF_CRABP1, Imm_MGE_MAF_NPY
  intersection2_genes <- compute_intersection_and_analyze(
    cell_type_names = c("Imm_CGE", "Imm_LGE", "Imm_MGE", "Imm_CGE_LGE", "Imm_MGE_MAF_CRABP1", "Imm_MGE_MAF_NPY"),
    intersection_name = "Imm_CGE_Imm_LGE_Imm_MGE_Imm_CGE_LGE_Imm_MGE_MAF_CRABP1_Imm_MGE_MAF_NPY",
    intersection_number = 2,
    output_dir = intersections_dir
  )
  
  # INTERSECTION 3: Previous cell types (from intersection 2) + Imm_Exc_2
  intersection3_genes <- compute_intersection_and_analyze(
    cell_type_names = c("Imm_CGE", "Imm_LGE", "Imm_MGE", "Imm_CGE_LGE", "Imm_MGE_MAF_CRABP1", "Imm_MGE_MAF_NPY", "Imm_Exc_2"),
    intersection_name = "Imm_CGE_Imm_LGE_Imm_MGE_Imm_CGE_LGE_Imm_MGE_MAF_CRABP1_Imm_MGE_MAF_NPY_Imm_Exc_2",
    intersection_number = 3,
    output_dir = intersections_dir
  )
  
  # INTERSECTION 4: Myel_Olig, Prem_Olig, and OPC
  intersection4_genes <- compute_intersection_and_analyze(
    cell_type_names = c("Myel_Olig", "Prem_Olig", "OPC"),
    intersection_name = "Myel_Olig_Prem_Olig_OPC",
    intersection_number = 4,
    output_dir = intersections_dir
  )
  
  # INTERSECTION 5: Ast_Rgl and Mcg
  intersection5_genes <- compute_intersection_and_analyze(
    cell_type_names = c("Ast_Rgl", "Mcg"),
    intersection_name = "Ast_Rgl_Mcg",
    intersection_number = 5,
    output_dir = intersections_dir
  )
  
  # INTERSECTION 6: Myel_Olig, Prem_Olig, OPC, Ast_Rgl, and Mcg
  intersection6_genes <- compute_intersection_and_analyze(
    cell_type_names = c("Myel_Olig", "Prem_Olig", "OPC", "Ast_Rgl", "Mcg"),
    intersection_name = "Myel_Olig_Prem_Olig_OPC_Ast_Rgl_Mcg",
    intersection_number = 6,
    output_dir = intersections_dir
  )
  
  # INTERSECTION 7: Intersection of Intersection 3 and Intersection 6
  # This computes the intersection of genes from intersection 3 and intersection 6
  if (!is.null(intersection3_genes) && !is.null(intersection6_genes)) {
    intersection7_genes <- intersect(intersection3_genes, intersection6_genes)
    
    if (length(intersection7_genes) > 0) {
      cat("\n", paste(rep("=", 70), collapse=""), "\n")
      cat("INTERSECTION 7: Intersection 3 + Intersection 6\n")
      cat("Common genes from Intersection 3 and Intersection 6:", length(intersection7_genes), "\n")
      cat(paste(rep("=", 70), collapse=""), "\n\n")
      
      # Save intersection to CSV
      safe_name <- "Intersection3_Intersection6"
      output_file <- file.path(intersections_dir, paste0("intersection_7_", safe_name, "_common_genes.csv"))
      write_csv(data.frame(Gene_name = intersection7_genes, stringsAsFactors = FALSE), output_file)
      cat("Saved common genes to:", basename(output_file), "\n")
      
      # Run enrichment analyses
      cat("\nRunning enrichment analyses...\n")
      output_name <- paste0("intersection_7_", safe_name)
      
      # Run enrichR
      if (enrichr_available) {
        enrichr_result <- run_enrichr_analysis(intersection7_genes, output_name, results_base_dir = intersections_dir)
      }
      
      # Run clusterProfiler
      if (clusterprofiler_available && orgdb_available) {
        clusterprofiler_result <- run_clusterprofiler_analysis(intersection7_genes, output_name, results_base_dir = intersections_dir)
      }
      
      cat("\n")
    } else {
      cat("\n", paste(rep("=", 70), collapse=""), "\n")
      cat("INTERSECTION 7: Intersection 3 + Intersection 6\n")
      cat("No common genes found between Intersection 3 and Intersection 6.\n")
      cat(paste(rep("=", 70), collapse=""), "\n\n")
      intersection7_genes <- NULL
    }
  } else {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("INTERSECTION 7: Intersection 3 + Intersection 6\n")
    cat("Cannot compute: One or both of Intersection 3 and Intersection 6 returned NULL.\n")
    cat(paste(rep("=", 70), collapse=""), "\n\n")
    intersection7_genes <- NULL
  }
  
} else {
  cat("No valid gene lists found.\n")
}

# ===== FINAL SUMMARY =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")
cat("Results saved:\n")
cat("  - cell_type_*_hits.txt (hits with Sig=1 and Shift_to_goal_end > 0 for each cell type)\n")
cat("  - cell_types_summary.csv (file summary)\n")
cat("  - Individual cell type enrichment results (enrichR and clusterProfiler)\n")
cat("  - Intersection 1: Imm_CGE, Imm_LGE, Imm_MGE (CSV, Venn diagram, enrichment)\n")
cat("  - Intersection 2: Imm_CGE, Imm_LGE, Imm_MGE, Imm_CGE_LGE, Imm_MGE_MAF_CRABP1, Imm_MGE_MAF_NPY (CSV, Venn diagram, enrichment)\n")
cat("  - Intersection 3: Previous (Intersection 2) + Imm_Exc_2 (CSV, Venn diagram, enrichment)\n")
cat("  - Intersection 4: Myel_Olig, Prem_Olig, and OPC (CSV, Venn diagram, enrichment)\n")
cat("  - Intersection 5: Ast_Rgl and Mcg (CSV, Venn diagram, enrichment)\n")
cat("  - Intersection 6: Myel_Olig, Prem_Olig, OPC, Ast_Rgl, and Mcg (CSV, Venn diagram, enrichment)\n")
cat("  - Intersection 7: Intersection of Intersection 3 and Intersection 6 (CSV, enrichment)\n")
cat("\n")

