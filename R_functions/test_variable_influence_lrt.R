# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

test_variable_influence_lrt <- function(count_matrix,
                                        coldata,
                                        variables = NULL,
                                        padj_cutoff = 0.05,
                                        lfc_cutoff = 0.58,
                                        max_combination_size = 3) {
  
  # Get all variables from coldata if not specified
  if (is.null(variables)) {
    # Try to identify likely variables (exclude sample names, IDs, etc.)
    variables <- colnames(coldata)
    # cat("Available variables in coldata:", paste(variables, collapse = ", "), "\n")
    # cat("Testing all variables\n")
  } else {
    # Check that all specified variables exist
    missing_vars <- dplyr::setdiff(variables, colnames(coldata))
    if (length(missing_vars) > 0) {
      stop("These variables not found in coldata: ", paste(missing_vars, collapse = ", "))
    }
  }
  
  # cat("\nVariables to test:", paste(variables, collapse = ", "), "\n")
  # cat("padj cutoff:", padj_cutoff, "\n")
  # cat("LFC cutoff:", lfc_cutoff, "\n\n")
  
  # Storage for results
  results_list <- list()
  skipped_tests <- list()
  
  # 1. Test each variable alone (reduced = ~1)
  # cat("=== Testing individual variables ===\n")
  for (var in variables) {
    # cat("Testing:", var, "vs ~1\n")
    
    # Check design validity
    formula_full <- as.formula(paste("~", var))
    rank_check <- check_design_rank(count_matrix, coldata, formula_full)
    
    if (!rank_check$is_valid) {
      # cat("  SKIPPED: Design matrix not full rank\n")
      skipped_tests[[paste(var, "vs", "~1")]] <- "Not full rank"
      next
    }
    
    tryCatch({
      # Create dds for this test
      dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = coldata,
        design = formula_full
      )
      
      # Run LRT
      dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1, quiet = TRUE)
      res <- results(dds_lrt)
      
      # Count DEGs
      n_deg_padj <- sum(res$padj < padj_cutoff, na.rm = TRUE)
      n_deg_both <- sum(res$padj < padj_cutoff & abs(res$log2FoldChange) >= lfc_cutoff, 
                        na.rm = TRUE)
      
      results_list[[var]] <- data.frame(
        full_model = var,
        reduced_model = "~1",
        test_type = "single_variable",
        tested_variable = var,
        n_variables = 1,
        n_deg_padj_only = n_deg_padj,
        n_deg_padj_and_lfc = n_deg_both,
        stringsAsFactors = FALSE
      )
      
      # cat("  DEGs (padj):", n_deg_padj, "| DEGs (padj + LFC):", n_deg_both, "\n")
      
    }, error = function(e) {
      cat("  ERROR:", e$message, "\n")
      skipped_tests[[paste(var, "vs", "~1")]] <- e$message
    })
  }
  
  # 2. Test pairs of variables (full = A + B, reduced = A or B)
  if (length(variables) >= 2 && max_combination_size >= 2) {
    # cat("\n=== Testing pairs of variables ===\n")
    
    for (i in 1:(length(variables) - 1)) {
      for (j in (i + 1):length(variables)) {
        var1 <- variables[i]
        var2 <- variables[j]
        
        # Test A + B vs A (to see effect of B)
        # cat("Testing:", var1, "+", var2, "vs", var1, "\n")
        
        formula_full <- as.formula(paste("~", var1, "+", var2))
        formula_reduced <- as.formula(paste("~", var1))
        
        # Check design validity
        rank_check_full <- check_design_rank(count_matrix, coldata, formula_full)
        rank_check_reduced <- check_design_rank(count_matrix, coldata, formula_reduced)
        
        if (!rank_check_full$is_valid || !rank_check_reduced$is_valid) {
          # cat("  SKIPPED: Design matrix not full rank\n")
          skipped_tests[[paste(var1, var2, "vs", var1)]] <- "Not full rank"
          next
        }
        
        tryCatch({
          dds <- DESeqDataSetFromMatrix(
            countData = count_matrix,
            colData = coldata,
            design = formula_full
          )
          
          dds_lrt <- DESeq(dds, test = "LRT", reduced = formula_reduced, quiet = TRUE)
          res <- results(dds_lrt)
          
          n_deg_padj <- sum(res$padj < padj_cutoff, na.rm = TRUE)
          n_deg_both <- sum(res$padj < padj_cutoff & abs(res$log2FoldChange) >= lfc_cutoff, 
                            na.rm = TRUE)
          
          results_list[[paste(var1, var2, "vs", var1, sep = "_")]] <- data.frame(
            full_model = paste(var1, "+", var2),
            reduced_model = var1,
            test_type = "added_variable",
            tested_variable = var2,
            n_variables = 2,
            n_deg_padj_only = n_deg_padj,
            n_deg_padj_and_lfc = n_deg_both,
            stringsAsFactors = FALSE
          )
          
          # cat("  DEGs (padj):", n_deg_padj, "| DEGs (padj + LFC):", n_deg_both, "\n")
          
        }, error = function(e) {
          # cat("  ERROR:", e$message, "\n")
          skipped_tests[[paste(var1, var2, "vs", var1)]] <- e$message
        })
        
        # Test A + B vs B (to see effect of A)
        # cat("Testing:", var1, "+", var2, "vs", var2, "\n")
        
        formula_reduced2 <- as.formula(paste("~", var2))
        rank_check_reduced2 <- check_design_rank(count_matrix, coldata, formula_reduced2)
        
        if (!rank_check_full$is_valid || !rank_check_reduced2$is_valid) {
          # cat("  SKIPPED: Design matrix not full rank\n")
          skipped_tests[[paste(var1, var2, "vs", var2)]] <- "Not full rank"
          next
        }
        
        tryCatch({
          dds <- DESeqDataSetFromMatrix(
            countData = count_matrix,
            colData = coldata,
            design = formula_full
          )
          
          dds_lrt <- DESeq(dds, test = "LRT", reduced = formula_reduced2, quiet = TRUE)
          res <- results(dds_lrt)
          
          n_deg_padj <- sum(res$padj < padj_cutoff, na.rm = TRUE)
          n_deg_both <- sum(res$padj < padj_cutoff & abs(res$log2FoldChange) >= lfc_cutoff, 
                            na.rm = TRUE)
          
          results_list[[paste(var1, var2, "vs", var2, sep = "_")]] <- data.frame(
            full_model = paste(var1, "+", var2),
            reduced_model = var2,
            test_type = "added_variable",
            tested_variable = var1,
            n_variables = 2,
            n_deg_padj_only = n_deg_padj,
            n_deg_padj_and_lfc = n_deg_both,
            stringsAsFactors = FALSE
          )
          
          # cat("  DEGs (padj):", n_deg_padj, "| DEGs (padj + LFC):", n_deg_both, "\n")
          
        }, error = function(e) {
          cat("  ERROR:", e$message, "\n")
          skipped_tests[[paste(var1, var2, "vs", var2)]] <- e$message
        })
      }
    }
  }
  
  # 3. Test triplets (full = A + B + C, reduced = A + B)
  if (length(variables) >= 3 && max_combination_size >= 3) {
    # cat("\n=== Testing triplets of variables ===\n")
    
    for (i in 1:(length(variables) - 2)) {
      for (j in (i + 1):(length(variables) - 1)) {
        for (k in (j + 1):length(variables)) {
          var1 <- variables[i]
          var2 <- variables[j]
          var3 <- variables[k]
          
          # cat("Testing:", var1, "+", var2, "+", var3, "vs", var1, "+", var2, "\n")
          
          formula_full <- as.formula(paste("~", var1, "+", var2, "+", var3))
          formula_reduced <- as.formula(paste("~", var1, "+", var2))
          
          # Check design validity
          rank_check_full <- check_design_rank(count_matrix, coldata, formula_full)
          rank_check_reduced <- check_design_rank(count_matrix, coldata, formula_reduced)
          
          if (!rank_check_full$is_valid || !rank_check_reduced$is_valid) {
            # cat("  SKIPPED: Design matrix not full rank\n")
            skipped_tests[[paste(var1, var2, var3, "vs", var1, var2)]] <- "Not full rank"
            next
          }
          
          tryCatch({
            dds <- DESeqDataSetFromMatrix(
              countData = count_matrix,
              colData = coldata,
              design = formula_full
            )
            
            dds_lrt <- DESeq(dds, test = "LRT", reduced = formula_reduced, quiet = TRUE)
            res <- results(dds_lrt)
            
            n_deg_padj <- sum(res$padj < padj_cutoff, na.rm = TRUE)
            n_deg_both <- sum(res$padj < padj_cutoff & abs(res$log2FoldChange) >= lfc_cutoff, 
                              na.rm = TRUE)
            
            results_list[[paste(var1, var2, var3, "vs", var1, var2, sep = "_")]] <- data.frame(
              full_model = paste(var1, "+", var2, "+", var3),
              reduced_model = paste(var1, "+", var2),
              test_type = "added_to_pair",
              tested_variable = var3,
              n_variables = 3,
              n_deg_padj_only = n_deg_padj,
              n_deg_padj_and_lfc = n_deg_both,
              stringsAsFactors = FALSE
            )
            
            # cat("  DEGs (padj):", n_deg_padj, "| DEGs (padj + LFC):", n_deg_both, "\n")
            
          }, error = function(e) {
            cat("  ERROR:", e$message, "\n")
            skipped_tests[[paste(var1, var2, var3, "vs", var1, var2)]] <- e$message
          })
        }
      }
    }
  }
  
  # Combine all results
  if (length(results_list) == 0) {
    cat("\nWARNING: No tests were successful!\n")
    results_df <- data.frame()
  } else {
    results_df <- do.call(rbind, results_list)
    rownames(results_df) <- NULL
  }
  
  # Summary
  # cat("\n=== SUMMARY ===\n")
  # cat("Successful tests:", length(results_list), "\n")
  # cat("Skipped/failed tests:", length(skipped_tests), "\n")
  
  if (length(skipped_tests) > 0) {
    cat("\nSkipped tests:\n")
    for (test_name in names(skipped_tests)) {
      cat("  ", test_name, ":", skipped_tests[[test_name]], "\n")
    }
  }
  
  return(list(
    results = results_df,
    skipped = skipped_tests
  ))
}