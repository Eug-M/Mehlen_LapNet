# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

test_variable_influence_lrt_sequential_duplicates_withVarInt <- function(count_matrix,
                                                                         coldata,
                                                                         duplicates,
                                                                         variables = NULL,
                                                                         var_interest = NULL,
                                                                         padj_cutoff = 0.05,
                                                                         lfc_cutoff = 0.58) {
  
  # Get all variables from coldata if not specified
  if (is.null(variables)) {
    variables <- colnames(coldata)
  } else {
    missing_vars <- dplyr::setdiff(variables, colnames(coldata))
    if (length(missing_vars) > 0) {
      stop("These variables not found in coldata: ", paste(missing_vars, collapse = ", "))
    }
  }
  duplicates_list <- do.call(c, duplicates)
  
  # Storage for all results
  all_rounds_results <- list()
  selected_variables <- var_interest
  remaining_variables <- variables
  
  round_num <- 1
  check_zero_deg <- 1
  
  # Continue until all variables are tested
  while (length(remaining_variables) > 0 & check_zero_deg == 1) {
    
    # cat("\n========================================\n")
    # cat("ROUND", round_num, "\n")
    # cat("========================================\n")
    
    if (round_num == 1) {
      # cat("Testing variables alone (reduced = ~1)\n\n")
      round_results <- list()
      reduced_formula <- "~1"
      full_formula <- as.formula(paste("~", var_interest))
      rank_check <- check_design_rank(count_matrix, coldata, full_formula)
      dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = coldata,
        design = full_formula
      )
      # Run LRT
      dds_lrt <- DESeq(dds, test = "LRT", 
                       reduced = as.formula(reduced_formula), 
                       quiet = TRUE)
      res <- results(dds_lrt)
      
      # Count DEGs
      n_deg_padj <- sum(res$padj < padj_cutoff, na.rm = TRUE)
      n_deg_both <- sum(res$padj < padj_cutoff & abs(res$log2FoldChange) >= lfc_cutoff, 
                        na.rm = TRUE)
      
      round_results[[var_interest]] <- data.frame(
        round = round_num,
        tested_variable = var_interest,
        current_design = "none",
        full_model = paste("~", var_interest),
        reduced_model = reduced_formula,
        n_deg_padj_only = n_deg_padj,
        n_deg_padj_and_lfc = n_deg_both,
        stringsAsFactors = FALSE
      )
      round_df <- do.call(rbind, round_results)
      rownames(round_df) <- NULL
      
      # Select best variable (most DEGs with padj + LFC)
      best_variable <- round_df$tested_variable[var_interest]
      best_n_deg <- round_df$n_deg_padj_and_lfc[var_interest]
      all_rounds_results[[paste0("round_", round_num)]] <- round_df
      
      round_num <- 2
      reduced_formula <- paste("~", var_interest)
      
    } else {
      # cat("Current design:", paste("~", paste(selected_variables, collapse = " + ")), "\n")
      # cat("Testing addition of remaining variables\n\n")
      reduced_formula <- paste("~", paste(selected_variables, collapse = " + "))
    }
    
    # Storage for this round
    round_results <- list()
    skipped_tests <- list()
    
    # Test each remaining variable
    for (var in remaining_variables) {
      
      full_formula_str <- paste("~", paste(c(var, selected_variables), collapse = " + "))
      # cat("Testing:", var, "added to current design\n")
      
      # Check design validity
      full_formula <- as.formula(full_formula_str)
      rank_check <- check_design_rank(count_matrix, coldata, full_formula)
      
      if (!rank_check$is_valid) {
        cat("  SKIPPED: Design matrix not full rank\n")
        skipped_tests[[var]] <- "Not full rank"
        next
      }
      
      tryCatch({
        # Create dds
        dds <- DESeqDataSetFromMatrix(
          countData = count_matrix,
          colData = coldata,
          design = full_formula
        )
        
        # Run LRT
        dds_lrt <- DESeq(dds, test = "LRT", 
                         reduced = as.formula(reduced_formula), 
                         quiet = TRUE)
        res <- results(dds_lrt)
        
        # Count DEGs
        n_deg_padj <- sum(res$padj < padj_cutoff, na.rm = TRUE)
        n_deg_both <- sum(res$padj < padj_cutoff & abs(res$log2FoldChange) >= lfc_cutoff, 
                          na.rm = TRUE)
        
        round_results[[var]] <- data.frame(
          round = round_num,
          tested_variable = var,
          current_design = ifelse(round_num == 1, "none", 
                                  paste(selected_variables, collapse = " + ")),
          full_model = full_formula_str,
          reduced_model = reduced_formula,
          n_deg_padj_only = n_deg_padj,
          n_deg_padj_and_lfc = n_deg_both,
          stringsAsFactors = FALSE
        )
        
        # cat("  DEGs (padj):", n_deg_padj, "| DEGs (padj + LFC):", n_deg_both, "\n")
        
      }, error = function(e) {
        cat("  ERROR:", e$message, "\n")
        skipped_tests[[var]] <- e$message
      })
    }
    
    # Combine results for this round
    if (length(round_results) == 0) {
      cat("\nNo successful tests in this round. Stopping.\n")
      break
    }
    
    round_df <- do.call(rbind, round_results)
    rownames(round_df) <- NULL
    
    # Select best variable (most DEGs with padj + LFC)
    best_idx <- order(round_df$n_deg_padj_and_lfc, round_df$n_deg_padj_only, 
                      decreasing = TRUE)[1] # which.max(round_df$n_deg_padj_and_lfc)
    best_variable <- round_df$tested_variable[best_idx]
    best_n_deg <- round_df$n_deg_padj_and_lfc[best_idx]
    
    if (sum(round_df$n_deg_padj_and_lfc) == 0 & sum(round_df$n_deg_padj_only) == 0) {
      cat("\nNo variable gave at least 1 DEG with the LRT method. Stopping.\n")
      check_zero_deg <- 0
    }
    
    # cat("\n--- ROUND", round_num, "SUMMARY ---\n")
    # cat("Best variable:", best_variable, "with", best_n_deg, "DEGs\n")
    
    # Add to results
    all_rounds_results[[paste0("round_", round_num)]] <- round_df
    
    # Update for next round
    selected_variables <- c(best_variable, selected_variables)
    if (best_variable %in% duplicates_list) {
      which_sublist <- which(sapply(duplicates, function(x) best_variable %in% x))
      remaining_variables <- dplyr::setdiff(remaining_variables, duplicates[[which_sublist]])
    }
    else {
      remaining_variables <- dplyr::setdiff(remaining_variables, best_variable)
    }
    
    round_num <- round_num + 1
  }
  
  # Combine all rounds
  final_results <- do.call(rbind, all_rounds_results)
  rownames(final_results) <- NULL
  
  # Summary
  # cat("\n========================================\n")
  # cat("FINAL SUMMARY\n")
  # cat("========================================\n")
  # cat("Variables selected in order:\n")
  for (i in seq_along(selected_variables)) {
    var_info <- final_results[final_results$round == i & 
                                final_results$tested_variable == rev(selected_variables)[i], ]
    # cat(i, ".", rev(selected_variables)[i], "- DEGs:",
    #     var_info$n_deg_padj_and_lfc, "\n")
  }
  
  cat("\nFinal design:", paste("~", paste(selected_variables, collapse = " + ")), "\n")
  
  return(list(
    results = final_results,
    selected_order = rev(selected_variables),
    selection_summary = data.frame(
      rank = 1:length(selected_variables),
      variable = rev(selected_variables),
      n_deg = sapply(rev(selected_variables), function(v) {
        final_results[final_results$tested_variable == v, "n_deg_padj_and_lfc"][1]
      })
    )
  ))
}