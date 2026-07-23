# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# Function to assess variable importance using PLS
assess_variable_pls <- function(coldata, rld, ntop = 500, ncomp = 2) {
  
  # Get expression data
  rv <- MatrixGenerics::rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  X <- t(assay(rld)[select, ])  # Samples x Genes
  # Prepare Y matrix - convert categorical to dummy variables
  Y_list <- list()
  variable_names <- colnames(coldata)
  
  for (var in variable_names) {
    var_data <- coldata[[var]]
    # Skip if all NA
    if (all(is.na(var_data))) next
    # Skip if ID or equivalent # on ne peut pas prendre les classifications, qui sont déterminées directement sur les gènes présents dans le RNAseq (tourne en rond)
    if (var == "SampleID" | var == "Date_traitement") next
    if (is.numeric(var_data) && !is.factor(var_data)) {
      next  # we want only the centered and reduced variables (factors)
    } else if (is.factor(var_data) && endsWith(var, '_cr')) {
      # Numeric variable - use as is (but centered and scaled)
      Y_list[[var]] <- as.numeric(as.vector(var_data))
    } else {
      # Categorical - convert to dummy variables
      var_factor <- base::as.factor(var_data)
      if (length(levels(var_factor)) >= 2) {
        # Create dummy variables (one-hot encoding)
        dummy <- model.matrix(~ var_factor - 1)
        colnames(dummy) <- paste0(var, "_", levels(var_factor))
        Y_list[[var]] <- dummy
      }
    }
    if (sum(is.na(var_data)) > 0) {
      if (is.matrix(Y_list[[var]])) {
        # Categorical variable stored as dummy matrix
        new_rows <- matrix(NA, 
                           nrow = length(which(is.na(var_data))), 
                           ncol = ncol(Y_list[[var]]),
                           dimnames = list(as.character(which(is.na(var_data))), colnames(Y_list[[var]]))
        )
        Y_list[[var]] <- rbind(Y_list[[var]], new_rows)
        Y_list[[var]] <- Y_list[[var]][order(as.numeric(rownames(Y_list[[var]]))), ]
      } else {
        # Numeric vector (_cr variables)
        Y_list[[var]][which(is.na(var_data))] <- NA
      }
    }
  }
  
  # Combine into Y matrix
  Y <- do.call(cbind, Y_list)
  
  # Remove rows with NA
  complete_cases <- complete.cases(X, Y)
  X_clean <- X[complete_cases, ]
  Y_clean <- Y[complete_cases, ]
  
  # Run PLS
  pls_result <- pls(X_clean, Y_clean, ncomp = ncomp, mode = "regression")
  
  # Calculate variable importance for each component
  importance_list <- list()
  
  for (comp in 1:ncomp) {
    # Get loadings for Y variables
    y_loadings <- pls_result$loadings$Y[, comp]
    
    # Calculate proportion of variance explained by each original variable
    var_importance <- data.frame(
      variable = character(),
      loading = numeric(),
      importance = numeric(),
      stringsAsFactors = FALSE
    )
    
    current_col <- 1
    for (var in variable_names) {
      if (var %in% names(Y_list)) {
        n_cols <- ncol(as.matrix(Y_list[[var]]))
        var_cols <- current_col:(current_col + n_cols - 1)
        
        # Sum of squared loadings for this variable
        var_loading <- sum(y_loadings[var_cols]^2)
        
        var_importance <- rbind(var_importance, 
                                data.frame(variable = var,
                                           loading = sqrt(var_loading),
                                           importance = var_loading * 100,
                                           stringsAsFactors = FALSE))
        
        current_col <- current_col + n_cols
      }
    }
    
    importance_list[[paste0("Comp", comp)]] <- var_importance
  }
  
  return(list(
    pls_result = pls_result,
    importance = importance_list,
    var_explained = pls_result$prop_expl_var$X,
    Y_matrix = Y_clean,
    X_matrix = X_clean
  ))
}