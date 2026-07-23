# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# Function to assess variable contribution to PCs
assess_variable_contribution <- function(rld, pca, ntop = 500, variables = NULL) {
  
  # If variables not specified, use all from colData
  if (is.null(variables)) { variables <- colnames(colData(rld)) }
  # # Calculate PCA # already done above, so added it to stdin
  # rv <- MatrixGenerics::rowVars(assay(rld))
  # select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  # pca <- prcomp(t(assay(rld)[select, ]))
  # Get PC coordinates
  pca_data <- data.frame(pca$x[, 1:2])
  pca_data <- cbind(pca_data, colData(rld))
  # Calculate variance explained
  percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100
  # For each variable, calculate correlation/association with PCs
  results_list <- list()
  
  for (var in variables) {
    # print(var)
    # Skip if variable doesn't exist or has only one level
    if (!var %in% colnames(pca_data)) next
    var_data <- pca_data[[var]]
    # Check if variable is numeric or categorical
    if (is.numeric(var_data) && !is.factor(var_data)) {
      next  # we want only the centered and reduced variables (factors)
      # # For numeric: calculate correlation
      # cor_pc1 <- cor(pca_data$PC1, var_data, use = "complete.obs")
      # cor_pc2 <- cor(pca_data$PC2, var_data, use = "complete.obs")
      # # R-squared as measure of variance explained
      # r2_pc1 <- cor_pc1^2 * 100
      # r2_pc2 <- cor_pc2^2 * 100
    } else if (is.factor(var_data) && endsWith(var, '_cr')) {
      var_data <- as.numeric(as.vector(pca_data[[var]]))
      cor_pc1 <- cor(pca_data$PC1, var_data, use = "complete.obs")
      cor_pc2 <- cor(pca_data$PC2, var_data, use = "complete.obs")
      # R-squared as measure of variance explained
      r2_pc1 <- cor_pc1^2 * 100
      r2_pc2 <- cor_pc2^2 * 100
    }
    else {
      # For categorical: ANOVA / calculate eta-squared
      var_data <- base::as.factor(var_data)
      # Skip if only one level
      if (length(levels(var_data)) < 2) next
      # ANOVA for PC1
      aov_pc1 <- aov(pca_data$PC1 ~ var_data)
      ss_total_pc1 <- sum((pca_data$PC1 - mean(pca_data$PC1))^2)
      ss_var_pc1 <- sum(aov_pc1$residuals^2)
      r2_pc1 <- (1 - ss_var_pc1/ss_total_pc1) * 100
      # ANOVA for PC2
      aov_pc2 <- aov(pca_data$PC2 ~ var_data)
      ss_total_pc2 <- sum((pca_data$PC2 - mean(pca_data$PC2))^2)
      ss_var_pc2 <- sum(aov_pc2$residuals^2)
      r2_pc2 <- (1 - ss_var_pc2/ss_total_pc2) * 100
    }
    
    results_list[[var]] <- data.frame(
      variable = var,
      PC1_variance_explained = r2_pc1,
      PC2_variance_explained = r2_pc2
    )
  }
  
  # Combine results
  contribution_df <- do.call(rbind, results_list)
  rownames(contribution_df) <- NULL
  # Sort by PC1 contribution
  contribution_df <- contribution_df[order(contribution_df$PC1_variance_explained, 
                                           decreasing = TRUE), ]
  
  return(list(
    contribution = contribution_df,
    percentVar = percentVar,
    pca_data = pca_data
  ))
}