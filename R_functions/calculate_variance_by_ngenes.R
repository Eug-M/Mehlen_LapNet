# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# Function to calculate variance explained by top N genes
calculate_variance_by_ngenes <- function(rld, gene_numbers = seq(100, 2000, by = 100)) {
  
  # Get all gene variances
  rv <- MatrixGenerics::rowVars(assay(rld))
  results <- data.frame(
    n_genes = integer(),
    variance_explained = numeric()
  )
  
  for (n in gene_numbers) {
    # Select top n genes
    select <- order(rv, decreasing = TRUE)[seq_len(min(n, length(rv)))]
    # Run PCA on these genes
    pca <- prcomp(t(assay(rld)[select, ]))
    # Calculate total variance explained by first 2 PCs
    percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100
    total_var_pc1_pc2 <- sum(percentVar[1:2])
    # Or you can use variance explained by all PCs (should be close to 100%)
    # For selection, typically we care about PC1+PC2
    
    results <- rbind(results, data.frame(
      n_genes = n,
      variance_explained = total_var_pc1_pc2
    ))
  }
  return(results)
}