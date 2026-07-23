# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later
apply_gene_symbol_mapping <- function(count_matrix, mapping_file, df_genelengths = NULL) {
  
  # cat("=== APPLYING GENE SYMBOL MAPPING TO COUNT MATRIX ===\n\n")
  
  # Read mapping if it's a file path
  if (is.character(mapping_file)) {
    mapping <- read.csv(mapping_file, stringsAsFactors = FALSE)
  } else {
    mapping <- mapping_file
  }
  
  # cat("Mapping contains", nrow(mapping), "gene symbols\n")
  # cat("Count matrix contains", nrow(count_matrix), "genes\n\n")
  
  # gene_name <- mapping$external_gene_name[i] 
  # Pre-compute which IDs exist in count_matrix 
  available_ids <- rownames(count_matrix)
  available_set <- as.list(setNames(rep(TRUE, length(available_ids)), available_ids))
  
  # Separate regular and merged
  # is_merged <- grepl("_merged", mapping$final_gene_id)
  is_merged <- mapping$resolution_strategy == 'summed'
  regular_genes <- mapping[!is_merged, ]
  merged_genes <- mapping[is_merged, ]
  
  # Process regular genes 
  regular_ids_to_keep <- regular_genes[match(available_ids, regular_genes$final_gene_id, nomatch = 0), 'final_gene_id']
  # regular_ids_to_keep <- regular_genes$final_gene_id[regular_genes$final_gene_id %in% available_ids]
  
  if (is.null(df_genelengths)) {
    new_matrix <- as.data.frame(count_matrix[regular_ids_to_keep, , drop = FALSE])
  } else {
    new_matrix <- as.data.frame(count_matrix[regular_ids_to_keep, , drop = FALSE])
    new_matrix$gene_lengths <- df_genelengths$gene_length[match(regular_ids_to_keep, rownames(df_genelengths))]
  }
  rownames(new_matrix) <- regular_genes[match(available_ids, regular_genes$final_gene_id, nomatch = 0), 'external_gene_name']
  # rownames(new_matrix) <- regular_genes$external_gene_name[regular_genes$final_gene_id %in% available_ids]
  
  # cat("Kept", length(regular_ids_to_keep), "regular genes\n")
  
  # Process merged genes (only if they exist)
  if (nrow(merged_genes) > 0) {
    # cat("Processing", nrow(merged_genes), "merged genes...\n")
    
    # Pre-allocate matrix
    n_samples <- ncol(count_matrix)
    if (is.null(df_genelengths)) {
      merged_matrix <- as.data.frame(matrix(0, nrow = nrow(merged_genes), ncol = n_samples,
                                            dimnames = list(merged_genes$external_gene_name,
                                                            colnames(count_matrix))))
    } else {
      merged_matrix <- as.data.frame(matrix(0, nrow = nrow(merged_genes), ncol = n_samples+1,
                                            dimnames = list(merged_genes$external_gene_name,
                                                            c(colnames(count_matrix), "gene_lengths"))))
    }
    
    # Loop only through merged genes 
    rows_to_remove <- list()
    for (i in seq_len(nrow(merged_genes))) {
      original_ids <- unlist(strsplit(merged_genes$ensembl_ids_to_use[[i]], ';'))
      
      # Fast check which IDs exist
      existing_ids <- original_ids[original_ids %in% available_ids]
      
      if (length(existing_ids) > 0) {
        if (is.null(df_genelengths)) {
          merged_matrix[i, ] <- colSums(count_matrix[existing_ids, , drop = FALSE])
        } else {
          merged_matrix[i, ] <- c(colSums(count_matrix[existing_ids, , drop = FALSE]), 
                                  sum(df_genelengths$gene_length[match(existing_ids, rownames(df_genelengths))]))
        }
      } else {
        rows_to_remove <- c(rows_to_remove, -i)
      }
    }
    rows_to_remove <- unlist(rows_to_remove)
    if (length(rows_to_remove) > 0) {
      merged_matrix <- merged_matrix[rows_to_remove,]
    }
    
    # Combine
    new_matrix <- rbind(new_matrix, merged_matrix)
  }
  
  # cat("Output:", nrow(new_matrix), "genes x", ncol(new_matrix), "samples\n")
  
  return(new_matrix)
}