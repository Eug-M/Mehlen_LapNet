# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# using enricher from clusterProfiler 
clustP_query <- function(genes, genesets, universe = NULL, adj_p_cutoff) { 
  result <- list() 
  if (is.vector(genes) && length(genes) > 0 && !all(genes == '')) 
  { 
    universe.gene.pathway <- universe 
    other.genes <- universe.gene.pathway[!universe.gene.pathway %in% genesets$gene_symbol] 
    other.genes <- na.omit(unique(other.genes)) 
    other.genes.df <- cbind(rep("other.genes", length(other.genes)), other.genes) 
    colnames(other.genes.df) <- colnames(genesets) 
    genes.sets.all <- rbind(genesets, other.genes.df)
    
    res <-
      enricher(
        genes,
        universe = universe,
        TERM2GENE = genes.sets.all,
        qvalueCutoff = adj_p_cutoff,
        maxGSSize=10000
      )
    result <- as.data.frame(res)
  } else {
    warning(paste0('genes must be a non-empty vector of gene names'))
  }
  return(result)
}