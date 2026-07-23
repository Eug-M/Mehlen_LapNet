# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

figures_ORA <- function(results, name_geneset, type_DEG) {
  ## from: https://biostatsquid.com/pathway-enrichment-analysis-plots/
  # Barplot
  # Most common method to visualize enriched terms. Shows enrichment scores (e.g. p values) and gene count or ratio as bar height and color 
  p1 <- barplot(results, showCategory = 12, 
                title = paste("MSigDB", name_geneset, "Enrichment -", type_DEG, "regulated genes"))
  p2 <- mutate(results, qscore = -log(p.adjust, base = 10)) %>% 
    barplot(x = "qscore")
  # Dotplot
  # Dot plot is similar to bar plot with the capability to encode another score as dot size.
  p3 <- enrichplot::dotplot(results, showCategory = 15)
  # cnetplot
  # Shows the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network.
  p4 <- cnetplot(results)
  # Heatplot
  # The heatplot is similar to cnetplot, while displaying the relationships as a heatmap. useful when there is a large number of significant terms
  p5 <- heatplot(results, showCategory = 12)
  # Treeplot
  # The treeplot() function performs hierarchical clustering of enriched terms. 
  enrichres2 <- enrichplot::pairwise_termsim(results) # calculate pairwise similarities of the enriched terms using Jaccard’s similarity index
  p6 <- enrichplot::treeplot(enrichres2)
  # Enrichment map 
  # organizes enriched pathways into a network with edges connecting overlapping gene sets. In this way, mutually overlapping gene sets are tend to cluster together, making it easy to identify functional module.
  p7 <- emapplot(enrichres2)
  # Upsetplot
  # alternative to cnetplot. emphasizes the gene overlapping among different gene sets.
  p8 <- enrichplot::upsetplot(results)
  
  return(list(plot1 = p1, plot2 = p2, plot3 = p3, plot4 = p4, plot5 = p5, plot6 = p6, plot7 = p7, plot8 = p8))
}