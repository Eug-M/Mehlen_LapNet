# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

get_kegg_plots <- function(x) {
  gseaplot(gseaKEGG, geneSetID = gseaKEGG_results$ID[x], title = gseaKEGG_results$Description[x])
  pathview(gene.data = stats, pathway.id = gseaKEGG_results$ID[x], species = "hsa", 
           limit = list(gene = 2, cpd = 1))
  knitr::include_graphics(paste0(gseaKEGG_results$ID[x], ".png"))
}