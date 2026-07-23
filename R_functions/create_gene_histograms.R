# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

create_gene_histograms <- function(gene_name, df1, df2) {
  gene_log <- paste0(gene_name, "_log")
  
  p1 <- ggplot(df1, aes(x = .data[[gene_log]])) + 
    geom_histogram() +
    expand_limits(x = 0, y = 0) +
    ggtitle(paste0("log2(", gene_name, "+1) LapNet"))
  
  p2 <- ggplot(df2, aes(x = .data[[gene_log]])) + 
    geom_histogram() +
    expand_limits(x = 0, y = 0) +
    ggtitle(paste0("log2(", gene_name, "+1) Nicolle"))
  
  # grid.arrange(p1, p2, nrow = 1)
  return(list(plot1 = p1, plot2 = p2))
}