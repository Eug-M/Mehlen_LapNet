# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

create_myPCAplots <- function(df, column, PCx, PCy, nber_top_genes) {
  pcaData <- plotPCA(rld, intgroup=column, pcsToUse=c(PCx,PCy), returnData=TRUE, ntop=nber_top_genes) #, ntop=1000)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  percent_x <- paste0(percentVar[1], "%")
  percent_y <- paste0(percentVar[2], "%")
  x <- paste0('PC', PCx)
  y <- paste0('PC', PCy)
  if (column == 'SampleID') {
    pca_plot <- ggplot(pcaData, aes(.data[[x]], .data[[y]], color=group, label=rownames(df))) +
      geom_point(size=3) +
      ggtitle(column) +
      # theme(plot.title = element_text(hjust=0.95, vjust=0.15, margin = margin(t = -20, b = 10, r = 5), size=11),
      theme(plot.title = element_text(margin = margin(t = 10, b = -20), size=11), legend.position="none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin = unit(c(1, 1, 1, 1), "mm")) +
      geom_text_repel(hjust="top", vjust="top") +
      geom_density2d(alpha=.5) +
      coord_fixed() 
    return(list("plot" = pca_plot, "PCx" = percent_x, "PCy" = percent_y))
  }
  else if (column %in% numeric_columns) {
    sc_color <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(df[,which(colnames(df) == column)]), max(df[,which(colnames(df) == column)])))
    pca_plot <- ggplot(pcaData, aes(.data[[x]], .data[[y]], color=group)) +
      geom_point(size=3) +
      ggtitle(column) +
      # theme(plot.title = element_text(hjust=0.95, vjust=0.15, margin = margin(t = -20, b = 10, r = 5), size=11),
      theme(plot.title = element_text(margin = margin(t = 10, b = -20), size=11), legend.position="none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin = unit(c(1, 1, 1, 1), "mm")) +
      geom_density2d(alpha=.5) +
      coord_fixed() +
      sc_color
    return(list("plot" = pca_plot, "PCx" = percent_x, "PCy" = percent_y))
  }
  else {
    pca_plot <- ggplot(pcaData, aes(.data[[x]], .data[[y]], color=group)) +
      geom_point(size=3) +
      ggtitle(column) +
      theme(plot.title = element_text(margin = margin(t = 10, b = -20), size=11), legend.position="none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin = unit(c(1, 1, 1, 1), "mm")) +
      geom_density2d(alpha=.5) +
      coord_fixed()
    return(list("plot" = pca_plot, "PCx" = percent_x, "PCy" = percent_y))
  }
}