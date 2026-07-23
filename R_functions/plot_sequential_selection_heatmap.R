# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# Heatmap showing DEGs for each variable at each selection round
plot_sequential_selection_heatmap <- function(lrt_results, criterion = "padj_and_lfc") {
  
  results_df <- lrt_results$results
  
  # Choose which DEG count to use
  if (criterion == "padj_only") {
    results_df$n_deg <- results_df$n_deg_padj_only
    legend_title <- "DEGs (padj only)"
  } else {
    results_df$n_deg <- results_df$n_deg_padj_and_lfc
    legend_title <- "DEGs (padj + LFC)"
  }
  
  # Order variables by selection order
  variable_order <- lrt_results$selected_order
  if (length(variable_order) < length(unique(results_df$tested_variable))) {
    for (var in unique(results_df$tested_variable)) {
      if (var %!in% variable_order) {
        variable_order <- c(variable_order, var)
      }
    }
  }
  results_df$tested_variable <- factor(results_df$tested_variable, 
                                       levels = variable_order)
  
  # Create the plot
  p <- ggplot(results_df, aes(x = factor(round), y = tested_variable, fill = n_deg)) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = n_deg), color = "white", size = 4, fontface = "bold") +
    scale_fill_viridis(option = "plasma", name = legend_title,
                       na.value = "grey90") +
    labs(title = "Sequential Variable Selection: DEGs per Round",
         subtitle = "Variables ordered by selection rank (down = selected first)",
         x = "Selection Round",
         y = "Variable") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "right",
          panel.grid = element_blank())
  
  return(p)
}