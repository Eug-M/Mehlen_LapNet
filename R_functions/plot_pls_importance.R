# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

plot_pls_importance <- function(pls_results) {
  # Component 1
  imp1 <- pls_results$importance$Comp1
  # imp1 <- imp1[order(imp1$importance, decreasing = TRUE), ]
  imp1$component <- "Component 1"
  # Component 2
  imp2 <- pls_results$importance$Comp2
  # imp2 <- imp2[order(imp2$importance, decreasing = TRUE), ]
  imp2$component <- "Component 2"
  
  imp_combined <- rbind(imp1, imp2)
  
  ggplot(imp_combined, aes(x = reorder(variable, importance), y = importance, fill = component)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    scale_fill_manual(values = c("Component 1" = "steelblue", 
                                 "Component 2" = "coral")) +
    labs(title = "Variable Importance - PLS Components 1 and 2",
         x = "Variable",
         y = "Importance (%)",
         fill = "Component") +
    theme_minimal()+
    theme(legend.position = "top")
  
}