# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

create_classif_boxplots <- function(classif, classif_col, df1, df2) {
  p1 <- ggplot(df1, aes(x = .data[[classif_col]], y=Hallmark_score, fill=.data[[classif_col]])) + 
    geom_boxplot(notch = TRUE) +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(position=position_jitterdodge(dodge.width=0.9), size=0.4) +
    theme_bw()+ 
    theme(axis.text.x = element_text(angle=-40, hjust=.1)) +
    ggtitle(paste("Hallmark", classif, "LapNet")) 
  
  p2 <- ggplot(df2, aes(x = .data[[classif_col]], y=Hallmark_score, fill=.data[[classif_col]])) + 
    geom_boxplot(notch = TRUE) +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(position=position_jitterdodge(dodge.width=0.9), size=0.4) +
    theme_bw()+ 
    theme(axis.text.x = element_text(angle=-40, hjust=.1)) +
    ggtitle(paste("Hallmark", classif, "Nicolle"))
  
  p3 <- ggplot(df1, aes(x = .data[[classif_col]], y=Mak_score, fill=.data[[classif_col]])) + 
    geom_boxplot(notch = TRUE) +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(position=position_jitterdodge(dodge.width=0.9), size=0.4) +
    theme_bw()+ 
    theme(axis.text.x = element_text(angle=-40, hjust=.1)) +
    ggtitle(paste("Mak", classif, "LapNet"))
  
  p4 <- ggplot(df2, aes(x = .data[[classif_col]], y=Mak_score, fill=.data[[classif_col]])) + 
    geom_boxplot(notch = TRUE) +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(position=position_jitterdodge(dodge.width=0.9), size=0.4) +
    theme_bw()+ 
    theme(axis.text.x = element_text(angle=-40, hjust=.1)) +
    ggtitle(paste("Mak", classif, "Nicolle"))
  
  grid.arrange(p1, p2, p3, p4, nrow = 2)
}