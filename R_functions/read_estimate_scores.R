# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

read_estimate_scores <- function(gct_file) {
  
  # Read GCT file (skip first 2 lines which are headers)
  gct_data <- read.delim(gct_file, skip = 2, header = TRUE, row.names = 1)
  
  # The scores are typically in rows: StromalScore, ImmuneScore, ESTIMATEScore
  # Transpose to get samples as rows
  scores_df <- as.data.frame(t(gct_data[,which(colnames(gct_data) != "Description")]))
  
  # Clean up column names if needed
  colnames(scores_df) <- gsub("\\.", "", colnames(scores_df))
  
  return(scores_df)
}