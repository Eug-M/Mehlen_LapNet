# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# Calculate tumor purity using the formula from the ESTIMATE paper: Purity = cos(0.6049872018 + 0.0001467884 * ESTIMATE score)
calculate_tumor_purity <- function(estimate_scores) {
  
  # Extract ESTIMATE score column
  if ("ESTIMATEScore" %in% colnames(estimate_scores)) {
    estimate_score <- as.numeric(estimate_scores$ESTIMATEScore)
  } else if ("ESTIMATE.score" %in% colnames(estimate_scores)) {
    estimate_score <- as.numeric(estimate_scores$ESTIMATE.score)
  } else {
    # Find column with "ESTIMATE" in name
    estimate_col <- grep("ESTIMATE", colnames(estimate_scores), ignore.case = TRUE, value = TRUE)[1]
    estimate_score <- as.numeric(estimate_scores[[estimate_col]])
  }
  
  # Calculate purity
  purity <- cos(0.6049872018 + 0.0001467884 * estimate_score)
  
  # Add to dataframe
  estimate_scores$TumorPurity <- purity
  
  return(estimate_scores)
}