# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# Function to check if a design matrix is full rank
check_design_rank <- function(count_matrix, coldata, formula) {
  
  tryCatch({
    # Create model matrix
    mm <- model.matrix(formula, data = coldata)
    # Check rank
    is_full_rank <- qr(mm)$rank == ncol(mm)
    
    return(list(
      is_valid = is_full_rank,
      rank = qr(mm)$rank,
      ncol = ncol(mm)
    ))
  }, error = function(e) {
    return(list(
      is_valid = FALSE,
      error = e$message
    ))
  })
}