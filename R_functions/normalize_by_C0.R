# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

normalize_by_C0 <- function(coldata,
                            size_cols = c("Size_lesion_C0", "Size_lesion_C4", 
                                          "Size_lesion_C8", "Size_lesion_C12"),
                            time_cols = c("Days_lesion_C0", "Days_lesion_C4", 
                                          "Days_lesion_C8", "Days_lesion_C12")) {
  
  coldata$Size_norm_C0 <- NA
  coldata$Size_norm_C4 <- NA
  coldata$Size_norm_C8 <- NA
  coldata$Size_norm_C12 <- NA
  
  for (i in 1:nrow(coldata)) {
    factor_norm <- 100 / coldata$Size_lesion_C0[i]
    coldata$Size_norm_C0[i] <- 100
    coldata$Size_norm_C4[i] <- coldata$Size_lesion_C4[i] * factor_norm
    coldata$Size_norm_C8[i] <- coldata$Size_lesion_C8[i] * factor_norm
    coldata$Size_norm_C12[i] <- coldata$Size_lesion_C12[i] * factor_norm
  }
  return(coldata)
}