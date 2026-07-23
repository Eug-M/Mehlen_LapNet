# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

calculate_tumor_spline <- function(coldata, 
                                   size_cols = c("Size_norm_C0", "Size_norm_C4", 
                                                 "Size_norm_C8", "Size_norm_C12"),
                                   time_cols = c("Days_lesion_C0", "Days_lesion_C4", 
                                                 "Days_lesion_C8", "Days_lesion_C12"),
                                   df = 2) {  # degrees of freedom for spline
  
  coldata$tumor_spline_slope <- NA
  coldata$tumor_spline_curvature <- NA
  
  for (i in 1:nrow(coldata)) {
    sizes <- as.numeric(coldata[i, size_cols])
    time_points <- as.numeric(coldata[i, time_cols])
    valid_idx <- !is.na(sizes)
    
    if (sum(valid_idx) >= 3) {  # Need at least 3 points for spline
      valid_sizes <- sizes[valid_idx]
      valid_times <- time_points[valid_idx]
      
      tryCatch({
        # Fit spline
        spline_fit <- lm(valid_sizes ~ ns(valid_times, df = min(df, sum(valid_idx) - 1)))
        
        # Extract overall trend (coefficient of first spline basis)
        coldata$tumor_spline_slope[i] <- coef(spline_fit)[2]
        
        # If enough df, get curvature (second derivative approximation)
        if (length(coef(spline_fit)) > 2) {
          coldata$tumor_spline_curvature[i] <- coef(spline_fit)[3]
        }
      }, error = function(e) {
        # Skip if spline fails
      })
    } else if (sum(valid_idx) == 2) {
      # Fall back to linear for 2 points
      valid_sizes <- sizes[valid_idx]
      valid_times <- time_points[valid_idx]
      lm_fit <- lm(valid_sizes ~ valid_times)
      coldata$tumor_spline_slope[i] <- coef(lm_fit)[2]
    }
  }
  
  return(coldata)
}