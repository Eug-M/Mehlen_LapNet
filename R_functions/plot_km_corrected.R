# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# Function to create KM plots OS & PFS with corrected p-value
plot_km_corrected <- function(data, origin, gene_column, type = 'OS', n_tests = 1) {
  
  col_level <- paste0(unlist(strsplit(gene_column, split='_', fixed=TRUE))[1], '_level_', origin)
  data[[col_level]] <- case_when(
    data[[gene_column]] > median(data[which(data$Origin == origin),gene_column]) ~ "High",
    data[[gene_column]] <= median(data[which(data$Origin == origin),gene_column]) ~ "Low")
  # Create formula
  if (type == 'OS') {
    formula_surv <- as.formula(paste("Surv(OS, Death_status) ~", col_level))
    y_lab = "Overall Survival"
  }
  else {
    formula_surv <- as.formula(paste("Surv(PFS, Prog_status) ~", col_level))
    y_lab = "Progression Free Survival"
  }
  if (origin == "LapNet") {
    break_time = 50
  }
  else {
    break_time = 100
  }
  tryCatch({
    # Fit KM
    km_fit <- surv_fit(formula_surv, data = data[which(data$Origin == origin),])
    # Calculate p-value
    logrank_test <- survdiff(formula_surv, data = data[which(data$Origin == origin),])
    pval_raw <- pchisq(logrank_test$chisq, 
                       df = length(logrank_test$n) - 1, 
                       lower.tail = FALSE)
    # Bonferroni correction
    pval_corrected <- min(pval_raw * n_tests, 1)
    # Create label
    pval_label <- paste0(
      "p = ", format.pval(pval_raw, digits = 3), 
      "\nBonf. adj. p = ", format.pval(pval_corrected, digits = 3)
    )
    # Legend order
    legend_order <- gsub(paste0(col_level, "="), "", names(km_fit$strata))
    # Plot
    p1 <- survminer::ggsurvplot(
      km_fit, 
      data = data[which(data$Origin == origin),], 
      conf.int = TRUE,
      surv.scale = "percent",
      break.time.by = break_time,
      xlab = "Follow-up (months)",
      ylab = y_lab,
      pval = pval_label,
      pval.coord = c(0, 0.15),
      risk.table = TRUE,
      legend.title = paste0(unlist(strsplit(gene_column, split='_', fixed=TRUE))[1], " level"),
      legend.labs = legend_order,
      font.legend = 10, 
      palette = "Dark2",
      surv.median.line = "hv",
      ggtheme = theme_light(),
      title = paste(type, origin)
    )
    # message(c('Médianes obtenues', type, ' ', origin, ' : \n', legend_order[1], ' : ', median(km_fit)[1], '\n', legend_order[2], ' : ', median(km_fit)[2]))
    # return(p1)
    msg <- c(
      paste('Médianes obtenues', type, origin, ': \n'),
      paste(legend_order[1], ':', median(km_fit)[1], '\n'),
      paste(legend_order[2], ':', median(km_fit)[2])
    )
    return(list(plot = p1, message = msg))
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(list(plot = NULL, message = e$message))
  })
}