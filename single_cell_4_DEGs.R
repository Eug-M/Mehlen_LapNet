# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

### from: https://ngs101.com/single-cell-rna-seq-beginners-guide-part-5-cell-type-specific-differential-expression-proportion-testing-and-functional-pathway-analysis/

#-----------------------------------------------
# STEP 1: Load required libraries
#-----------------------------------------------

# Core single-cell analysis
library(Seurat)
library(SeuratObject)

# Differential expression
library(DESeq2)
library(edgeR)
library(limma)

# Proportion analysis
library(speckle)

# Functional analysis
library(fgsea)
library(msigdbr)

# Visualization
library(ggplot2)
library(pheatmap)
library(ggridges)
library(ggrepel)
library(EnhancedVolcano)
library(patchwork)
library(dplyr)
library(tidyr)

# Set working directory
setwd("/home/eugenie-modolo/Documents/Gloriana/differential_analysis")

# Create output directories
dir.create("plots", showWarnings = FALSE)
dir.create("plots/proportions", showWarnings = FALSE)
dir.create("plots/DE_celllevel", showWarnings = FALSE)
dir.create("plots/DE_pseudobulk", showWarnings = FALSE)
dir.create("plots/functional_scoring", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Set random seed for reproducibility
set.seed(42)

# Configure plotting defaults
theme_set(theme_classic(base_size = 12))



#-----------------------------------------------
# STEP 2: Load annotated Seurat object from Part 4
#-----------------------------------------------

# Load the annotated object
seurat_obj <- readRDS("/home/eugenie-modolo/Documents/Gloriana/cell_type_annotation/annotations/integrated_annotated_seurat.rds")



#-----------------------------------------------
# STEP 3: Calculate cell type proportions per sample
#-----------------------------------------------

# Calculate proportions for each sample
proportion_data <- seurat_obj@meta.data %>%
  group_by(sample_id, condition, final_annotation) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(
    total_cells = sum(count),
    proportion = count / total_cells,
    percentage = proportion * 100
  ) %>%
  ungroup()

# Save proportion data
write.csv(proportion_data, "results/cell_type_proportions_per_sample.csv", row.names = FALSE)

# Display summary
proportion_summary <- proportion_data %>%
  group_by(condition, final_annotation) %>%
  summarise(
    mean_percentage = mean(percentage),
    sd_percentage = sd(percentage),
    .groups = "drop"
  )



#-----------------------------------------------
# STEP 4: Visualize cell type proportions
#-----------------------------------------------

# Stacked barplot showing composition of each sample
p_stacked <- ggplot(proportion_data, 
                    aes(x = sample_id, y = percentage, fill = final_annotation)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ condition, scales = "free_x") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(title = "Cell Type Composition by Sample",
       x = "Sample", y = "Percentage (%)",
       fill = "Cell Type") +
  scale_fill_brewer(palette = "Set3")

ggsave("plots/proportions/01_stacked_composition.png", p_stacked,
       width = 12, height = 6, dpi = 300)

# Grouped barplot with error bars
p_grouped <- ggplot(proportion_summary, 
                    aes(x = final_annotation, y = mean_percentage, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = mean_percentage - sd_percentage, ymax = mean_percentage + sd_percentage),
                position = position_dodge(0.9), width = 0.2) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mean Cell Type Proportions by Condition",
       x = "Cell Type", y = "Mean Percentage (%)") +
  scale_fill_manual(values = c("Control" = "#2E86AB", "NP137" = "#F18F01"))

ggsave("plots/proportions/02_grouped_proportions.png", p_grouped,
       width = 12, height = 6, dpi = 300)



#-----------------------------------------------
# STEP 5: Statistical testing of proportion differences
#-----------------------------------------------

# Prepare data for propeller
cluster_ids <- seurat_obj$final_annotation
sample_ids <- seurat_obj$sample_id

# Create group vector (condition per sample)
sample_to_condition <- unique(seurat_obj@meta.data[, c("sample_id", "condition")])
rownames(sample_to_condition) <- sample_to_condition$sample_id
group <- sample_to_condition[sample_ids, "condition"]

# Run propeller test
propeller_results <- propeller(
  clusters = cluster_ids,
  sample = sample_ids,
  group = group
)

# Save results
write.csv(propeller_results, "results/propeller_proportion_test.csv", row.names = FALSE)

# Conclusions : significant FDR for Pneumocyte type II (increase in NP137), fibroblast (decrease), 
# macrophage (decrease), lung club cell (increase), immune (decrease).



#-----------------------------------------------
# STEP 6: Cell-level DE analysis with FindMarkers
#-----------------------------------------------

# Get unique cell types
cell_types <- unique(seurat_obj$final_annotation)

# Set default assay and identify for comparison
DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- "condition"

# Function to run FindMarkers for one cell type
run_findmarkers <- function(celltype, seurat_obj) {
  cat("  Processing:", celltype, "\n")
  
  # Subset to this cell type
  seurat_obj_sub <- subset(seurat_obj, subset = final_annotation == celltype)
  
  # Run FindMarkers comparing conditions
  markers <- FindMarkers(
    seurat_obj_sub,
    ident.1 = "NP137",
    ident.2 = "Control",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.1,
    verbose = FALSE
  )
  
  # Add gene names and cell type
  markers$gene <- rownames(markers)
  markers$cell_type <- celltype
  
  return(markers)
}

# Run FindMarkers for all cell types
findmarkers_list <- lapply(cell_types, run_findmarkers, seurat_obj = seurat_obj)
names(findmarkers_list) <- cell_types

# Combine results
findmarkers_all <- do.call(rbind, findmarkers_list)
rownames(findmarkers_all) <- NULL

# Add significance flag
findmarkers_all$significant <- findmarkers_all$p_val_adj < 0.05 & 
  abs(findmarkers_all$avg_log2FC) > 0.25

# Save results
write.csv(findmarkers_all, "results/findmarkers_celllevel_results.csv", row.names = FALSE)

# Summary of significant genes per cell type
sig_summary <- findmarkers_all %>%
  filter(significant) %>%
  group_by(cell_type) %>%
  summarise(
    total_genes_tested = n(),
    sig_genes = sum(significant),
    sig_up = sum(avg_log2FC > 0.25 & p_val_adj < 0.05),
    sig_down = sum(avg_log2FC < -0.25 & p_val_adj < 0.05),
    .groups = "drop"
  )



#-----------------------------------------------
# STEP 7: Visualize FindMarkers results
#-----------------------------------------------

# Volcano plots for each cell type
for (ct in cell_types) {
  # Get data for this cell type
  ct_data <- findmarkers_all %>% filter(cell_type == ct)
  
  # Identify top genes to label
  top_genes <- ct_data %>%
    filter(significant) %>%
    arrange(p_val_adj) %>%
    head(10)
  
  # Create volcano plot
  p_volcano <- ggplot(ct_data, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = significant), alpha = 0.5, size = 1) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    theme_classic(base_size = 12) +
    labs(title = paste("Volcano Plot:", ct),
         x = "log2 Fold Change (NP137 vs Control)",
         y = "-log10(FDR)",
         color = "Significant") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "blue")
  
  # Add labels for top genes
  if (nrow(top_genes) > 0) {
    p_volcano <- p_volcano +
      geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 10)
  }
  
  ggsave(paste0("plots/DE_celllevel/volcano_", gsub(" ", "_", ct), ".png"),
         p_volcano, width = 8, height = 7, dpi = 300)
}



#-----------------------------------------------
# STEP 8: Create pseudobulk count matrices using PseudobulkExpression
#-----------------------------------------------

# Function to create pseudobulk for one cell type
create_pseudobulk_seurat5 <- function(celltype, seurat_obj) {
  # Subset to this cell type
  seurat_subset <- subset(seurat_obj, subset = final_annotation == celltype)
  
  # Use PseudobulkExpression to aggregate counts by sample
  pseudobulk_result <- PseudobulkExpression(
    seurat_subset,
    assays = "RNA",
    group.by = "sample_id",     # Aggregate by sample
    layer = "counts",            # Use raw counts
    method = "aggregate"         # Sum counts
  )
  
  # Extract the RNA counts matrix
  pseudobulk_matrix <- pseudobulk_result$RNA
  
  # Convert column names back to match original sample_id format
  colnames(pseudobulk_matrix) <- gsub("-", "_", colnames(pseudobulk_matrix))
  
  # Return as regular matrix (required for DESeq2)
  return(as.matrix(pseudobulk_matrix))
}

# Create pseudobulk for all cell types
pseudobulk_list <- lapply(cell_types, function(ct) {
  cat("  Processing:", ct, "\n")
  create_pseudobulk_seurat5(ct, seurat_obj)
})

names(pseudobulk_list) <- cell_types



#-----------------------------------------------
# STEP 9: Run DESeq2 pseudobulk analysis
#-----------------------------------------------

# Create sample metadata
sample_metadata <- unique(seurat_obj@meta.data[, c("sample_id", "condition", "patient_id")])

sample_metadata <- as.data.frame(sample_metadata)
rownames(sample_metadata) <- sample_metadata$sample_id

# Function to run DESeq2 for one cell type
run_deseq2_celltype <- function(celltype, pseudobulk_counts, sample_metadata) {
  
  # Get sample names from count matrix
  samples_in_counts <- colnames(pseudobulk_counts)
  
  # Verify that sample_metadata rownames match count matrix colnames
  if (!all(samples_in_counts %in% rownames(sample_metadata))) {
    stop("ERROR: Sample names in count matrix don't match sample_metadata rownames!")
  }
  
  # Order metadata to match count matrix columns
  sample_metadata_ordered <- sample_metadata[samples_in_counts, , drop = FALSE]
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = pseudobulk_counts,
    colData = sample_metadata_ordered,
    design = ~ condition
  )
  
  # Filter low-count genes (improve power and speed)
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  
  # Run DESeq2 analysis
  dds <- DESeq(dds, quiet = TRUE)
  
  # Extract results
  res <- results(
    dds,
    contrast = c("condition", "NP137", "Control"),
    alpha = 0.05
  )
  
  # Convert to data frame
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$cell_type <- celltype
  
  # Add significance flag
  res_df$significant <- !is.na(res_df$padj) & 
    res_df$padj < 0.05 & 
    abs(res_df$log2FoldChange) > 0.25
  
  return(res_df)
}

# Run DESeq2 for all cell types
deseq2_list <- lapply(cell_types, function(ct) {
  run_deseq2_celltype(ct, pseudobulk_list[[ct]], sample_metadata)
})

names(deseq2_list) <- cell_types

# Combine results
deseq2_all <- do.call(rbind, deseq2_list)
rownames(deseq2_all) <- NULL

# Save results
write.csv(deseq2_all, "results/deseq2_pseudobulk_results.csv", row.names = FALSE)

# Summary of significant genes per cell type
deseq2_summary <- deseq2_all %>%
  group_by(cell_type) %>%
  summarise(
    total_genes_tested = n(),
    sig_genes = sum(significant, na.rm = TRUE),
    sig_up = sum(significant & log2FoldChange > 0.25, na.rm = TRUE),
    sig_down = sum(significant & log2FoldChange < -0.25, na.rm = TRUE),
    .groups = "drop"
  )



#-----------------------------------------------
# STEP 10: limma-voom pseudobulk analysis
#-----------------------------------------------

# Function to run limma-voom for one cell type
limma_paired_analysis <- function(
    counts,
    meta,
    design_formula = ~0 + condition,
    contrast_str = "NP137 - Control",
    min_count = 10,
    voom_with_quality_weights = TRUE
) {
  
  # Order metadata to match count matrix columns
  meta <- meta[match(colnames(counts), meta$sample_id), ]
  rownames(meta) <- meta$sample_id
  
  # Convert factors
  meta$condition <- factor(meta$condition, levels = c("Control", "NP137"))
  
  # Filter low expressed genes
  keep <- rowSums(counts) >= min_count
  counts <- counts[keep, ]
  
  # Create design matrix
  design <- model.matrix(design_formula, data = meta)
  colnames(design) <- gsub("condition", "", colnames(design))
  
  # Voom transformation with quality weights
  if (voom_with_quality_weights) {
    v <- voomWithQualityWeights(counts, design, plot = FALSE)
  } else {
    v <- voom(counts, design, plot = FALSE)
  }
  
  # Fit linear model
  fit <- lmFit(v, design)
  
  # Define contrast
  contrast.matrix <- makeContrasts(
    contrasts = contrast_str,
    levels = colnames(design)
  )
  
  # Fit contrasts
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  
  # Get results
  res <- topTable(fit2, number = Inf, sort.by = "none", adjust.method = "BH")
  
  return(res)
}

# Run limma for all cell types
limma_list <- lapply(cell_types, function(ct) {
  cat("\nProcessing:", ct, "\n")
  limma_paired_analysis(pseudobulk_list[[ct]], sample_metadata)
})

names(limma_list) <- cell_types

# Process limma results
limma_all <- lapply(names(limma_list), function(ct) {
  res <- limma_list[[ct]]
  res$gene <- rownames(res)
  res$cell_type <- ct
  res$significant <- res$adj.P.Val < 0.05 & abs(res$logFC) > 0.25
  return(res)
}) %>% bind_rows()

# Save results
write.csv(limma_all, "results/limma_pseudobulk_results.csv", row.names = FALSE)

# Summary of significant genes per cell type
limma_summary <- limma_all %>%
  group_by(cell_type) %>%
  summarise(
    total_genes_tested = n(),
    sig_genes = sum(significant, na.rm = TRUE),
    sig_up = sum(significant & logFC > 0.25, na.rm = TRUE),
    sig_down = sum(significant & logFC < -0.25, na.rm = TRUE),
    .groups = "drop"
  )



#-----------------------------------------------
# STEP 11: Download gene sets from MSigDB
#-----------------------------------------------

# Get Hallmark gene sets for humans
hallmark_gene_sets <- msigdbr(
  species = "Mus musculus",
  collection = "H"  # Hallmark gene sets
)

# Create named list format for scoring
# We'll use HALLMARK_INFLAMMATORY_RESPONSE 
hallmark_list <- hallmark_gene_sets %>%
  filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>%
  pull(gene_symbol)



#-----------------------------------------------
# STEP 12: Score cells for inflammatory response pathway
#-----------------------------------------------

# Use Seurat's AddModuleScore to calculate pathway scores
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = list(hallmark_list),
  name = "Inflammatory_Response",
  nbin = 24,
  ctrl = 100
)

# The score is added as a new metadata column
# Rename for clarity (Seurat adds "1" to the name)
colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == "Inflammatory_Response1"] <- "Inflammatory_Response_Score"



#-----------------------------------------------
# STEP 13: Visualize inflammatory response scores
#-----------------------------------------------

# Ridge plot showing score distributions by cell type and condition
p_ridge <- ggplot(seurat_obj@meta.data, 
                  aes(x = Inflammatory_Response_Score, 
                      y = final_annotation, 
                      fill = condition)) +
  geom_density_ridges(alpha = 0.7) +
  theme_classic(base_size = 12) +
  labs(title = "Inflammatory Response Pathway Scores",
       x = "Pathway Score",
       y = "Cell Type",
       fill = "Condition") +
  scale_fill_manual(values = c("Control" = "#2E86AB", "NP137" = "#F18F01"))

ggsave("plots/functional_scoring/ridge_inflammatory_response.png", p_ridge,
       width = 10, height = 8, dpi = 300)



#-----------------------------------------------
# STEP 14: Pseudobulk statistical testing of pathway scores
#-----------------------------------------------

# Aggregate scores by sample and cell type
score_aggregated <- seurat_obj@meta.data %>%
  group_by(sample_id, condition, final_annotation) %>%
  summarise(
    mean_score = mean(Inflammatory_Response_Score),
    .groups = "drop"
  )

# Function to test one cell type
test_pathway_celltype <- function(celltype, score_data) {
  # Filter to this cell type
  ct_data <- score_data %>% filter(final_annotation == celltype)
  
  if (nrow(ct_data) < 6) {
    return(data.frame(
      cell_type = celltype,
      mean_healthy = NA,
      mean_post = NA,
      diff = NA,
      p_value = NA,
      significant = FALSE
    ))
  }
  
  # Separate by condition
  healthy_scores <- ct_data %>% filter(condition == "Control") %>% pull(mean_score)
  post_scores <- ct_data %>% filter(condition == "NP137") %>% pull(mean_score)
  
  # Statistical test (t-test)
  test_result <- t.test(post_scores, healthy_scores)
  
  return(data.frame(
    cell_type = celltype,
    mean_healthy = mean(healthy_scores),
    mean_post = mean(post_scores),
    diff = mean(post_scores) - mean(healthy_scores),
    p_value = test_result$p.value
  ))
}

# Test all cell types
pathway_tests <- lapply(cell_types, test_pathway_celltype, 
                        score_data = score_aggregated)
pathway_results <- do.call(rbind, pathway_tests)

# Add FDR correction
pathway_results$FDR <- p.adjust(pathway_results$p_value, method = "BH")

# Save results
write.csv(pathway_results, "results/inflammatory_response_pathway_test.csv", 
          row.names = FALSE)
# Positive diff: Pathway score higher in NP137 (pathway activated)
# Conclusions: pathway significantly suppressed in endothelium, fibroblast, epithelial, lung club cell, Pneumocyte type II

