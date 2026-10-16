# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

### from: https://ngs101.com/how-to-analyze-single-cell-rna-seq-data-complete-beginners-guide-part-4-cell-type-identification/

# SingleR + celldex: Reference-based annotation with curated cell type atlases
# scCATCH: Tissue-specific marker databases for automated annotation
# scType dependencies: Tools for gene symbol validation and Excel file reading
# Visualization packages: For creating comparison plots and Sankey diagrams

#-----------------------------------------------
# STEP 1: Load required libraries
#-----------------------------------------------

# Core single-cell analysis
library(Seurat)
library(SeuratObject)

# Automated annotation methods
library(SingleR)
library(celldex)
library(scCATCH)

# Data structures
library(SingleCellExperiment)

# Marker processing
library(HGNChelper)

# Visualization and data manipulation
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(patchwork)
library(scales)
library(viridis)

# Set working directory
setwd("/home/eugenie-modolo/Documents/Gloriana/cell_type_annotation")

# Create output directories
dir.create("plots", showWarnings = FALSE)
dir.create("plots/manual_annotation", showWarnings = FALSE)
dir.create("plots/automated_annotation", showWarnings = FALSE)
dir.create("plots/method_comparison", showWarnings = FALSE)
dir.create("annotations", showWarnings = FALSE)

# Set random seed for reproducibility
set.seed(42)

# Configure plotting defaults
theme_set(theme_classic(base_size = 12))



#-----------------------------------------------
# STEP 2: Load integrated data from previous step
#-----------------------------------------------

# Path to integrated Seurat object from previous step
integrated_path <- "/home/eugenie-modolo/Documents/Gloriana/single_cell_integration_analysis/integrated_data/integrated_clustered_seurat.rds"

# Load the object
integrated_seurat <- readRDS(integrated_path)

# Determine which UMAP to use (from best integration method in previous step)
available_reductions <- names(integrated_seurat@reductions)
if ("umap.harmony" %in% available_reductions) {
  umap_reduction <- "umap.harmony"
} else if ("umap.cca" %in% available_reductions) {
  umap_reduction <- "umap.cca"
} else if ("umap.rpca" %in% available_reductions) {
  umap_reduction <- "umap.rpca"
} else if ("umap.mnn" %in% available_reductions) {
  umap_reduction <- "umap.mnn"
} else {
  umap_reduction <- "umap"
}



#-----------------------------------------------
# STEP 3: Visualize clusters before annotation
#-----------------------------------------------

# Create multi-panel overview
p1 <- DimPlot(integrated_seurat, reduction = umap_reduction,
              group.by = "seurat_clusters", label = TRUE, label.size = 5,
              pt.size = 0.1, raster=FALSE) +
  ggtitle("Clusters (Pre-Annotation)") +
  theme(plot.title = element_text(face = "bold", size = 14))

p2 <- DimPlot(integrated_seurat, reduction = umap_reduction,
              group.by = "sample_id", pt.size = 0.1, raster=FALSE) +
  ggtitle("Samples") +
  theme(legend.text = element_text(size = 7))

p3 <- DimPlot(integrated_seurat, reduction = umap_reduction,
              group.by = "condition", pt.size = 0.1, raster=FALSE) +
  ggtitle("Condition") +
  scale_color_manual(values = c("Control" = "#2E86AB", 
                                "NP137" = "#F18F01"))

# Combine
p_overview <- (p1 | p2) / (p3 | plot_spacer())

ggsave("plots/00_starting_clusters.png", p_overview, 
       width = 14, height = 10, dpi = 300)



#-----------------------------------------------
# STEP 4: Normalize and scale data (Critical: FindAllMarkers, FeaturePlot, DotPlot, and DoHeatmap all require normalized data)
#-----------------------------------------------

# Normalize data (log-normalization)
integrated_seurat <- NormalizeData(integrated_seurat, 
                                   normalization.method = "LogNormalize",
                                   scale.factor = 10000, 
                                   verbose = FALSE)

# Scale all genes (required for DoHeatmap and marker analysis)
integrated_seurat <- ScaleData(integrated_seurat, 
                               features = rownames(integrated_seurat), 
                               verbose = FALSE)



#-----------------------------------------------
# STEP 5: Define and check canonical markers
#-----------------------------------------------

# Define comprehensive marker panel for PBMC cell types
canonical_markers <- list(
  "immune_cells" = c("Ptprc", "Cd86", "Cd4", "Cd8a", "Cd28", "Cd68"),
  "lymphocytes" = c("Cd247", "Cd3g", "Cd3d", "Cd3e", "Cd2"),
  "dendrits_and_monocytes" = c("Itgax", "Cd14", "Mafb"),
  "Treg_lymphocytes" = c("Cd2", "Foxp3", "Cdh5", "Cd8a"),
  "activated_lymphocytes" = c("Gzmb", "Gzma"),
  "granzymes" = c("Fasl", "Tgfbr1", "Itgb8", "Cd4"),
  "Epithelial" = c("Epcam"),
  "Pericyte" = c("Acta2", "Cspg4"),
  "Endothelium" = c("Cdh5", "Pecam1"),
  "Endothelium_lymphatic" = c("Prox1", "Lyve1"),
  # "" = c(),
  "T_cells" = c("Cd3d", "Cd3e", "Cd3g"),
  "CD4_T" = c("Cd3d", "Cd4", "Il7r"),
  "CD8_T" = c("Cd3d", "Cd8a", "Cd8b1"),
  "B_cells" = c("Cd79a", "Ms4a1", "Cd19"),
  "NK_cells" = c("Nkg7", "Gnly", "Ncam1"),
  "Monocytes" = c("Cd14", "Lyz3", "S100a8", "S100a9"),
  "Classical_Mono" = c("Cd14", "S100a8", "Fcgr4"),
  "Non_Classical_Mono" = c("Fcgr4", "Ms4a7"),
  # "Dendritic_cells" = c("Fcer1a", "Cst3"),
  # "Plasmacytoid_DC" = c("Il3ra", "Gzmb", "Serpinf1"),
  "Platelets" = c("Ppbp", "Pf4", "Gp9"),
  "Mesenchymal" = c("Vim"),
  "Alveolar Fibroblast" = c("Pdgfra", "Npnt", "Ces1d", "Tcf21", "Scube2"),
  "Lipofibroblast" = c("Plin2", "Tcf21"),
  "Myofibroblast" = c("Acta2", "Postn", "Cthrc1", "Mmp11", "Col1a1", "Col3a1")
)

# Check which markers are present in dataset
all_genes <- rownames(integrated_seurat)

for (cell_type in names(canonical_markers)) {
  markers <- canonical_markers[[cell_type]]
  present <- markers %in% all_genes
  
  cat(sprintf("%-20s: %d/%d markers present\n", 
              cell_type, sum(present), length(markers)))
  
  if (!all(present)) {
    missing <- markers[!present]
    cat(sprintf("  Missing: %s\n", paste(missing, collapse = ", ")))
  }
}



#-----------------------------------------------
# STEP 6: Visualize canonical markers
#-----------------------------------------------

# Function to create violin plots for a marker set
plot_violin_markers <- function(markers, cell_type_name, filename_suffix) {
  present_markers <- markers[markers %in% rownames(integrated_seurat)]
  
  if (length(present_markers) == 0) {
    cat("  ⚠ No markers present for", cell_type_name, "\n")
    return(NULL)
  }
  
  cat("  •", cell_type_name, ":", length(present_markers), "markers\n")
  
  p <- VlnPlot(integrated_seurat,
               features = present_markers,
               group.by = "seurat_clusters", split.by = "condition",
               pt.size = 0,
               ncol = 3, raster=FALSE) &
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 11, face = "bold"))
  
  n_rows <- ceiling(length(present_markers) / 3)
  height <- max(4, n_rows * 2.5)
  
  ggsave(paste0("plots/manual_annotation/violin_", filename_suffix, ".png"),
         p, width = 14, height = height, dpi = 300)
  
  return(p)
}

# Function to create UMAP feature plots for a marker set
plot_umap_markers <- function(markers, cell_type_name, filename_suffix) {
  present_markers <- markers[markers %in% rownames(integrated_seurat)]
  
  if (length(present_markers) == 0) {
    return(NULL)
  }
  
  p <- FeaturePlot(integrated_seurat,
                   features = present_markers,
                   reduction = umap_reduction,
                   ncol = 4,
                   pt.size = 0.1, raster=FALSE) &
    theme(plot.title = element_text(size = 10),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right")
  
  n_rows <- ceiling(length(present_markers) / 4)
  height <- max(4, n_rows * 3)
  
  ggsave(paste0("plots/manual_annotation/umap_", filename_suffix, ".png"),
         p, width = 16, height = height, dpi = 300)
  
  return(p)
}

# Generate plots for ALL cell types
for (cell_type_key in names(canonical_markers)) {
  markers <- canonical_markers[[cell_type_key]]
  cell_type_name <- gsub("_", " ", cell_type_key)
  
  plot_violin_markers(markers, cell_type_name, cell_type_key)
  plot_umap_markers(markers, cell_type_name, cell_type_key)
}

# Create comprehensive DotPlot showing all key markers
all_markers <- unique(unlist(canonical_markers))
all_markers_present <- all_markers[all_markers %in% rownames(integrated_seurat)]

p_dotplot <- DotPlot(integrated_seurat,
                     features = all_markers_present,
                     group.by = "seurat_clusters") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 8)) +
  labs(title = "All Canonical Markers by Cluster",
       subtitle = paste(length(all_markers_present), "markers across all cell types"))

ggsave("plots/manual_annotation/dotplot_all_canonical_markers.png", p_dotplot,
       width = 12, height = max(8, length(all_markers_present) * 0.15), dpi = 300)



#-----------------------------------------------
# STEP 7: Find cluster-specific marker genes
#-----------------------------------------------

# Find markers for each cluster vs all other cells
cluster_markers <- FindAllMarkers(
  integrated_seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

# Remove ribosomal and mitochondrial genes
cluster_markers_filtered <- cluster_markers %>%
  dplyr::filter(!grepl("^Rp[sl]", gene)) %>%
  dplyr::filter(!grepl("^mt-", gene))

# Get top 5 markers per cluster
top_markers <- cluster_markers_filtered %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  arrange(cluster, dplyr::desc(avg_log2FC))

cat("\nTop 5 markers per cluster:\n")
print(top_markers %>% dplyr::select(cluster, gene, avg_log2FC, pct.1, pct.2))

# Save all markers
write.csv(cluster_markers_filtered, 
          "annotations/cluster_markers_all.csv",
          row.names = FALSE)

write.csv(top_markers,
          "annotations/top5_markers_per_cluster.csv",
          row.names = FALSE)

# Create heatmap of top markers
top_genes <- top_markers$gene

p_heatmap <- DoHeatmap(
  integrated_seurat,
  features = top_genes,
  group.by = "seurat_clusters",
  size = 3
) +
  theme(axis.text.y = element_text(size = 6)) +
  labs(title = "Top 5 Markers per Cluster")

ggsave("plots/manual_annotation/06_heatmap_top_markers.png", p_heatmap,
       width = 12, height = 14, dpi = 300)



#-----------------------------------------------
# STEP 8: Manual cell type assignment
#-----------------------------------------------

# MAPPING - CUSTOMIZED using the figures from Step 7
cluster_to_celltype <- c(
  "0" = "endothelium",
  "1" = "epithelial",
  "2" = "epithelial",
  "3" = "lung activated fibroblast",
  "4" = "epithelial",
  "5" = "immune",
  "6" = "lung alveolar fibroblast",
  "7" = "B cells",
  "8" = "epithelial",
  "9" = "T cells",
  "10" = "monocyte",
  "11" = "endothelium",
  "12" = "epithelial",
  "13" = "endothelium",
  "14" = "immune",
  "15" = "epithelial",
  "16" = "endothelium",
  "17" = "lung activated fibroblast & pericyte",
  "18" = "immune",
  "19" = "immune",
  "20" = "immune",
  "21" = "epithelial",
  "22" = "endothelium"
)

# Apply to all cells
cluster_ids <- as.character(integrated_seurat$seurat_clusters)
integrated_seurat$manual_annotation <- unname(cluster_to_celltype[cluster_ids])

# Visualize manual annotation
p_manual <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "manual_annotation",
  label = TRUE,
  label.size = 4,
  pt.size = 0.1,
  repel = TRUE,
  raster=FALSE
) +
  ggtitle("Manual Cell Type Annotation") +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("plots/manual_annotation/07_manual_annotation_umap.png", p_manual,
       width = 10, height = 8, dpi = 300)

# Split by condition
p_manual_split <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "manual_annotation",
  split.by = "condition",
  pt.size = 0.1,
  ncol = 2,
  raster=FALSE
) +
  ggtitle("Manual Annotation by Condition")

ggsave("plots/manual_annotation/08_manual_annotation_by_condition.png",
       p_manual_split, width = 14, height = 6, dpi = 300)



#-----------------------------------------------
# STEP 9: SingleR with multiple reference datasets
#-----------------------------------------------

# Convert Seurat to SingleCellExperiment for SingleR
sce <- as.SingleCellExperiment(integrated_seurat)

# Reference 1: MouseRNAseqData (Broad)
mouse_ref <- celldex::MouseRNAseqData()

singler_mouse <- SingleR(
  test = sce,
  ref = mouse_ref,
  labels = mouse_ref$label.main,
  assay.type.test = "logcounts"
)

integrated_seurat$singler_mouse <- singler_mouse$labels

# Reference 2: Immgen (Mouse microarray)
immgen_ref <- celldex::ImmGenData()

singler_immgen <- SingleR(
  test = sce,
  ref = immgen_ref,
  labels = immgen_ref$label.main,
  assay.type.test = "logcounts"
)

integrated_seurat$singler_immgen <- singler_immgen$labels


#-----------------------------------------------
# STEP 10: Compare SingleR references
#-----------------------------------------------

# Visualize all three references
p_mouse <- DimPlot(integrated_seurat, reduction = umap_reduction,
                  group.by = "singler_mouse", pt.size = 0.1, raster=FALSE) +
  labs(title = "Mouse RNASeq") +
  theme(legend.text = element_text(size = 7))

p_immgen <- DimPlot(integrated_seurat, reduction = umap_reduction,
                    group.by = "singler_immgen", pt.size = 0.1, raster=FALSE) +
  labs(title = "Mouse microarray") +
  theme(legend.text = element_text(size = 7))

p_manual_ref <- DimPlot(integrated_seurat, reduction = umap_reduction,
                        group.by = "manual_annotation", pt.size = 0.1, raster=FALSE) +
  labs(title = "Manual") +
  theme(legend.text = element_text(size = 7))

p_ref_comparison <- (p_manual_ref | p_mouse) / (p_immgen)

ggsave("plots/automated_annotation/10_singler_reference_comparison.png",
       p_ref_comparison, width = 16, height = 12, dpi = 300)

# Choosing the annotation /!\ à adapter !
integrated_seurat$singler_annotation <- integrated_seurat$singler_mouse

# Visualize chosen reference
p_singler_final <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "singler_annotation",
  label = TRUE,
  label.size = 3,
  pt.size = 0.1,
  repel = TRUE,
  raster=FALSE
) +
  ggtitle("SingleR Annotation (Manual)") +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("plots/automated_annotation/11_singler_final_annotation.png",
       p_singler_final, width = 11, height = 8, dpi = 300)



#-----------------------------------------------
# STEP 11: Load scType functions
#-----------------------------------------------

# Load scType functions from GitHub
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R") # not useful when using prepare_gene_sets_mouse function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")



#-----------------------------------------------
# STEP 12: Load scType marker database
#-----------------------------------------------


# Run scripts from file "scType_mouse_lung.R" !

# open both db files
library(openxlsx)
db_sctype_wolung <- read.xlsx("ScTypeDB_mouse_wolung.xlsx")
db_sctype_withlung <- read.xlsx("ScTypeDB_mouse_withlung.xlsx")
db_cellmarker <- read.xlsx("cellmarker_sctype_format.xlsx")

# Prepare gene sets WITHOUT using gene_sets_prepare() from scType  (keeps mouse case, no HGNC check)
prepare_gene_sets_mouse <- function(db, tissue) {
  d <- db[db$tissueType %in% tissue, ]
  parse <- function(x) lapply(x, function(s) {
    if (is.na(s) || s == "") character(0) else unique(trimws(strsplit(s, ",")[[1]]))
  })
  gs  <- parse(d$geneSymbolmore1); names(gs)  <- d$cellName
  gs2 <- parse(d$geneSymbolmore2); names(gs2) <- d$cellName
  list(gs_positive = gs, gs_negative = gs2)
}

gs_list_sctype_wolung <- prepare_gene_sets_mouse(db_sctype_wolung, unique(db_sctype_wolung$tissueType))
gs_list_sctype_withlung <- prepare_gene_sets_mouse(db_sctype_withlung, unique(db_sctype_withlung$tissueType))
gs_list_cellmarker <- prepare_gene_sets_mouse(db_cellmarker, unique(db_cellmarker$tissueType))

# checking all went well
gs_list_sctype_withlung$gs_positive[["Alveolar type 2 cells (AT2)"]]   # should be Sftpc, Sftpb, ...

cat("Available cell types in database sctype_wolung:\n")
print(names(gs_list_sctype_wolung$gs_positive))
cat("Available cell types in database sctype_withlung:\n")
print(names(gs_list_sctype_withlung$gs_positive))
cat("Available cell types in database cellmarker:\n")
print(names(gs_list_cellmarker$gs_positive))



#-----------------------------------------------
# STEP 13: Calculate scType scores
#-----------------------------------------------

# Get scaled data
scaled_data <- as.matrix(LayerData(integrated_seurat, layer = "scale.data"))

## sctype_wolung
# Calculate scores
es_max_sctype_wolung <- sctype_score(
  scRNAseqData = scaled_data,
  scaled = TRUE,
  gs = gs_list_sctype_wolung$gs_positive,
  gs2 = gs_list_sctype_wolung$gs_negative,
  gene_names_to_uppercase = FALSE # to use only if using the prepare_gene_sets_mouse function
)
# Assign cell types to clusters
clusters <- integrated_seurat$seurat_clusters
cL_results_sctype_wolung <- do.call("rbind", lapply(unique(clusters), function(cl) {
  es_max_cl <- sort(rowSums(es_max_sctype_wolung[, clusters == cl]), decreasing = TRUE)
  
  top_score <- es_max_cl[1]
  top_type <- names(es_max_cl)[1]
  
  data.frame(
    cluster = cl,
    type = top_type,
    scores = top_score,
    ncells = sum(clusters == cl)
  )
}))
# Filter low-confidence assignments
score_threshold_sctype_wolung <- quantile(cL_results_sctype_wolung$scores, 0.6)
cL_results_sctype_wolung$type[cL_results_sctype_wolung$scores < score_threshold_sctype_wolung] <- "Unknown"
# Create mapping and apply to cells
cluster_to_sctype_sctype_wolung <- setNames(cL_results_sctype_wolung$type, cL_results_sctype_wolung$cluster)
integrated_seurat$sctype_wolung_annotation <- unname(cluster_to_sctype_sctype_wolung[
  as.character(integrated_seurat$seurat_clusters)
])
# Visualize
p_sctype_sctype_wolung <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "sctype_wolung_annotation",
  label = TRUE,
  label.size = 4,
  pt.size = 0.1,
  repel = TRUE,
  raster=FALSE
) +
  ggtitle("scType: wolung") +
  theme(plot.title = element_text(face = "bold", size = 14))
ggsave("plots/automated_annotation/12_sctype_annotation_sctype_wolung.png", p_sctype_sctype_wolung,
       width = 11, height = 8, dpi = 300)
# rm(es_max_sctype_wolung)
# rm(p_sctype_sctype_wolung)


## sctype_withlung
es_max_sctype_withlung <- sctype_score(
  scRNAseqData = scaled_data,
  scaled = TRUE,
  gs = gs_list_sctype_withlung$gs_positive,
  gs2 = gs_list_sctype_withlung$gs_negative,
  gene_names_to_uppercase = FALSE # to use only if using the prepare_gene_sets_mouse function
)
cL_results_sctype_withlung <- do.call("rbind", lapply(unique(clusters), function(cl) {
  es_max_cl <- sort(rowSums(es_max_sctype_withlung[, clusters == cl]), decreasing = TRUE)
  
  top_score <- es_max_cl[1]
  top_type <- names(es_max_cl)[1]
  
  data.frame(
    cluster = cl,
    type = top_type,
    scores = top_score,
    ncells = sum(clusters == cl)
  )
}))
score_threshold_sctype_withlung <- quantile(cL_results_sctype_withlung$scores, 0.6)
cL_results_sctype_withlung$type[cL_results_sctype_withlung$scores < score_threshold_sctype_withlung] <- "Unknown"
cluster_to_sctype_sctype_withlung <- setNames(cL_results_sctype_withlung$type, cL_results_sctype_withlung$cluster)
integrated_seurat$sctype_withlung_annotation <- unname(cluster_to_sctype_sctype_withlung[
  as.character(integrated_seurat$seurat_clusters)
])
p_sctype_sctype_withlung <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "sctype_withlung_annotation",
  label = TRUE,
  label.size = 4,
  pt.size = 0.1,
  repel = TRUE,
  raster=FALSE
) +
  ggtitle("scType: withlung") +
  theme(plot.title = element_text(face = "bold", size = 14))
ggsave("plots/automated_annotation/12_sctype_annotation_sctype_withlung.png", p_sctype_sctype_withlung,
       width = 11, height = 8, dpi = 300)
# rm(es_max_sctype_withlung)
# rm(p_sctype_sctype_withlung)


## cellmarker
es_max_cellmarker <- sctype_score(
  scRNAseqData = scaled_data,
  scaled = TRUE,
  gs = gs_list_cellmarker$gs_positive,
  gs2 = gs_list_cellmarker$gs_negative,
  gene_names_to_uppercase = FALSE # to use only if using the prepare_gene_sets_mouse function
)
cL_results_cellmarker <- do.call("rbind", lapply(unique(clusters), function(cl) {
  es_max_cl <- sort(rowSums(es_max_cellmarker[, clusters == cl]), decreasing = TRUE)
  
  top_score <- es_max_cl[1]
  top_type <- names(es_max_cl)[1]
  
  data.frame(
    cluster = cl,
    type = top_type,
    scores = top_score,
    ncells = sum(clusters == cl)
  )
}))
score_threshold_cellmarker <- quantile(cL_results_cellmarker$scores, 0.6)
cL_results_cellmarker$type[cL_results_cellmarker$scores < score_threshold_cellmarker] <- "Unknown"
cluster_to_sctype_cellmarker <- setNames(cL_results_cellmarker$type, cL_results_cellmarker$cluster)
integrated_seurat$cellmarker_annotation <- unname(cluster_to_sctype_cellmarker[
  as.character(integrated_seurat$seurat_clusters)
])
p_sctype_cellmarker <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "cellmarker_annotation",
  label = TRUE,
  label.size = 4,
  pt.size = 0.1,
  repel = TRUE,
  raster=FALSE
) +
  ggtitle("scType: from CellMarker2.0") +
  theme(plot.title = element_text(face = "bold", size = 14))
ggsave("plots/automated_annotation/12_sctype_annotation_cellmarker.png", p_sctype_cellmarker,
       width = 11, height = 8, dpi = 300)
# rm(es_max_cellmarker)
# rm(p_sctype_cellmarker)


# Save scType results (! à modifier en indiquant le clustering choisi)
integrated_seurat$sctype_annotation <- integrated_seurat$sctype_wolung_annotation
write.csv(cL_results_sctype_wolung, "annotations/sctype_wolung_cluster_scores.csv", row.names = FALSE)



#-----------------------------------------------
# STEP 14: scCATCH annotation
#-----------------------------------------------

# Prepare data for scCATCH
obj_sccatch <- createscCATCH(
  data = LayerData(integrated_seurat, layer = "data"),
  cluster = as.character(integrated_seurat$seurat_clusters)
)

# Find marker genes for each cluster 
obj_sccatch <- findmarkergene(
  object = obj_sccatch,
  species = "Mouse",
  marker = cellmatch,
  tissue = "Lung", # from: https://github.com/ZJUFanLab/scCATCH/wiki/mouse_tissues
  cancer = "Normal", 
  cell_min_pct = 0.25,
  logfc = 0.25,
  pvalue = 0.05
)

# Find cell types based on markers
obj_sccatch <- findcelltype(obj_sccatch)

# Extract results
sccatch_celltype <- obj_sccatch@celltype

# Create cluster to cell type mapping
cluster_to_sccatch <- setNames(
  sccatch_celltype$cell_type,
  sccatch_celltype$cluster
)

# Apply to all cells
integrated_seurat$sccatch_annotation <- unname(cluster_to_sccatch[
  as.character(integrated_seurat$seurat_clusters)
])

# Handle unmapped clusters
integrated_seurat$sccatch_annotation[
  is.na(integrated_seurat$sccatch_annotation)
] <- "Unknown"

# Visualize
p_sccatch <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "sccatch_annotation",
  label = TRUE,
  label.size = 4,
  pt.size = 0.1,
  repel = TRUE,
  raster=FALSE
) +
  ggtitle("scCATCH: Tissue-Specific Database") +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("plots/automated_annotation/13_sccatch_annotation.png", p_sccatch,
       width = 11, height = 8, dpi = 300)

# Save detailed results
write.csv(sccatch_celltype, 
          "annotations/sccatch_cluster_assignments.csv",
          row.names = FALSE)

# Save marker evidence
marker_evidence <- obj_sccatch@markergene
write.csv(marker_evidence,
          "annotations/sccatch_marker_evidence.csv",
          row.names = FALSE)



#-----------------------------------------------
# STEP 15: Generate confusion matrices 
#-----------------------------------------------

# Function to create confusion matrix
create_confusion_matrix <- function(method1_labels, method2_labels,
                                    method1_name, method2_name) {
  conf_table <- table(
    Method1 = method1_labels,
    Method2 = method2_labels
  )
}

# Compare each automated method to manual
conf_manual_vs_singler <- create_confusion_matrix(
  integrated_seurat$manual_annotation,
  integrated_seurat$singler_annotation,
  "Manual", "SingleR"
)

conf_manual_vs_sctype <- create_confusion_matrix(
  integrated_seurat$manual_annotation,
  integrated_seurat$sctype_annotation,
  "Manual", "scType"
)

conf_manual_vs_sccatch <- create_confusion_matrix(
  integrated_seurat$manual_annotation,
  integrated_seurat$sccatch_annotation,
  "Manual", "scCATCH"
)

# Save confusion matrices
write.csv(conf_manual_vs_singler,
          "annotations/confusion_matrix_manual_vs_singler.csv")
write.csv(conf_manual_vs_sctype,
          "annotations/confusion_matrix_manual_vs_sctype.csv")
write.csv(conf_manual_vs_sccatch,
          "annotations/confusion_matrix_manual_vs_sccatch.csv")



#-----------------------------------------------
# STEP 16: Create Sankey diagrams 
#-----------------------------------------------

# Function to prepare Sankey data
prepare_sankey_data <- function(method1_labels, method2_labels,
                                method1_name, method2_name) {
  sankey_df <- data.frame(
    Method1 = method1_labels,
    Method2 = method2_labels
  ) %>%
    group_by(Method1, Method2) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(Comparison = paste(method1_name, "vs", method2_name))

  return(sankey_df)
}

# Prepare data for each comparison
sankey_manual_singler <- prepare_sankey_data(
  integrated_seurat$manual_annotation,
  integrated_seurat$singler_annotation,
  "Manual", "SingleR"
)

sankey_manual_sctype <- prepare_sankey_data(
  integrated_seurat$manual_annotation,
  integrated_seurat$sctype_annotation,
  "Manual", "scType"
)

sankey_manual_sccatch <- prepare_sankey_data(
  integrated_seurat$manual_annotation,
  integrated_seurat$sccatch_annotation,
  "Manual", "scCATCH"
)

# Function to create Sankey plot
create_sankey_plot <- function(sankey_data, method1_name, method2_name) {
  p <- ggplot(sankey_data,
              aes(axis1 = Method1, axis2 = Method2, y = Count)) +
    geom_alluvium(aes(fill = Method1), width = 1/12, alpha = 0.7) +
    geom_stratum(width = 1/12, fill = "white", color = "grey") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              size = 3, fontface = "bold") +
    scale_x_discrete(limits = c(method1_name, method2_name),
                     expand = c(0.05, 0.05)) +
    labs(title = paste(method1_name, "vs", method2_name),
         subtitle = "Straight flows = agreement | Crossed flows = disagreement",
         y = "Number of Cells") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(face = "bold"))

  return(p)
}

# Create Sankey plots
p_sankey_singler <- create_sankey_plot(sankey_manual_singler, "Manual", "SingleR")
ggsave("plots/method_comparison/14_sankey_manual_vs_singler.png",
       p_sankey_singler, width = 12, height = 10, dpi = 300)

p_sankey_sctype <- create_sankey_plot(sankey_manual_sctype, "Manual", "scType")
ggsave("plots/method_comparison/15_sankey_manual_vs_sctype.png",
       p_sankey_sctype, width = 12, height = 10, dpi = 300)

p_sankey_sccatch <- create_sankey_plot(sankey_manual_sccatch, "Manual", "scCATCH")
ggsave("plots/method_comparison/16_sankey_manual_vs_sccatch.png",
       p_sankey_sccatch, width = 12, height = 10, dpi = 300)



#-----------------------------------------------
# STEP 17: Create comprehensive UMAP comparison
#-----------------------------------------------

# Create individual plots with consistent styling
p_manual_comp <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "manual_annotation",
  pt.size = 0.1,
  raster=FALSE
) +
  labs(title = "Manual") +
  theme(legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold"))

p_singler_comp <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "singler_annotation",
  pt.size = 0.1,
  raster=FALSE
) +
  labs(title = "SingleR (Mouse)") +
  theme(legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold"))

p_sctype_comp <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "sctype_annotation",
  pt.size = 0.1,
  raster=FALSE
) +
  labs(title = "scType (Markers)") +
  theme(legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold"))

p_sccatch_comp <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "sccatch_annotation",
  pt.size = 0.1,
  raster=FALSE
) +
  labs(title = "scCATCH (Tissue DB)") +
  theme(legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold"))

# Combine all four methods
# p_all_methods <- (p_manual_comp | p_singler_comp) /
#   (p_sctype_comp | p_sccatch_comp)
p_all_methods <- (p_manual_comp | p_sccatch_comp)

ggsave("plots/method_comparison/17_all_methods_comparison_sccatch.png",
       p_all_methods, width = 18, height = 14, dpi = 300)



#-----------------------------------------------
# STEP 18: Export annotations for manual comparison
#-----------------------------------------------

# Create comprehensive annotation table
annotation_comparison <- data.frame(
  cell_barcode = colnames(integrated_seurat),
  cluster = integrated_seurat$seurat_clusters,
  manual = integrated_seurat$manual_annotation,
  # singler = integrated_seurat$singler_annotation,
  # sctype = integrated_seurat$sctype_annotation,
  sccatch = integrated_seurat$sccatch_annotation,
  sample_id = integrated_seurat$sample_id,
  condition = integrated_seurat$condition,
  stringsAsFactors = FALSE
)

# Export to CSV for manual review
write.csv(annotation_comparison,
          "annotations/all_methods_comparison_sccatch.csv",
          row.names = FALSE)

# Create cluster-level summary (easier to review)
cluster_summary <- annotation_comparison %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    manual_type = names(sort(table(manual), decreasing = TRUE))[1],
    # singler_type = names(sort(table(singler), decreasing = TRUE))[1],
    # sctype_type = names(sort(table(sctype), decreasing = TRUE))[1],
    sccatch_type = names(sort(table(sccatch), decreasing = TRUE))[1],
    .groups = "drop"
  )

write.csv(cluster_summary,
          "annotations/cluster_level_comparison_sccatch.csv",
          row.names = FALSE)



#-----------------------------------------------
# STEP 19: Import consensus annotations
#-----------------------------------------------

# Create mapping from cluster to consensus annotation
cluster_to_consensus <- setNames(
  c("endothelium", "Pneumocyte type II", "epithelial", "fibroblast", "lung club cell", "macrophage", "fibroblast", "B cell", 
    "epithelial", "T cell", "macrophage", "endothelium", "fibroblast", "endothelium", "immune", 
    "epithelial", "endothelium", "fibroblast", "immune", "immune", "immune", "epithelial", "endothelium"),  
  c(0:22)
)

# Apply consensus to all cells
integrated_seurat$final_annotation <- unname(cluster_to_consensus[
  as.character(integrated_seurat$seurat_clusters)
])

# Visualize final consensus
p_consensus <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "final_annotation",
  label = TRUE,
  label.size = 4,
  pt.size = 0.1,
  repel = TRUE,
  raster=FALSE
) +
  ggtitle("Final Consensus Annotation") +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("plots/method_comparison/18_final_consensus_annotation.png",
       p_consensus, width = 11, height = 8, dpi = 300)

# Split by condition to check biological validity
p_consensus_split <- DimPlot(
  integrated_seurat,
  reduction = umap_reduction,
  group.by = "final_annotation",
  split.by = "condition",
  pt.size = 0.1,
  ncol = 2,
  raster=FALSE
) +
  ggtitle("Final Consensus by Condition")

ggsave("plots/method_comparison/19_consensus_by_condition.png",
       p_consensus_split, width = 14, height = 6, dpi = 300)



#-----------------------------------------------
# STEP 20: Save final annotated Seurat object
#-----------------------------------------------

# Save the complete annotated Seurat object
saveRDS(integrated_seurat, "annotations/integrated_annotated_seurat.rds")

# Create metadata export
metadata_export <- integrated_seurat@meta.data %>%
  dplyr::select(
    sample_id, condition, patient_id,
    seurat_clusters,
    manual_annotation,
    # singler_annotation,
    # sctype_annotation,
    # sccatch_annotation,
    final_annotation
  )

write.csv(metadata_export,
          "annotations/cell_metadata_final.csv",
          row.names = TRUE)




### Code Nico pour trouver les cellules tumorales, une fois les cellules immunitaires trouvées 
#https://github.com/broadinstitute/infercnv/wiki/infercnv-10x
library(infercnv)
rawcount <- as.matrix(LayerData(obj, assay = "RNA", layer = "counts"))

#sampleannotations
sample.annotation <- obj@meta.data[, c("celltype", "orig.ident")]
sample.annotation$orig.ident <- rownames(sample.annotation)
head(sample.annotation)
colnames(sample.annotation) <- NULL


geneorder <- read.table("/Users/nicolasrama/Desktop/mouse_gencode.GRCm39.vM32.basic.annotation.by_gene_name.infercnv_positions.txt", header = T)
head(geneorder)
rownames(geneorder) <- geneorder$GeneId
geneorder <- geneorder[, -1]
head(geneorder)

table(colnames(rawcount) == rownames(sample.annotation))

infercnv_obj = CreateInfercnvObject(raw_counts_matrix = rawcount, annotations_file = sample.annotation, gene_order_file = geneorder, ref_group_names = c("Neutrophils", "Mono/macrophage", "lymphocyteT", "LymphocyteB"))



infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=tempfile(),
                             cluster_by_groups=TRUE,
                             denoise=TRUE,
                             HMM=TRUE, num_threads = 1)


seurat_obj = infercnv::add_to_seurat(infercnv_output_path='/var/folders/j2/0rkf_l2x10z866tvy9bss1_h0000gn/T/RtmpB4z3ET/file58e1e24729f/',
                                     seurat_obj=obj, # optional
                                     top_n=10
)


FeaturePlot(seurat_obj, features = "proportion_scaled_cnv_chr19", order = T, label=F, split.by = "orig.ident") #
FeaturePlot(seurat_obj, features = "mitoRatio", order = T, label=T, min.cutoff = "q5") #

Idents(obj) <- obj@meta.data$SCT_snn_res.0.2
DimPlot(obj, split.by = "orig.ident")