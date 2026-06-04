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
setwd("/home/eugenie-modolo/Documents/Lapnet/donnees_KRAS/GSE303105_RAW/differential_analysis")

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





