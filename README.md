<!--
SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>

SPDX-License-Identifier: AGPL-3.0-or-later
-->

# Mehlen_LapNet
R scripts from project LapNet, carried out at CRCL (Centre de Recherche en Cancérologie de Lyon) facilities, using data from biotech Netris pharma.


## Sub-project 1: bulk RNA-Seq data from clinical study (phase Ib)

The 4 files cited hereafter were used to study the 38 samples obtained from 30 patients (of which 6 samples are from 3 negative-control patients) during the LapNet clinical trial. 

The goal of this trial was to assess if the addition of NP137 (Netrin-1 inhibitor) to the standard-of-care (modified FOLFIRINOX) would increase survival in first-line patients with locally advanced pancreatic cancer (LAPC). For this, samples where collected after receiving the treatment on 27 patients and 3 controls.

This clinical trial led to the publication of a paper in Nature: [Netrin1 blockade alleviates resistance to chemotherapy in pancreatic cancer](https://www.nature.com/articles/s41586-026-10436-4)

- rapport_LapNet_1_QC.Rmd: QC (Quality Control) checks on the 38 samples (PCA, hierarchical clustering heatmaps), in order to see if there are samples behaving incoherently, and decide if they should be removed from further analyses or not. The secondary goal is to find which variables better explain the two variables of interest of this study: survival and difference of expressions before and after receiving the treatment.
As there is no consensus on the number of most-variable genes that should be used for the PCA, I decided to implement an elbow method on the variance (PC1 + PC2) explained by the number of top variable genes. 
As the ultimate goal will be to choose which variable(s) to include in the DESeq2 design (cf. files rapport_LapNet_3_choixDesign.Rmd & rapport_LapNet_4_DEGs.Rmd), I also performed PLS (Partial Least Square) and the top variables contributing to the PCA variances, in order to have a first idea of which clinical variables most influence the gene expression in the samples.

- rapport_LapNet_2_survival.Rmd: study of the link between expression of Netrin-1, and its 15 main known receptors as well as global EMT (Epithelial-Mesenchymal Transition) levels, and the patients' survival (OS & PFS). This study notably uses Kaplan-Meier and ssGSEA statistical tests. These results are also compared to the samples from a comparable study: [Predictive genomic and transcriptomic analysis on endoscopic ultrasound-guided fine needle aspiration materials from primary pancreatic adenocarcinoma: a prospective multicentre study](https://pubmed.ncbi.nlm.nih.gov/39383608/). 

- rapport_LapNet_3_choixDesign.Rmd: In addition to the PLS studies performed in the file rapport_LapNet_1_QC.Rmd, I performed here LRT (Likelihood Ratio Tests) with DESeq2 on the two different biological questions (difference in survival and difference of expressions before and after receiving the treatment), in order to automatically choose (i) the variable of interest when several clinical variables could be used; (ii) the best variables to include in the design to explain the noise in the data, so that more DEGs (Differentially-Expressed Genes) could be found.
In addition, regarding the samples that had a low QC, I checked the number of genes with a null expression, in order to definitely decide which samples to exclude from the study.

- rapport_LapNet_4_DEGs.Rmd: study of the DEGs found for both biological questions. I first choose the best DESeq2 design based on the heatmaps and Venn diagram, then perform cross-validation (re-run DESeq2 with the same design on all the samples minus one, for each sample) in order to keep only the DEGs validated by at least 75% of the samples. Finally I run the classical bio-informatics study on the DEGs found (ORA - Over-Enrichment Analysis, GSEA - Gene Set Enrichment Analysis).



## Sub-project 2: bulk & single-cell RNA-Seq data from external studies

The files look into the data from two studies on the effect of Daraxonrasib (or RMC-6236, KRAS inhibitor), on the gene expression of human or mouse cancerous samples treated by Daraxonrasib.

- rapport_LapNet_postDaraxonrasib.Rmd: study the Netrin-1 expression, as well as it main receptors, on bulk RNA-Seq sequencing of 2 human and 2 mouse cell lines. Data from [Differential mRNA expression analysis of multiple human and mouse pancreatic cancer cells and tumors after DMSO, tazemetostat, RMC-6236, or the combination treatment](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315130).

- single_cell_*.R: study the Netrin-1 expression, as well as it main receptors, on RNA-Seq sequencing of 2 mouse single-cell samples (1 that received the treatment and 1 that did not). The different files are the various steps to process single-cell data. Data from [Prolonged KRAS-MAPK Inhibition Induces Interferon-mediated Epithelial-to-Mesenchymal Transition and Reveals Therapeutic Opportunities](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303105).

