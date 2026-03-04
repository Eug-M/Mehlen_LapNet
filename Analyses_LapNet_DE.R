# SPDX-FileCopyrightText: 2026 Eugenie Modolo <eugenie.modolo@lyon.unicancer.fr>
#
# SPDX-License-Identifier: AGPL-3.0-or-later

### Remarques
# pour prendre en compte la time-series de taille de tumeur : attention, l'écart entre date de screen et C4 n'est pas constant
# (contrairement aux écarts entre C4, C8 et C12 qui le sont à peu près) -> prendre en compte les dates dans le modèle ?


### Libraries
#BiocManager::install("DESeq2")
library("DESeq2")
library("SummarizedExperiment")
library("pheatmap")
library("RColorBrewer")
library("dplyr")
library("ggplot2")
library("apeglm")
library("ggrepel")
library("DEGreport")
#library("pcaExplorer") #incompatible avec bcp de packages ...



### Creating colData

coldata <- read.csv('/home/eugenie-modolo/Documents/Lapnet/colData_LapNet.csv', sep=";")
coldata$Best_response_2 <- case_when(
  coldata$Best_response == "SD" ~ "SD_PD",
  coldata$Best_response == "PD" ~ "SD_PD", 
  coldata$Best_response == "PR" ~ "PR",
  TRUE ~ NA_character_  # for any unexpected values
)

# scale and center the numeric columns
numeric_columns <- c('Age', 'Height', 'Weight', 'BMI', 'Size_lesion_C0', 'Size_lesion_C4', 'Size_lesion_C8', 'Size_lesion_C12', 
                     'best_lesion_perc', 'OS', 'PFS', 'Surgery', 'QC_duplication_perc', 'QC_aligned_M', 'QC_uniq_aligned')
coldata[numeric_columns] <- lapply(coldata[numeric_columns], scale)

# convert to factors the non-numeric columns
'%!in%' <- function(x,y)!('%in%'(x,y))
factor_columns <- colnames(coldata) %!in% numeric_columns
coldata[factor_columns] <- lapply(coldata[factor_columns], factor)



### Generating countData from Salmon transcript counts (less biased than Star_salmon counts)
counts_transcripts <- read.table("/home/eugenie-modolo/Documents/Lapnet/results_rnaseq_NFcore/salmon.merged.transcript_counts.tsv", 
                           header = TRUE, 
                           row.names = 1, 
                           sep = "\t",
                           check.names = FALSE)
count_columns <- colnames(counts_transcripts)[colnames(counts_transcripts) != "gene_id"]
counts_transcripts[count_columns] <- round(counts_transcripts[count_columns])
count_matrix <- as.matrix(counts_transcripts[, count_columns])

basic_coldata <- data.frame(coldata, row.names = count_columns)



### Generating countData from Salmon gene counts 
counts_genes <- read.table("/home/eugenie-modolo/Documents/Lapnet/results_rnaseq_NFcore/salmon.merged.gene_counts.tsv", 
                     header = TRUE, 
                     row.names = 1, 
                     sep = "\t",
                     check.names = FALSE)
count_columns <- colnames(counts_genes)[colnames(counts_genes) != "gene_name"]
# besoin de prendre les integers pour DESeq2
counts_genes[count_columns] <- round(counts_genes[count_columns])
count_matrix <- as.matrix(counts_genes[, count_columns])
rownames(count_matrix) <- counts_genes$gene_name  # Gene names become row names
# # on modifie tous les noms de gènes dupliqués, en les remplaçant par leur ENSG code, pour éviter le warning DESeq2
# dup_genes <- counts_genes[duplicated(counts_genes$gene_name),'gene_name']
# counts_genes$feature <- ifelse((is.na(counts_genes$gene_name) | counts_genes$gene_name %in% dup_genes), rownames(counts_genes), counts_genes$gene_name)
# count_columns <- colnames(counts_genes)[colnames(counts_genes) != "gene_name" & colnames(counts_genes) != "feature"]
# count_matrix_lapnet <- as.matrix(counts_genes[, count_columns])
# rownames(count_matrix_lapnet) <- counts_genes$feature

basic_coldata <- data.frame(coldata, row.names = count_columns)

# vérif pour DESeq2 : the rows of the metadata table and the columns of the expression matrix must have the same name and same order
sum(rownames(basic_coldata) == colnames(count_matrix)) == length(rownames(basic_coldata))
# all(colnames(basic_coldata) == rownames(count_matrix))

# /!\ la plateforme conseille de : retirer les gènes presque jamais exprimés, et faire shrink /!\



### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58



### Generating the different DESeqDataSets (according to the design formulae)
# /!\ Note: In order to benefit from the default settings of the package, you should put the variable of interest at the end of the formula and make sure the control level is the first level.
# This genotype-condition interaction example is examined in further detail in Example 3 in the help page for results, which can be found by typing ?results. 
# In particular, we show how to test for differences in the condition effect across genotype, and we show how to obtain the condition effect for non-reference genotypes.
dds_pre <- DESeqDataSetFromMatrix(
  countData = count_matrix[, basic_coldata$Time == "pre"],
  colData = basic_coldata[basic_coldata$Time == "pre", ],
  design = ~ Best_response_2
)
# /!\ 4307 duplicate rownames were renamed by adding numbers /!\

dds_pre <- DESeq(dds_pre)
#dds_prepost = DESeqDataSet(se, design = ~ xxx) # bien verifier que le pre est le premier dans le tableau (avant post), car cest notre reference
norm.count <- as.data.frame(counts(dds_pre, normalized=TRUE))

res_pre <- results(dds_pre, alpha=padj.cutoff) #, contrast=c("type", "single", "paired")) # a modifier
summary(res_pre)

resOrdered <- res_pre[order(res_pre$pvalue),] # 25 gènes ont padj<5% (mais je devrais regarder FDR de 10% ??)



### Data overview
mcols(res_pre)$description # checks the variables and tests which were used

## checking the DESeq2 dispersion (doit faire ~elbow)
plotDispEsts(dds_pre)

## looking at the number of outliers
# papier DESeq2 : "Cook’s distance is defined within each gene for each sample as the scaled distance that the coefficient vector
# of a linear model or GLM would move if the sample were removed and the model refit."
# Il faut que ce soit <0. Rq : ça dépend du design, car il utilise l'info de replicats par condition (il en faut au moins 3)
boxplot(log10(assays(dds_pre)[["cooks"]]), range=0, las=2)

## Plot counts for genes of interest (Netrin, ...) for the condition
d <- plotCounts(dds_pre, gene=which.min(res_pre$padj), intgroup="Best_response_2", returnData=TRUE) # a changer : on ne veut pas le gene avec la plus petite p-value
ggplot(d, aes(x=Best_response_2, y=count, color=Best_response_2)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("xxx") +
  theme(plot.title = element_text(hjust = 0.5))

## Plot counts of several genes (e.g. all DE genes)
# papier DESeq2 : "the sample covariate information (e.g. treatment or control) is not used, so that all samples are 
# treated equally" -> c'est indépendant du design
#res_pre_shrunken <- lfcShrink(dds_pre, contrast=c("Sex", "F", "M"), res=res) # 'M' is the control in the log2FC obtained
res_pre_shrunken <- lfcShrink(dds_pre, coef='Age')
normalized_counts <- counts(dds_pre, normalized=TRUE)
res_table_tb <- res_pre_shrunken %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
top20_sig_genes <- res_table_tb %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  pull(gene) %>% 		#Extract character vector of ordered genes
  head(n=20) 		#Extract the first 20 genes
top20_sig_norm <- normalized_counts %>%
  filter(gene %in% top20_sig_genes)
gathered_top20_sig <- top20_sig_norm %>%
  gather(colnames(top20_sig_norm)[2:9], key = "xxx", value = "normalized_counts") # 2:9 à changer ? / xxx = colonne du design
# vérif
View(gathered_top20_sig)
gathered_top20_sig <- inner_join(basic_coldata_pre, gathered_top20_sig)
ggplot(gathered_top20_sig) +
  geom_point(aes(x = gene, y = normalized_counts, color = xxx)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

## MA-plots: use Glimma / Volcano plots: use ? 
plotMA(res, ylim=c(-2,2))
# plotMA(res_pre_shrunken, ylim=c(-2,2))

## PCA
# /!\ Rq Laurie : il vaut mieux ne regarder que les gènes les plus variants (car y'en a qui peuvent tirer la variance pour rien) ;
# les outliers qu'on remarque : c'est normal avec des échantillons patients, ne rien exclure basé sur les PCA a priori, mais 
# regarder ensuite l'expression par échantillon des gènes DE et si un échantillon qui était outlier sur la PCA porte bcp de la 
# variabilité d'expression qui fait que le gène a été classé DE, alors refaire l'étude sans cet échantillon pour voir 
rld <- rlog(dds_pre, blind=TRUE)
# rld <- vst(dds_pre, blind=TRUE) # plus rapide, OK pour plotPCA()
# PC1 & PC2 :
plotPCA(rld, intgroup="Age", ntop=1000) # default ntop: 500
# 06-008 & 06-009 sont extrèmes ...
# Autres PC, méthode 1 :
pcaData <- plotPCA(rld, intgroup=c("Best_response_2", "Sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Best_response_2, shape=Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
# Autres PC, méthode 2 :
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
df <- cbind(basic_coldata_pre, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = SampleID))

## Hierarchical Clustering Heatmap
# Samples below 0.80 may indicate an outlier in your data and/or sample contamination
rld_mat <- assay(rld) 
rld_cor <- cor(rld_mat, method="spearman")
#head(rld_cor)
pheatmap(rld_cor)
# heat.colors <- brewer.pal(6, "Blues")
# pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10, 
#          fontsize_row = 10, height=20)


## Volcano plots
# Méthode 1
DEGreport::degPlot(dds = dds, res = res, n = 20, xs = "type", group = "condition") # dds object is output from DESeq2
DEGreport::degVolcano(
  data.frame(res[,c("log2FoldChange","padj")]), # table - 2 columns
  plot_text = data.frame(res[1:10,c("log2FoldChange","padj","id")])) # table to add names
DEGreport::degPlotWide(dds = dds, genes = row.names(res)[1:5], group = "condition")

# Methode 2
# Obtain logical vector where TRUE values denote padj values < padj.cutoff and fold change > 1.5 in either direction
res_table_tb <- res_tableOE_tb %>% 
  mutate(threshold = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
# Create a column to indicate which genes to label on the plot
res_tableOE <- res_tableOE %>% arrange(padj) %>% mutate(genelabels = "")
res_tableOE$genelabels[1:10] <- res_tableOE$gene[1:10] # top 10 genes (lowest padj)
# plot
ggplot(res_table_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("xxx") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  



### Data transformation and exploration
# /!\ Rq Laurie : ne PAS utiliser les données de survie pour déterminer les gènes DE, car ces données vont être utilisées ensuite
# pour voir si les gènes signature DE trouvés (selon Best_response, pre/post ou autre) ont bien un lien avec la survie -> ça 
# validera notre hypothèse
rld <- rlog(dds_pre, blind=FALSE)

# Heatmap of the count matrix
select <- order(rowMeans(counts(dds_pre,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_pre)[,c("Best_response_2", "Sex")])
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$Best_response_2, rld$Sex, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Heatmap of the expression of all the significant genes
norm_sig <- normalized_counts[,c(1,4:9)] %>% 
  filter(gene %in% sig$gene) %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 
annotation <- basic_coldata_pre %>% 
  select(xxx, sampletype) %>% 
  data.frame(row.names = "xxx")
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(norm_sig, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", # plotting Z-scores: computed (after the clustering, to improve the visualization) on a gene-by-gene basis by subtracting the mean and then dividing by the standard deviation
         fontsize_row = 10, 
         height = 20)



### Pairwise comparisons with Likelihood ratio test
# Note on LRT tests: The p-values are determined solely by the difference in deviance between the ‘full’ and ‘reduced’ model formula (not log2 fold changes).
# Generally, this test will result in a larger number of genes than the individual pair-wise comparisons. -> we can lower padj.cutoff (e.g. 0.001)
# one should not expect it to be exactly equal to the union of sets of genes using Wald tests (although we do expect a majority overlap).
# [...] even though there are fold changes present in the results table, they are not directly associated with the 
# actual hypothesis test. Thus, when filtering significant genes from the LRT we use only the FDR as our threshold.

## 23 échantillons pré : étude des likelihood ratio test selon le design
dds_pre <- DESeqDataSetFromMatrix(
  countData = count_matrix[, basic_coldata$Time == "pre"],
  colData = basic_coldata[basic_coldata$Time == "pre", ],
  design = ~ Best_response_2
)
dds_pre <- DESeq(dds_pre, test="LRT", reduced=~1)
res_lrt <- results(dds_pre, alpha=padj.cutoff)
summary(res_lrt)

## 22 échantillons pré : étude des likelihood ratio test selon le design
dds_pre <- DESeqDataSetFromMatrix(
  countData = count_matrix[, basic_coldata$Time == "pre" & basic_coldata$SampleID != '06-006'],
  colData = basic_coldata[basic_coldata$Time == "pre" & basic_coldata$SampleID != '06-006', ],
  design = ~ Best_response_2
)
dds_pre <- DESeq(dds_pre, test="LRT", reduced=~1)
res_lrt <- results(dds_pre, alpha=padj.cutoff)
summary(res_lrt)

## pré vs post
# 23 échantillons
dds_prepost <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = basic_coldata,
  design = ~ Time + SampleID
)
dds_prepost <- DESeq(dds_prepost, test="LRT", reduced=~Time)
res_lrt <- results(dds_prepost, alpha=padj.cutoff)
summary(res_lrt)

# 22 échantillons
dds_prepost <- DESeqDataSetFromMatrix(
  countData = count_matrix[, basic_coldata$SampleID != '06-006'],
  colData = basic_coldata[basic_coldata$SampleID != '06-006', ],
  design = ~ Time + SampleID
)
dds_prepost <- DESeq(dds_prepost, test="LRT", reduced=~Time)
res_lrt <- results(dds_prepost, alpha=padj.cutoff)
summary(res_lrt)

# 4 paires
pairs_ID <- c('01-007', '02-001', '02-003', '02-014')
dds_prepost <- DESeqDataSetFromMatrix(
  countData = count_matrix[, basic_coldata$SampleID %in% pairs_ID],
  colData = basic_coldata[basic_coldata$SampleID %in% pairs_ID, ],
  design = ~ Time + SampleID
)
dds_prepost <- DESeq(dds_prepost, test="LRT", reduced=~Time)
res_lrt <- results(dds_prepost, alpha=padj.cutoff)
summary(res_lrt)



### Other interesting ideas to explore
## Identifying gene clusters exhibiting particular patterns across samples
# exemple : gènes qui sont bas dans condition A, moyen dans condition B et hauts dans condition C
#library(DEGreport)
sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)
clustering_sig_genes <- sig_res_LRT %>%
  arrange(padj)
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "sampletype", col=NULL)
# chercher à faire la figure comme dans le cours 8, en regardant class(clusters) et head(clusters$df)
# regarder les gènes présents dans le(s) groupe(s) qui nous intéresse(nt) : 
cluster_groups <- clusters$df
group1 <- clusters$df %>%
  filter(cluster == 1)
# puis perform functional analysis to explore associated functions.

## Time course analyses with LRT
# This analysis will not return genes where the treatment effect does not change over time, even though the genes 
# may be differentially expressed between groups at a particular time point. The significant DE genes will represent
# those genes that have differences in the effect of treatment over time
# termes à modifier /!\
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ genotype + treatment + time + treatment:time)
dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ genotype + treatment + time)
# determine the significant genes with threshold of padj < 0.05 
clusters <- degPatterns(cluster_rlog, metadata = meta, time="time", col="treatment")
# extract the groups of genes associated with the patterns of interest similar to the actions performed previously, 
# then move on to functional analysis for each of the gene groups of interest



### Functional analysis (to be done in another R script !)
# Use DEvis ?

## Gene set selection: MSigDB

## ORA (Over-Representation Analysis): Gene Ontology with clusterProfiler
# 3 different GO ontology to look at: BP (Biological process), CC (Cellular component), MF (Molecular function) 
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)
# convert gene symbols into Ensembl IDs
# à modifier ! https://hbctraining.github.io/DGE_workshop/img/orgdb_annotation_databases.png
idx <- grch37$symbol %in% rownames(res_tableOE)
ids <- grch37[idx, ]
non_duplicates <- which(duplicated(ids$symbol) == FALSE)
ids <- ids[non_duplicates, ] 
res_ids <- inner_join(res_table_tb, ids, by=c("gene"="symbol"))  
all_genes <- as.character(res_ids$ensgene) # background dataset for hypergeometric testing using all genes tested for significance in the results
sig <- filter(res_ids, padj < 0.05) # Extract significant results
sig_genes <- as.character(sigOE$ensgene)

# Run GO enrichment analysis 
ego <- enrichGO(gene = sig_genes, 
                universe = all_genes,
                keyType = "ENSEMBL", # si ça marche pas, taper : keytype
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Visualizing results
dotplot(ego, showCategory=50) # displays the top 50 genes by gene ratio (= # genes related to GO term / total number of sig genes)
# To save the figure, click on the Export button in the RStudio Plots tab and Save as PDF. In the pop-up window, change:
# Orientation: to Landscape / PDF size to 8 x 14 to give a figure of appropriate size for the text labels
emapplot(ego, showCategory = 50) # enrichment GO plot: relationship between the top 50 most significantly enriched GO terms (padj.), by grouping similar terms together
# To save the figure, idem above but PDF size to 24 x 32
cnetplot() # category netplot, cf. cours 9. relationships between the genes associated with the top five most significant GO terms and the fold changes of the significant genes associated with these terms 

# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "results/xxx.csv")

#gProfileR with REVIGO
# cf. cours : https://hbctraining.github.io/DGE_workshop/lessons/functional_analysis_other_methods.html


## Functional class scoring tools (FCS): GSEA 
# plateforme: tester GSEA Preranked (en plus de fGSEA), car ça marche mieux lorsque peu de samples
# This type of analysis can be particularly helpful if the differential expression analysis only outputs a small list of significant DE genes.
# lfcShrink() to be used when we want to subset the significant genes based on fold change for further evaluation & for functional analysis tools (such as GSEA) which require fold change values as input
# plateforme: When using fGSEA, keep in mind to change the parameters MaxSise and minSize, to 1000 and 10 respectively, to avoid losing pathways
# plateforme: Be careful if you have very few differentially expressed genes in your analysis. GSEA can give you strong results that will be incoherent with the differential analysis
res_pre_shrunken # à transformer en foldchanges, cf. cours 9
biocLite("GSEABase")
library(GSEABase)
c2 <- read.gmt("/data/xxx.txt") # fichier MSigDB
msig <- GSEA(foldchanges, TERM2GENE=c2, verbose=FALSE)
msig_df <- data.frame(msig)



### Pathway Topology (PT) -> cf. cours 9



### Export results
resSig <- subset(resOrdered, padj < 0.1) # Exporting only the results which pass an adjusted p value threshold
write.csv(as.data.frame(resSig), 
          file="condition_treated_results.csv")



### Output the versions of all tools used in the DE analysis
sessionInfo()



#### Ancient code

### Creating coldata
# samplesID <- dir('/home/eugenie-modolo/Documents/Lapnet', pattern="0.-0.._")
# #samplesID_pre <- dir('/home/eugenie-modolo/Documents/Lapnet', pattern="0.-0.._pre")
# #samplesID_post <- dir('/home/eugenie-modolo/Documents/Lapnet', pattern="0.-0.._post")
# 
# # col_time = unlist( lapply(samplesID,function(elem) strsplit(elem,"_")[[1]][2]) )
# # #col_time = unlist( lapply(samplesID_pre,function(elem) strsplit(elem,"_")[[1]][2]) )
# # #col_time = unlist( lapply(samplesID_post,function(elem) strsplit(elem,"_")[[1]][2]) )
#coldata$OS_modif <- ifelse(is.na(coldata$OS), 900, coldata$OS)
#coldata_pre <- read.csv('/home/eugenie/Documents/Lapnet/colData_LapNet_pre.csv', sep=";")
#coldata_post <- read.csv('/home/eugenie/Documents/Lapnet/colData_LapNet_post.csv', sep=";")

### Merging with the tximeta output from NF-CORE
# # Convert the counts assay to a numeric matrix
# se <- readRDS("/home/eugenie-modolo/Documents/Lapnetsalmon.merged.gene.SummarizedExperiment.rds")
# counts_matrix <- as.matrix(assay(se, "counts"))
# storage.mode(counts_matrix) <- "numeric"
# counts_matrix <- round(assay(counts_matrix))
# assay(se, "counts") <- counts_matrix
# colData(se) <- DataFrame(basic_coldata)
# #se_pre <- se[, !colnames(se) %in% samplesID_post]
# se_pre <- se[, colData(se)$time == "pre"]
# se_post <- se[, colData(se)$time == "post"]
# dds_pre <- DESeqDataSet(se_pre, design = ~ Best_response_2) # marche pas !!

### Generating countData from Salmon transcript counts (less biased than Star_salmon counts)
# # a modifier selon la vignette tximport / tximeta !!!
# dir <- system.file("extdata", package="tximportData")
# samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
# samples$condition <- factor(rep(c("A","B"),each=3))
# rownames(samples) <- samples$run
# samples[,c("pop","center","run","condition")]
# files <- file.path(dir,"salmon", samples$run, "quant.sf") # "quant.sf.gz" accepted
# names(files) <- samples$run
# tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
# txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
# # OU
# coldata <- samples
# coldata$files <- files
# coldata$names <- coldata$run
# se <- tximeta(coldata)
# ddsTxi <- DESeqDataSet(se, design = ~ condition)
# # pour utiliser les noms des samples pour creer mes conditions, modifier ce code :
# levels(ddsMF$type) <- sub("-.*", "", levels(ddsMF$type)) # modifie "paired-end"  "single-read" en "paired" "single"
#dds_prepost = DESeqDataSetFromTximport(txi, colData = basic_coldata, design = ~ xxx)
