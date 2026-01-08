
### Libraries
#BiocManager::install("DESeq2")
library("RColorBrewer")
library("dplyr")
library("ggplot2")
library(gridExtra)
library("tidyverse")
#library("hrbrthemes")
library("viridis")
library('survival')
library(ranger)
library(ggfortify)
library(survminer)
library(KMEANS.KNN)
library("DESeq2")
library(factoextra)
library(GSVA)
library(mclust)



### Creating colData
coldata <- read.csv('/home/eugenie-modolo/Documents/Lapnet/colData_LapNet.csv', sep=";")
coldata$Best_response_2 <- case_when(
  coldata$Best_response == "SD" ~ "SD_PD",
  coldata$Best_response == "PD" ~ "SD_PD", 
  coldata$Best_response == "PR" ~ "PR",
  TRUE ~ NA_character_  # for any unexpected values
)
#coldata$OS_modif <- ifelse(is.na(coldata$OS), 900, coldata$OS)



### Generating countData from Salmon transcript counts (less biased than Star_salmon counts)
counts_transcripts <- read.table("/home/eugenie-modolo/Documents/Lapnet/results_rnaseq_NFcore/salmon.merged.transcript_tpm.tsv", 
                           header = TRUE, 
                           row.names = 1, 
                           sep = "\t",
                           check.names = FALSE)
count_columns <- colnames(counts_transcripts)[colnames(counts_transcripts) != "gene_id"]
counts_transcripts[count_columns] <- round(counts_transcripts[count_columns])
count_matrix_trans <- as.matrix(counts_transcripts[, count_columns])



### Generating countData from Salmon gene counts 
counts_genes <- read.table("/home/eugenie-modolo/Documents/Lapnet/results_rnaseq_NFcore/salmon.merged.gene_tpm.tsv", 
                     header = TRUE, 
                     row.names = 1, 
                     sep = "\t",
                     check.names = FALSE)
count_columns <- colnames(counts_genes)[colnames(counts_genes) != "gene_name"]
#counts_genes[count_columns] <- round(counts_genes[count_columns])
count_matrix <- as.matrix(counts_genes[, count_columns])
rownames(count_matrix) <- counts_genes$gene_name  # Gene names become row names

basic_coldata <- data.frame(coldata, row.names = count_columns)

# verif que les genes d interet sont presents dans le count_matrix
c("NTN1", "NEO1", "UNC5B") %in% counts_genes$gene_name

# verif que les genes d interet sont presents dans le count_matrix_trans
ensg_netrin <- rownames(which(counts_genes == "NTN1", arr.ind = T))
enst_netrin <- rownames(which(counts_transcripts == ensg_netrin, arr.ind = T))
ensg_neogenin <- rownames(which(counts_genes == "NEO1", arr.ind = T))
enst_neogenin <- rownames(which(counts_transcripts == ensg_neogenin, arr.ind = T))
ensg_uncb <- rownames(which(counts_genes == "UNC5B", arr.ind = T))
enst_uncb <- rownames(which(counts_transcripts == ensg_uncb, arr.ind = T))

count_matrix_trans[which(row.names(count_matrix_trans) %in% enst_netrin),]
count_matrix_trans[which(row.names(count_matrix_trans) %in% enst_neogenin),]
count_matrix_trans[which(row.names(count_matrix_trans) %in% enst_uncb),]



### Data overview

# Plot counts for genes of interest (Netrin, Neogenin, ...) for the condition
expression_genes <- count_matrix[which(row.names(count_matrix) %in% c('NTN1', 'NEO1', 'UNC5B')),]
t_expression_genes <- data.frame(t(expression_genes))
basic_coldata$NTN1 <- t_expression_genes$NTN1
basic_coldata$NTN1_log <- log2(t_expression_genes$NTN1 + 1)
basic_coldata$NEO1 <- t_expression_genes$NEO1
basic_coldata$NEO1_log <- log2(t_expression_genes$NEO1 + 1)


# verif que NEO1_log suit une gaussienne
ggplot(basic_coldata[which(basic_coldata$Time == 'pre'),], 
       aes(x=NEO1_log)) + 
  geom_histogram()

# p1 <- ggplot(basic_coldata, aes(x=Best_response_2, y=NTN1, fill=Best_response_2)) + 
#   geom_violin() +
#   ggtitle("Violin chart NTN1 all samples")
# p2 <- ggplot(basic_coldata, aes(x=Best_response_2, y=NEO1, fill=Best_response_2)) + 
#   geom_violin() +
#   ggtitle("Violin chart NEO1 all samples")

p1 <- ggplot(basic_coldata[which(basic_coldata$Time == 'pre'),], 
       aes(x=Best_response_2, y=NTN1, fill=Best_response_2)) + 
  geom_violin() +
  ggtitle("Violin chart NTN1 pre samples")
p2 <- ggplot(basic_coldata[which(basic_coldata$Time == 'pre'),], 
             aes(x=Best_response_2, y=NTN1_log, fill=Best_response_2)) + 
  geom_violin() +
  ggtitle("Violin chart log2(NTN1+1) pre samples")
p3 <- ggplot(basic_coldata[which(basic_coldata$Time == 'pre'),], 
       aes(x=Best_response_2, y=NEO1, fill=Best_response_2)) + 
  geom_violin() +
  ggtitle("Violin chart NEO1 pre samples")
p4 <- ggplot(basic_coldata[which(basic_coldata$Time == 'pre'),], 
             aes(x=Best_response_2, y=NEO1_log, fill=Best_response_2)) + 
  geom_violin() +
  ggtitle("Violin chart log2(NEO1+1) pre samples")
grid.arrange(p1, p2, p3, p4, nrow = 2)

p5 <- ggplot(basic_coldata[which(basic_coldata$Time == 'pre'),], 
             aes(x=Best_response_2, y=NEO1_log, fill=Best_response_2)) + 
  geom_boxplot(notch = TRUE) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggtitle("Boxplot log2(NEO1+1) pre samples")
p6 <- ggplot(basic_coldata, 
             aes(x=Best_response_2, y=NEO1_log, fill=Best_response_2)) + 
  geom_boxplot(notch = TRUE) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggtitle("Boxplot log2(NEO1+1) all samples")
grid.arrange(p5, p6, nrow = 1)

# Choosing visually the high/low limit (should use Hector's technique)
basic_coldata$NEO1_level <- case_when(
  basic_coldata$NEO1_log >= 2.75 ~ "High",
  basic_coldata$NEO1_log < 2.75 ~ "Low")
sum(basic_coldata$NEO1_level == "High")  # 14/29
sum(basic_coldata[which(basic_coldata$Time == 'pre'),]$NEO1_level == "High")  # 13/23

# Choosing the high/low limit with KMeans
km.out <- kmeans(basic_coldata[which(basic_coldata$Time == 'pre'),'NEO1_log'], centers = 2)
km.out # valide la limite entre 2.42 et 2.92
# Rq: il faut retirer le sample 06-006, qui a une QC très mauvaise et qui n'exprime pas NEO1 alors qu'il est classé PR
median(basic_coldata[which(basic_coldata$Time == 'pre'),'NEO1_log']) # 3.08



### Survival plots: OS
basic_coldata_pre <- basic_coldata[which(basic_coldata$Time == 'pre' & basic_coldata$SampleID != '06-006'),]
kaplan_meier <- with(basic_coldata_pre, Surv(OS, Death_status))
# summary(data.frame(kaplan_meier)$futime) # marche pas ...
km_fit_os_all <- survfit(Surv(OS, Death_status) ~ 1, data=basic_coldata_pre)
# median: 509

km_fit_os <- survfit(Surv(OS, Death_status) ~ NEO1_level, data=basic_coldata_pre)
# median High: 621, median Low: 501
summary(km_fit_os, times = c(1,50,100*(1:10)))
print(km_fit_os, print.rmean = TRUE) # restricted mean survival time (RMST)
autoplot(km_fit_os) #conf.int=FALSE

cox_os <- coxph(Surv(OS, Death_status) ~ NEO1_level, data=basic_coldata_pre)
summary(cox_os)
cox_fit_os <- survfit(cox_os)
autoplot(cox_fit_os)

survdiff(Surv(OS, Death_status) ~ NEO1_level, data=basic_coldata_pre)

survminer::ggsurvplot(
  km_fit_os, 
  data = basic_coldata_pre,          # again specify the data used to fit linelistsurv_fit_sex 
  conf.int = TRUE,              # do not show confidence interval of KM estimates
  surv.scale = "percent",        # present probabilities in the y axis in %
  break.time.by = 50,            # present the time axis with an increment of 10 days
  xlab = "Follow-up days",
  ylab = "Survival Probability",
  pval = T,                      # print p-value of Log-rank test 
  pval.coord = c(40,.21),        # print p-value at these plot coordinates
  risk.table = T,                # print the risk table at bottom 
  legend.title = "Neogenin level",       # legend characteristics
  legend.labs = c("High","Low"),
  font.legend = 10, 
  palette = "Dark2",             # specify color palette 
  surv.median.line = "hv",       # draw horizontal and vertical lines to the median survivals
  ggtheme = theme_light()        # simplify plot background
)

# NEO1_log en tant que variable linéaire: mauvaises p-values, ne vaut pas le coup
coxph(Surv(OS, Death_status) ~ NEO1_log, data=basic_coldata_pre)
coxph(Surv(PFS, Prog_status) ~ NEO1_log, data=basic_coldata_pre)
coxph(Surv(OS, Death_status) ~ NTN1_log, data=basic_coldata_pre)
coxph(Surv(PFS, Prog_status) ~ NTN1_log, data=basic_coldata_pre)

# comparaison 2 cat survie
km_fit_os_bestres <- survfit(Surv(OS, Death_status) ~ Best_response_2, data=basic_coldata_pre)
# median High: 621, median Low: 501
summary(km_fit_os_bestres, times = c(1,50,100*(1:10)))
print(km_fit_os_bestres, print.rmean = TRUE) # restricted mean survival time (RMST)
autoplot(km_fit_os_bestres) #conf.int=FALSE

# comparaison cats target %
ggplot(basic_coldata[which(basic_coldata$Time == 'pre'),], 
       aes(x=best_target_prec)) + 
  geom_histogram()
km.out <- kmeans(basic_coldata[which(basic_coldata$Time == 'pre'),'best_target_prec'], centers = 2)
km.out # séparation à -32
basic_coldata_pre$best_target_level <- case_when(
  basic_coldata_pre$best_target_prec >= -32 ~ "Normal",
  basic_coldata_pre$best_target_prec < -32 ~ "Shrinked")
km_fit_os_bestres2 <- survfit(Surv(OS, Death_status) ~ best_target_level, data=basic_coldata_pre)
autoplot(km_fit_os_bestres2) #conf.int=FALSE



### Survival plots: PFS
km_fit_pfs_all <- survfit(Surv(PFS, Prog_status) ~ 1, data=basic_coldata_pre)
# median: 399
km_fit_pfs <- survfit(Surv(PFS, Prog_status) ~ NEO1_level, data=basic_coldata_pre)
# median High: 331, median Low: 442

summary(km_fit_pfs, times = c(1,50,100*(1:10)))
print(km_fit_pfs, print.rmean = TRUE) # restricted mean survival time (RMST)
autoplot(km_fit_pfs) #conf.int=FALSE

cox_pfs <- coxph(Surv(PFS, Prog_status) ~ NEO1_level, data=basic_coldata_pre)
summary(cox_pfs)
cox_fit_pfs <- survfit(cox_pfs)
autoplot(cox_fit_pfs)

survdiff(Surv(PFS, Prog_status) ~ NEO1_level, data=basic_coldata_pre)

survminer::ggsurvplot(
  km_fit_pfs, 
  data = basic_coldata_pre,          # again specify the data used to fit linelistsurv_fit_sex 
  conf.int = TRUE,              # do not show confidence interval of KM estimates
  surv.scale = "percent",        # present probabilities in the y axis in %
  break.time.by = 50,            # present the time axis with an increment of 10 days
  xlab = "Follow-up days",
  ylab = "Non-Progression Probability",
  pval = T,                      # print p-value of Log-rank test 
  pval.coord = c(40,.21),        # print p-value at these plot coordinates
  risk.table = T,                # print the risk table at bottom 
  legend.title = "Neogenin level",       # legend characteristics
  legend.labs = c("High","Low"),
  font.legend = 10, 
  palette = "Dark2",             # specify color palette 
  surv.median.line = "hv",       # draw horizontal and vertical lines to the median survivals
  ggtheme = theme_light()        # simplify plot background
)



### ssGSEA pre
tmp <- read.csv('/home/eugenie-modolo/Documents/Reference_files/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2025.1.Hs.tsv', sep="\t")
list_genes <- unlist(strsplit(tmp[17,2], ","))
# expression_genes <- count_matrix[which(row.names(count_matrix) %in% list_genes),]
# /!\ 3 gènes présents en doublon : COL3A1, GREM1, PRSS2, dont les expressions varient beaucoup entre les deux gènes identiques ...
# COL3A1 : ENSG00000168542.19 et ENSG00000291748.1 -> le deuxième ne ressort pas dans les convertisseurs (2 sens)
# GREM1 : ENSG00000166923.13 et ENSG00000276886.4 -> le deuxième ne ressort pas dans les convertisseurs (2 sens)
# PRSS2 : ENSG00000275896.7 et ENSG00000282049.1 -> le deuxième ne ressort pas dans les convertisseurs (2 sens)
counts_genes[which(row.names(counts_genes)=='ENSG00000291748.1'),'gene_name'] <- 'COL3A1_2' # absent de counts_genes$gene_name
counts_genes[which(row.names(counts_genes)=='ENSG00000276886.4'),'gene_name'] <- 'GREM1_2' # absent de counts_genes$gene_name
counts_genes[which(row.names(counts_genes)=='ENSG00000282049.1'),'gene_name'] <- 'PRSS2_2' # absent de counts_genes$gene_name
count_columns <- colnames(counts_genes)[colnames(counts_genes) != "gene_name"]
count_matrix <- as.matrix(counts_genes[, count_columns])
rownames(count_matrix) <- counts_genes$gene_name

# code pour retirer les noms de gènes dupliqués, en les remplaçant par leur ENSG code
dup_genes <- counts_genes[duplicated(counts_genes$gene_name),'gene_name']
counts_genes$feature <- ifelse((is.na(counts_genes$gene_name) | counts_genes$gene_name %in% dup_genes), rownames(counts_genes), counts_genes$gene_name)
count_columns <- colnames(counts_genes)[colnames(counts_genes) != "gene_name" & colnames(counts_genes) != "feature"]
count_matrix <- as.matrix(counts_genes[, count_columns])
rownames(count_matrix) <- counts_genes$feature
count_matrix_pre <- count_matrix[,c("01-001_pre","01-003_pre","01-004_pre",
                        "01-006_pre","01-007_pre","01-011_pre","01-012_pre","02-001_pre","02-002_pre","02-003_pre",
                        "02-006_pre","02-010_pre","02-014_pre","03-002_pre","03-003_pre","03-006_pre","06-001_pre",
                        "06-007_pre","06-008_pre","06-009_pre","08-001_pre","08-002_pre")]

# ssGSEA
list_genes_ssgsea <- list(list_genes, "")
names(list_genes_ssgsea) <- c("EMT", "")
#res <- igsva() # pour avoir l'app visuelle
gsvaparams <- gsvaParam(count_matrix_pre,
                        list_genes_ssgsea)
result_gsva <- gsva(gsvaparams)#, verbose = TRUE)
scores_EMT <- result_gsva[1,]
basic_coldata_pre <- basic_coldata[which(basic_coldata$Time == 'pre' & basic_coldata$SampleID != '06-006'),]
# verif que les échantillons sont dans le bon ordre :
# sum(rownames(basic_coldata_pre) == names(scores_EMT)) == length(rownames(basic_coldata_pre))
basic_coldata_pre$EMT_score <- scores_EMT
#km.out <- kmeans(basic_coldata_pre$EMT_score, centers = 2)
#km.out # il coupe à 0
basic_coldata_pre$EMT_s_level <- case_when(
  basic_coldata_pre$EMT_score >= 0 ~ "High",
  basic_coldata_pre$EMT_score < 0 ~ "Low")

# Heatmap of the scores
library(gplots)
expression_genes <- count_matrix_pre[which(row.names(count_matrix_pre) %in% list_genes),]
distcor <- function(m){
  as.dist(1-cor(t(m))  )
  }
#heatmap.2(as.matrix(log2(expression_genes+1)), distfun = distcor, hclustfun = function(d)hclust(d,method="ward.D2"), scale="none", trace="none", dendrogram="col", Rowv = F,  col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(50))
heatmap.2(as.matrix(log2(expression_genes+1)), distfun = distcor  , hclustfun = function(d)hclust(d,method="ward.D2"), scale="row", trace="none", dendrogram="both",  col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(50))
target <- c("02-006_pre", "06-009_pre", "06-008_pre", "01-003_pre", 
            "06-007_pre", "03-002_pre", "03-006_pre", "02-002_pre", "01-012_pre", "02-003_pre", "06-001_pre", 
            "08-002_pre", "01-006_pre", "01-007_pre", "01-004_pre", "03-003_pre", "08-001_pre", "01-011_pre", 
            "02-001_pre", "02-010_pre", "01-001_pre", "02-014_pre")
basic_coldata_pre$row.names <- rownames(basic_coldata_pre)
basic_coldata_pre %>% mutate(row.names=fct_relevel(row.names, target)) %>% ggplot(aes(x=row.names, y=EMT_score)) +
  geom_point() +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1))

# Kaplan-Meier maps according to scores
km_fit_os <- survfit(Surv(OS, Death_status) ~ EMT_s_level, data=basic_coldata_pre)
autoplot(km_fit_os)
km_fit_pfs <- survfit(Surv(PFS, Prog_status) ~ EMT_s_level, data=basic_coldata_pre)
autoplot(km_fit_pfs)



### ssGSEA pre-post
target <- c("01-007_pre","01-007_post","02-001_pre","02-001_post","02-003_pre",
            "02-003_post","02-014_pre","02-014_post")
basic_coldata_prepost <- basic_coldata[which(rownames(basic_coldata) %in% target),]
count_matrix_prepost <- count_matrix[,target]
gsvaparams <- gsvaParam(count_matrix_prepost,
                        list_genes_ssgsea)
result_gsva <- gsva(gsvaparams)#, verbose = TRUE)
scores_EMT <- result_gsva[1,]
basic_coldata_prepost$EMT_score <- scores_EMT
hist(basic_coldata_prepost$EMT_score)




### PCA & ACM of all samples
#dim(count_matrix) # 82059    29
count_matrix_filt <- count_matrix[rowSums(count_matrix) != 0,]
#dim(count_matrix_filt) # 59229    29

## log2(count+1)
data_f_pca <- prcomp(t(log2(count_matrix_filt + 1)))
# fviz_pca_ind(data_f_pca,
#              geom = c("point","text"),
#              geom.ind = basic_coldata$Best_response_2,
#              col.ind = basic_coldata$Time
# )
pca_data <- data.frame(
  PC1 = data_f_pca$x[, 1],
  PC2 = data_f_pca$x[, 2]
)
ggplot(pca_data, aes(x = PC1, y = PC2, 
                     color = basic_coldata$Best_response_2,  
                     shape = basic_coldata$Time)) + 
  geom_point(size = 3) +
  geom_text(aes(label = basic_coldata$SampleID), 
            vjust = -1.2, hjust = 0.5, size = 3, 
            show.legend = FALSE)
  scale_shape_manual(values = c(16, 17)) +  # Circle and triangle; adjust as needed
  labs(title = "PCA Plot",
       x = paste0("PC1 (", round(summary(data_f_pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(data_f_pca)$importance[2,2]*100, 1), "%)")) +
  theme_minimal() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
# Scree plot
pc_var <- data_f_pca$sdev^2 / (ncol(data_f_pca$x) - 1)
tibble(
  pc = 1:29,
  var = pc_var / sum(pc_var)
) %>%
  ggplot() +
  geom_point(aes(x = pc, y = var))
# Variables contribution
# fviz_pca_biplot(data_f_pca) # marche pas car noms de gènes sont dupliqués
count_columns <- colnames(counts_genes)[colnames(counts_genes) != "gene_name"]
count_matrix_ENSG <- as.matrix(counts_genes[, count_columns])
count_matrix_ENSG_filt <- count_matrix_ENSG[rowSums(count_matrix_ENSG) != 0,]
data_f_pca2 <- prcomp(t(log2(count_matrix_ENSG_filt + 1)))
fviz_pca_biplot(data_f_pca2)
# ENSG00000276788.1	SNORD26 Small Nucleolar RNA
# ENSG00000168925.12	CTRB1 a member of the serine protease family of enzymes and forms a principal precursor of the pancreatic proteolytic enzymes


## VST
count_matrix_vst <- vst(round(count_matrix_filt))


## anscombe
count_matrix_ans <- 2*sqrt(count_matrix_filt + 3/8)



### Export results
write.csv(as.data.frame(xx), 
          file="xx.csv")


#### Ancient code

