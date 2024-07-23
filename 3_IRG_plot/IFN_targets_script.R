library("edgeR")
library(DESeq2)
library(ggplot2)
library(xlsx)
library(dplyr)
# 

# ###############
# # Functions for edgeR
# ###############

## function to read and do quick processing of count data

readCountData <- function(filename, sampTab){
  count_data <- read.csv(filename, sep = "\t", header = TRUE)
  rownames(count_data) <- count_data$Genes
  count_data <- count_data[,-1]
  colnames(count_data) <- rownames(sampleTable)
  
  return(count_data)
}
# 
## a function to generate a table of gene hits with statistics and normalized counts
getDiffGenesTable_edgeR <- function(results, dataset, alpha = 0.05, lfcThreshold = 1) {
  
  # perform filtering of the data
  res <- as.data.frame(results[results$table$FDR < alpha & 
                                 abs(results$table$logFC) >= lfcThreshold, ])
  
  # obtain normalized counts
  norm_counts <- dataset$counts/dataset$samples$norm.factors
  
  # generating a final table of regulated genes
  resSig_full <- res[ order(-res$logFC), ]
  resSig_full$Fold_change <- 2^resSig_full$logFC
  
  norm_counts_sig <- norm_counts[rownames(resSig_full),]
  
  # combine significant hits table and normalized counts matrix
  output <- cbind(resSig_full, norm_counts_sig)
  
  # output the results
  return(output)
  
}


# function to perform all steps of analysis using edgeR
performAnalysis <- function(counts, samples, genes, design, alpha, lfcThreshold){
  
  # creating edgeR analysis object and performing some analysis steps
  edge_object <- DGEList(counts=counts, samples=samples, genes=genes)
  edge_object <- calcNormFactors(edge_object)
  edge_object <- estimateDisp(edge_object, design)
  
  # getting results
  fit <- glmFit(edge_object, design)
  lrt <- glmLRT(fit, coef=ncol(design))
  
  tt.all <- topTags(lrt, n=nrow(edge_object), sort.by="none")
  
  # getDiffGenesTable_edgeR <- function(results, dataset, alpha = 0.05, lfcThreshold = 1)
  results <- getDiffGenesTable_edgeR(tt.all, edge_object, alpha = alpha, lfcThreshold = lfcThreshold)
  
  return(results)
  
}


# ##########################
# # End Functions definition
# ##########################

############################
# Reading the data
############################


sampleTable <- read.csv("sampleData.txt", sep = "\t")
rownames(sampleTable) <- sampleTable$sampleName 
rownames(sampleTable)
sampleTable$genotype <- factor(sampleTable$genotype, levels = c("wildtype","mir34a_deletion"))
sampleTable$treatment <- factor(sampleTable$treatment, levels = c("DMSO","CPT"))
str(sampleTable)
sampleTable

# load the raw count data
count_data <- readCountData("counts.txt", sampleTable)
nrow(count_data)
# [1] 35197

# filter counts with 0s 
count_data <- count_data[rowSums(count_data == 0) < 4, ]
count_data <- count_data[rowSums(count_data) > 24,]
nrow(count_data)
# [1] 25626

#######################################
## Full dataset
#######################################


# edgeR differential analysis

genetable <- data.frame(gene = rownames(count_data))
edge_full <- DGEList(counts=count_data, samples=sampleTable, genes=genetable)
names(edge_full)

design_full <- model.matrix(~0 + genotype + treatment, edge_full$samples)

edge_full <- calcNormFactors(edge_full)
edge_full <- estimateDisp(edge_full, design_full)

# getting results
fit <- glmFit(edge_full, design_full)
lrt <- glmLRT(fit, coef=ncol(design_full))

tt.all <- topTags(lrt, n=nrow(edge_full), sort.by="none")

# identify final DEGs from this dataset
results_full <- getDiffGenesTable_edgeR(tt.all, edge_full, alpha = 0.05, lfcThreshold = 1)
dim(results_full)


# get IFN target genes

ifn_genes <- read.csv("ifn_genes.csv")
ifn_genes <- as.vector(ifn_genes$Gene_name)


# add gene symbols to the count data table and 
# convert v4.3.2 IDs to gene symbols
ID_symbol_map <- read.csv("v4_3_2geneinfo.txt", sep = '\t')
row.names(ID_symbol_map) <- ID_symbol_map$LLgeneID

results_full$symbol <- ID_symbol_map[row.names(results_full), ]$LLgeneSymbol

results_full <- results_full[, c(20,1:19)]
head(results_full)


# 2. subset this list based on whether a gene was differentially expressed in the full dataset
# under conditions of fold change > 2 and P-value < 0.05

IFN_regulated_inData <- intersect(results_full$symbol, ifn_genes)
IFN_regulated_inData

write.xlsx(results_full[results_full$symbol %in% IFN_regulated_inData,], "regulated_IFNtargets_data.xlsx", sheetName="Sheet1", col.names=TRUE, row.names=TRUE, append=FALSE)

################## number of IFN genes in the whole dataset ###################

count_data$symbol <- ID_symbol_map[row.names(count_data), ]$LLgeneSymbol

count_data <- count_data[, c(13,1:12)]
head(count_data)

IFN_regulated_all <- intersect(count_data$symbol, ifn_genes)
IFN_regulated_all


########################################
# Heatmap for the full dataset genes
########################################

# Heatmaps
library(genefilter)
library(pheatmap)

# convert the estimates to counts 
saved_rn <- rownames(count_data)
count_data <- as.data.frame(sapply(count_data, ceiling))
rownames(count_data) <- saved_rn

dds_full <- DESeqDataSetFromMatrix(countData = count_data, colData = sampleTable, design = ~ genotype)
vsd_full <- vst(dds_full)

# heatmap for IFN target genes
mat  <- assay(vsd_full)
mat  <- mat - rowMeans(mat)
mat <-as.data.frame(mat)

mat$symbol <- ID_symbol_map[row.names(mat), ]$LLgeneSymbol

mat <- mat[mat$symbol %in% IFN_regulated_inData,]
rownames(mat) <- mat$symbol

mat <- as.matrix(mat[,c(1:12)])

anno <- as.data.frame(colData(vsd_full)[, c("genotype","treatment")])

pheatmap(mat, method="complete", 
         annotation_col = anno, show_rownames = T, 
         annotation_legend = TRUE, legend=T, 
         cluster_cols=TRUE,
         treeheight_row = 0, treeheight_col = 0, fontsize = 8,
         color=colorRampPalette(c("navy", "white", "red"))(50))


# average columns
mat_every3 <- t(apply(mat, 1, tapply, gl(4, 3), mean))

colnames(mat_every3) <- c("wt - DMSO", "wt - CPT", "miR-34a_del - DMSO", "miR-34a_del - CPT" )
colnames(mat_every3)


annoTable <- read.csv("anno_table.csv", sep = "\t")
rownames(annoTable) <- annoTable$name
rownames(annoTable)
annoTable$genotype <- factor(annoTable$genotype, levels = c("wildtype","mir34a_deletion"))
annoTable$treatment <- factor(annoTable$treatment, levels = c( "DMSO", "CPT"))
anno <- annoTable[, c( "treatment", "genotype")]

anno <- anno %>% arrange(treatment)

mat_every3 <- mat_every3[order(mat_every3[,2], decreasing = TRUE),]

# splits used for the figure
pheatmap(mat_every3[1:55,c("wt - DMSO", "miR-34a_del - DMSO",   "wt - CPT",  "miR-34a_del - CPT")], 
         cellwidth = 25, cellheight = 10, treeheight_col = 0, scale = 'row',
         cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, annotation_col = anno, fontsize = 8,
         color=colorRampPalette(c("navy", "white", "red"))(50))

pheatmap(mat_every3[56:114,c("wt - DMSO", "miR-34a_del - DMSO",   "wt - CPT",  "miR-34a_del - CPT")], 
         cellwidth = 25, cellheight = 10, treeheight_col = 0, scale = 'row',
         cluster_cols = FALSE, cluster_rows = FALSE, treeheight_row = 0, annotation_col = anno, fontsize = 8,
         color=colorRampPalette(c("navy", "white", "red"))(50))



