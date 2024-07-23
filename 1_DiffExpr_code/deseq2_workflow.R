# load libraries
#BiocManager::install("DESeq2")
# install.packages("xlsx")

library(DESeq2)
library(ggplot2)
library(xlsx)


###############
# Functions
###############

# function to read and do quick processing of 
readCountData <- function(filename, sampTab){
  count_data <- read.csv(filename, sep = "\t", header = TRUE)
  rownames(count_data) <- count_data$Genes
  count_data <- count_data[,-1]
  colnames(count_data) <- rownames(sampleTable)
  
  return(count_data)
}

# a function to generate a table of gene hits with statistics and normalized counts
getDiffGenesTable <- function(dataset, alpha = 0.05,lfcThreshold = 1) {
  # perform filtering of the data
  res.05 <- results(dataset, alpha = alpha, lfcThreshold = lfcThreshold)
  
  # obtain normalized counts  
  dds <- estimateSizeFactors(dataset)
  norm_counts <- counts(dataset, normalized=TRUE)

  # generating a final table of regulated genes
  resSig <- subset(res.05, padj < alpha)
  resSig_full <- resSig[ order(-resSig$log2FoldChange), ]
  resSig_full$Fold_change <- 2^resSig_full$log2FoldChange
  
  norm_counts_sig <- norm_counts[rownames(resSig_full),]

  # combine significant hits table and normalized counts matrix
  output <- cbind(resSig_full, norm_counts_sig)
  
  # output the results
  return(output)
}

##########################
# End Functions definition
##########################

# load and process the sample table - custom for each table
sampleTable <- read.csv("sampleData.txt", sep = "\t")
rownames(sampleTable) <- sampleTable$sampleName 
rownames(sampleTable)
sampleTable$genotype <- factor(sampleTable$genotype, levels = c("wildtype","miR34a mutant"))
sampleTable$treatment <- factor(sampleTable$treatment, levels = c("DMSO","CPT"))
str(sampleTable)
sampleTable

# load the raw count data
count_data <- readCountData("counts.txt", sampleTable)

# add gene symbols to the count data table and 
# convert v4.3.2 IDs to gene symbols
ID_symbol_map <- read.csv("v4_3_2geneinfo.txt", sep = '\t')
row.names(ID_symbol_map) <- ID_symbol_map$LLgeneID

count_data_symbols <- count_data

count_data_symbols$symbol <- ID_symbol_map[row.names(count_data_symbols), ]$LLgeneSymbol

count_data_symbols <- count_data_symbols[, c(13,1:12)]
head(count_data_symbols)

write.table(count_data_symbols, "Zebrafish_28h_cpt_experiment.txt", sep = '\t')

# convert the estimates to counts 
saved_rn <- rownames(count_data)
count_data <- as.data.frame(sapply(count_data, ceiling))
rownames(count_data) <- saved_rn


## Full dataset 

# 1. Read in the count data into a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data, 
                              colData = sampleTable, 
                              design = ~ treatment)

# Pre-filtering the data

# check the number of rows before filtering
nrow(dds)
# [1] 35197

dds <- dds[rowSums(counts(dds) == 0) < 3, ]
nrow(dds)
# [1] 25182

dds <- dds[rowSums(counts(dds)) > 24, ]
nrow(dds)
# [1] 25063

dds <- dds[rowSums(counts(dds) >= 10) >= 3,]
nrow(dds)
# [1] 23060

# 2. Perform vst transform to scale and stabilize variance 
vsd <- vst(dds)


# 3.2 Plot MDS
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(vsd)))
ggplot(mds, aes(X1,X2,color=treatment,shape=genotype)) + geom_point(size=3)

ggsave("all_data_MDS-plot.png", dpi = 300)

# 4. Perform differential expression analysis 
dds <- DESeq(dds)

# compile the results into a single table
complete_results_table <- getDiffGenesTable(dds, alpha = 0.05, lfcThreshold = 1)
nrow(complete_results_table)
# [1] 2894

# add a symbol column
complete_results_table$symbol <- ID_symbol_map[row.names(complete_results_table), ]$LLgeneSymbol
ncol(complete_results_table)

complete_results_table <- complete_results_table[, c(20,1:19)]
head(complete_results_table)


write.table(complete_results_table, "regulated-genes_cpt_dataset.tsv", sep = "\t")

# check the number of up-regulated genes
nrow(complete_results_table[complete_results_table$log2FoldChange > 0,])
# [1] 671

# check the number of down-regulated genes
nrow(complete_results_table[complete_results_table$log2FoldChange < 0,])
# [1] 2223



################################## MA plot #####################################

#install.packages("ggpubr")

library("apeglm")
resultsNames(dds)

resultsMA <- results(dds, lfcThreshold = 1, alpha = 0.01)

library(ggpubr)

ggmaplot(resultsMA, 
         fdr = 0.05, fc = 2, size = 0.8,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         legend = "top", top = 80,
         genenames = ID_symbol_map[rownames(resultsMA), ]$LLgeneSymbol,
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = c("bold", 14),
         font.main = c("bold", 12),
         ggtheme = ggplot2::theme_minimal()) + theme(text = element_text(size=14))

ggsave("MAplot.png", dpi = 600)

# Heatmaps
library(genefilter)
library(pheatmap)


# A heatmap of ALL regulated genes
vsd_all <- vsd[row.names(complete_results_table),]
mat_all <- assay(vsd_all)
mat_all <- mat_all - rowMeans(mat_all)
anno <- as.data.frame(colData(vsd_all)[, c("genotype","treatment")])

pheatmap(mat_all,
         method="complete", 
         annotation_col = anno, show_rownames = F, 
         annotation_legend = TRUE, legend=T, 
         cluster_cols=TRUE,scale="row", 
         treeheight_row = 0, treeheight_col = 0,
         color=colorRampPalette(c("navy", "white", "red"))(50))

