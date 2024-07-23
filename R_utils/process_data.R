library(readr)
library(stringr)


# get all files
base_dir <- getwd()
results_genes <- list()
dir_files <- list.files(base_dir)

# get gene files
gene_files <- dir_files[endsWith(dir_files, "2nd_pass_.genes.results")]

for(fname in gene_files){
  # read data from a file
  data <- read_tsv(fname, col_names = TRUE)
  
  # make sample name
  name_parts <- unlist(strsplit(fname,  '_2nd_'))
  sample_name <- name_parts[1]
  
  # store the data
  results_genes[[ sample_name ]] = data$expected_count
}


results_genes_df <- as.data.frame(results_genes)
row.names(results_genes_df) <- data$gene_id

# correct the names
colnames(results_genes_df) <- sapply(as.vector(colnames(results_genes_df)), function(x){str_replace(x, "mir34a", "miR34a" )}, USE.NAMES = FALSE)

# reorder columns
results_genes_df <- results_genes_df[, c("wt_DMSO_1", "wt_DMSO_2", "wt_DMSO_3", "wt_CPT_1", "wt_CPT_2", "wt_CPT_3", "miR34a_del_DMSO_1", "miR34a_del_DMSO_2", "miR34a_del_DMSO_3", "miR34a_del_CPT_1", "miR34a_del_CPT_2",  "miR34a_del_CPT_3")]


head(results_genes_df)


write.csv(results_genes_df, "28hpf-miR34a_del_STAR-RSEM_genes_results.csv")


# get isoform files
results_isoforms <- list()

isoform_files <- dir_files[endsWith(dir_files, "2nd_pass_.isoforms.results")]

for(fname in isoform_files){
  # read data from a file
  data <- read_tsv(fname, col_names = TRUE)
  
  # make sample name
  name_parts <- unlist(strsplit(fname,  '_2nd_'))
  sample_name <- name_parts[1]
  
  # store the data
  results_isoforms[[ sample_name ]] = data$expected_count
}


results_isoforms_df <- as.data.frame(results_isoforms)
row.names(results_isoforms_df) <- data$transcript_id

# correct the names
colnames(results_isoforms_df) <- sapply(as.vector(colnames(results_isoforms_df)), function(x){str_replace(x, "mir34a", "miR34a" )}, USE.NAMES = FALSE)

# reorder columns
results_isoforms_df <- results_isoforms_df[, c("wt_DMSO_1", "wt_DMSO_2", "wt_DMSO_3", "wt_CPT_1", "wt_CPT_2", "wt_CPT_3", "miR34a_del_DMSO_1", "miR34a_del_DMSO_2", "miR34a_del_DMSO_3", "miR34a_del_CPT_1", "miR34a_del_CPT_2",  "miR34a_del_CPT_3")]


head(results_isoforms_df)

write.csv(results_isoforms_df, "28hpf-miR34a_del_STAR-RSEM_isoforms_results.csv")

