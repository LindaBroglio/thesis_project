library(DESeq2)
library(ggplot2)

#source("/home/lbroglio/rna_seq/R_output_files/R_model_comparison/p_adj_extract.R")

perform_deseq2_analysis <- function(counts_file, metadata_file) {
  
  counts_data <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t")
  metadata <- read.table(metadata_file, header = TRUE, row.names = 1, sep = "\t")
  print(dim(counts_data))
  
  # filters
  counts_data <- counts_data[rowSums(counts_data) > 0.0, ]
  counts_data <- counts_data[rowSums(counts_data >= 10) >= 3, ]
  print(dim(counts_data))

  # sample_sums <- colSums(counts_data)
  # print(sample_sums)

  metadata_filtered <- metadata[metadata$status_pro %in% c("healthy_twin", "affected_twin"), ]
  counts_data_filtered <- counts_data[, rownames(metadata_filtered)]
  
  dds <- DESeqDataSetFromMatrix(countData = counts_data_filtered,
                                colData = metadata_filtered,
                                design = ~ status_pro + family)
  
  dds <- DESeq(dds)

  return(dds)
}

library(dplyr)

produce_and_save_table <- function(dds, output_path, order_by = NULL) {
  
  norm_counts <- counts(dds, normalized = TRUE)
  
  # cpm <- t(t(norm_counts) / colSums(norm_counts)) * 1e6
  res <- results(dds, contrast = c("status_pro", "affected_twin", "healthy_twin"))
  
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df) 
  
  # cpm_df <- as.data.frame(cpm)
  # cpm_df$gene <- rownames(cpm_df)
  # #CPM and DESeq2 results 
  # combined_table <- inner_join(cpm_df, res_df, by = "gene")

  norm_counts_df <- as.data.frame(norm_counts)
  norm_counts_df$gene <- rownames(norm_counts_df)
  
  combined_table <- inner_join(norm_counts_df, res_df, by = "gene")

  # Final table with gene column at the front
  final_table <- combined_table %>% select(gene, everything())
  final_table <- final_table %>% mutate_if(is.numeric, ~ round(., 5)) 

  write.table(final_table, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(final_table)
}

counts_file <- "/home/lbroglio/rna_seq/gene_counts_table.tsv"  
metadata_file <- "/home/lbroglio/rna_seq/metadata.tsv"   
output_path <- "/home/lbroglio/rna_seq/DGE/gene_centric_normcounts_and_stats.tsv"

de_results <- perform_deseq2_analysis(counts_file, metadata_file)

final_table <- produce_and_save_table(de_results, output_path)

# extract_top_padj(
#   input_path = output_path,
#   output_path = output_path_table_padj
# )

