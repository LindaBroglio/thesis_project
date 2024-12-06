
library(DESeq2)
library(plyr)
library(dplyr)

run_DESeq_analysis <- function(counts_path, metadata_path, factor_col, factor_levels) {
  # load data
  counts <- read.table(counts_path, header = TRUE, row.names = 1, sep = "\t")
  colData <- read.table(metadata_path, header = TRUE, row.names = 1, sep = "\t")
  
  # filter colData, rows where the factor column has one of input levels
  colData <- colData[colData[[factor_col]] %in% factor_levels, ]
  colData <- colData[factor_col]
  
  # filter counts to match colData
  counts <- counts[, rownames(colData)]
  
  # convert the factor column to a factor with input levels
  colData[[factor_col]] <- factor(colData[[factor_col]], levels = factor_levels)
  
  print(head(counts))
  print(colData)

  # create DESeq dataset
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = as.formula(paste("~", factor_col)))
  
  # run DESeq analysis
  dds <- DESeq(dds)
  res <- results(dds)

  return(list(
  norm_counts = counts(dds, normalized = TRUE), 
  res = results(dds)
  ))
}


produce_and_save_table <- function(deseq_output, output_path) {
  # get results from the analysis output
  norm_counts <- deseq_output$norm_counts
  res <- deseq_output$res

  # CPM normalized counts
  cpm <- t(t(norm_counts) / colSums(norm_counts)) * 1e6
  
  # combine CPM with DESeq2 results 
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df) 
  
  cpm_df <- as.data.frame(cpm)
  cpm_df$gene <- rownames(cpm_df)

  combined_table <- join(cpm_df, res_df)
  final_table <- combined_table %>% select(gene, everything())
  
  # save final table as a tsv
  write.table(final_table, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(final_table)
}


deseq_output <- run_DESeq_analysis("/home/lbroglio/rna_seq/gene_counts_table.tsv",
                                   "/home/lbroglio/rna_seq/metadata.tsv",
                                   "status_pro",
                                   c("healthy_twin", "affected_twin"))


produce_and_save_table(deseq_output, "/home/lbroglio/rna_seq/R_output_files/h_twin_vs_a_twin.tsv")


# summary(res)
# png("/home/lbroglio/rna_seq/R_output_files/plots/MA_plot.png", width = 800, height = 600)
# plotMA(res, ylim = c(-2, 2), cex = 0.5, colSig = "red")
# top_genes <- rownames(res)[which(res$padj < 0.1)[1:10]]
# with(res[top_genes, ], text(baseMean, log2FoldChange, labels = top_genes, pos = 4, cex = 0.7, col = "blue"))
# dev.off()
# res_sorted <- res[order(res$padj), ]
# top_genes <- head(res_sorted, 10)
# print(top_genes)

