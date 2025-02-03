library(biomaRt)

data <- read.csv("/home/lbroglio/rna_seq/DGE/gene_centric_normcounts_and_stats.tsv", header = TRUE, sep = "\t")

mart <- useEnsembl(biomart = 'genes', 
                   dataset = 'hsapiens_gene_ensembl',
                   version = 110)

annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "description"), 
  filters = "ensembl_gene_id",
  values = data$gene,
  mart = mart
)

data_with_annotations <- merge(data, annotations, by.x = "gene", by.y = "ensembl_gene_id", all.x = TRUE)

data_with_annotations <- data_with_annotations[match(data$gene, data_with_annotations$gene), ]

data_with_annotations <- data_with_annotations[, c(names(data_with_annotations)[1], "external_gene_name", 
                                                   names(data_with_annotations)[-1][names(data_with_annotations)[-1] != "external_gene_name"])]

output_path <- "/home/lbroglio/rna_seq/DGE/gene_centric_normcounts_stats_names.tsv"
write.table(data_with_annotations, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Annotated table saved to:", output_path, "\n")
