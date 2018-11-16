
meta_filter <- function(input, meta_file_path) {

  
  metadata <- read.table(file = meta_file_path, sep = "\t", header= TRUE, row.names = 1)
  metadata <- metadata[, c("TNF", "DEM")]
  colnames(metadata) <- c("INFLAM", "OX")

#metadata$UPR <- sapply(sapply(metadata$UPR, strsplit, "_"), "[[",2)
metadata$OX <- sapply(sapply(metadata$OX, strsplit, "_"), "[[",2)
#metadata$DDR <- sapply(sapply(metadata$DDR, strsplit, "_"), "[[",2)
metadata$INFLAM <- sapply(sapply(metadata$INFLAM, strsplit, "_"), "[[",2)

all_genes <- metadata %>%
  gather(key = pathway, value = gene) 


gene_metadata <- data.frame( gene = unique(all_genes$gene) ) # unique set of genes as not to multiply data in the join operation



gene_metadata <- gene_metadata %>% mutate(OX = if_else(gene %in% metadata$OX, 1, 0),
                         INFLAM = if_else(gene %in% metadata$INFLAM, 1, 0)#,
                         #DDR = if_else(gene %in% metadata$DDR, 1, 0),
                         #UPR = if_else(gene %in% metadata$UPR, 1, 0)
                         )



data_annotated <-left_join(input, gene_metadata, by = c("GeneSymbol" = "gene"))



data_annotated %>%
    filter( !is.na(OX)   )

}



