
gene_FC_filter <- function(input) {
#input <- sel_filter_data

iterate_gene_set <- function(input, pathway_ind, ind)  { 
  
  sel_ind <- pathway_ind[-ind]
    expr <- lazyeval::interp(quote(x == y), x = as.name(pathway_ind[ind]), y = 1)
  input %>%
  filter_(expr) %>%
  group_by(CELL_ID, pathway_ind[ind], GeneSymbol ) %>%
  summarise( mean_l2FC = max(abs(log2FoldChange[is.finite(log2FoldChange)]), na.rm = TRUE)) %>%
  top_n(n = 30, wt = abs(mean_l2FC))
  
}

  pathway_ind <- c("OX", "INFLAM")#, "DDR", "UPR")

  sel_genes <- lapply(1:length(pathway_ind),  function(i) {
    iterate_gene_set(input = sel_filter_data, pathway_ind=pathway_ind, ind = i)
    }
    )

  sel_genes <- do.call('rbind', sel_genes)
  sel_genes <- unique(sel_genes$GeneSymbol)
  
input %>% filter(GeneSymbol %in% sel_genes)

}




