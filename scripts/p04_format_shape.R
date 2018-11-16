

probe_summarise_format_data <- function(input) {

 input %>% 
  select(log2FoldChange, GeneSymbol, CELL_ID, TREATMENT, CONCENTRATION, OX, INFLAM) %>%
  group_by(GeneSymbol, CELL_ID, TREATMENT, CONCENTRATION, OX, INFLAM) %>%  # summarise the multiple probe ids
  summarise(log2FoldChange = mean(log2FoldChange, na.rm = TRUE) ) %>% 
  unite(col_id, CELL_ID, TREATMENT, CONCENTRATION) %>%
  spread(key = col_id, value = log2FoldChange)
}



