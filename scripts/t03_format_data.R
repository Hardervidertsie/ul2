

format_data <- function(sel_data) {

# create dose levels


dose_level_table <- sel_data %>% group_by(cell_line,treatment, dose_uM) %>%
  summarise() %>%
  mutate(dose_level = 1: n())

dose_level_table <- dose_level_table %>% mutate(dose_level = dose_level + 1) %>% 
  mutate(dose_level = if_else(treatment == "DMSO", 1, dose_level))


sel_data <- left_join(sel_data, dose_level_table, by = c("cell_line", "treatment", "dose_uM"))

sel_data <- sel_data %>% mutate(row_id = paste(cell_line, treatment, control, dose_level))



sel_data$row_id <- factor(sel_data$row_id, levels = 
                            c("Nrf2 DMSO DMEM 1", "Nrf2 DEM DMEM 2", 
                              "Nrf2 DEM DMEM 3",  "Nrf2 DEM DMEM 4",
                              "SRXN1 DMSO DMEM 1", "SRXN1 DEM DMEM 2",
                              "SRXN1 DEM DMEM 3", "SRXN1 DEM DMEM 4",
                              "A20 DMSO TNF 1", "A20 DEM TNF 2",
                              "A20 DEM TNF 3", "A20 DEM TNF 4", 
                              "ICAM1 DMSO TNF 1",   "ICAM1 DEM TNF 2",   
                              "ICAM1 DEM TNF 3",  "ICAM1 DEM TNF 4"))

sel_data$timeID <- factor(sel_data$timeID, levels = unique(sel_data$timeID))

sel_data$siRNA <- gsub("siCtrl$", "siCtrl1", sel_data$siRNA)

sel_data$siRNA <- factor(sel_data$siRNA, levels = c("siNrf2", "siRelA", "siCtrl1"))


sel_data %>%  select(row_id, siRNA, timeID, Integrated.Intensity) %>% 
  mutate(siRNA_time = paste(siRNA, timeID)) %>%
  select(-timeID,-siRNA) %>% group_by(row_id, siRNA_time) %>%
  summarise(GFP_signal = mean(Integrated.Intensity, na.rm = TRUE)) %>% 
  spread( key = siRNA_time, value = GFP_signal)

}
