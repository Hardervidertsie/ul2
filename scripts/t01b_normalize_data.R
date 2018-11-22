normalize_data <- function(input_data) {
  
  # min max normalize Integrated.Intensity column per cell line repID block
  
  minmaxFun <- function(x) {
    (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
  }
  
  
  # create dose levels
  
  
  dose_level_table <- input_data %>% group_by(cell_line,treatment, dose_uM) %>%
    summarise() %>%
    mutate(dose_level = 1: n())
  
  dose_level_table <- dose_level_table %>% mutate(dose_level = dose_level + 1) %>% 
    mutate(dose_level = if_else(treatment == "DMSO", 1, dose_level))
  
  
  input_data <- left_join(input_data, dose_level_table, by = c("cell_line", "treatment", "dose_uM"))
  
  # normalize
  
  input_data %>% group_by(cell_line, repID) %>%
    mutate(Integrated.Intensity = minmaxFun(Integrated.Intensity))
  
  
  
  
}