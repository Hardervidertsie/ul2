select_data <- function(data_in) {

data_in %>% 
 
    filter(siRNA %in% c("siNrf2", "siRelA", "siCtrl1", "siCtrl")) %>%
    filter(!(cell_line %in% c("Nrf2", "SRXN1") & control == "TNF")) 
    

}

