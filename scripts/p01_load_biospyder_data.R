

load_data <- function(){ 
file_paths <- dir("../data_BioSpyder", full.names = TRUE)
raw_data <- lapply(file_paths, read.delim)
do.call('rbind', raw_data)
}

