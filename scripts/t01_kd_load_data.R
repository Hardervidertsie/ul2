


load_data <- function(rela) {

  file_paths <- dir("../data_kd", full.names = TRUE, pattern = ".txt$")
  if(!rela){
    file_paths <- file_paths[!grepl("RelA", file_paths)]
  }
raw_data <- lapply(file_paths, read.delim, header = TRUE)
raw_data <- lapply(raw_data, function(x) {
  x[, !grepl("X[.]{0,1}[0-9]{0,12}", colnames(x))]
  } )

lapply(raw_data, head)

do.call('rbind', raw_data)
}



