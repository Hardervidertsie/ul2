




load_data <- function(rootdir, debug = FALSE){
  options(stringsAsFactors = FALSE)
    #trace(utils:::unpackPkgZip, edit=TRUE)
  require(gdata)
  require(stringr)
  if(debug ){
    rootdir <- "J:/Workgroups/FWN/LACDR/TOX/data steven wink/Unilever2/imaging_data_files"
      }
  # load data
  l1files <- list.files(rootdir, full.names = TRUE)
  l1files <- l1files[!grepl("~", l1files)]
  FResults <- lapply(l1files, function(filename) {
      
    
    cell_line <- str_extract(string = filename,
                pattern = "(/[A-Z a-z 0-9 ]{3,15}.xlsx)$")
    
    cell_line <- gsub("/", "", cell_line)
    cell_line <- gsub(".xlsx", "", cell_line)
    
    sheet_names <- sheetNames(filename) 
    sheet_names <- sheet_names[!sheet_names %in% "Comments"]
      
       lapply( sheet_names, function(sn) { 
        tmp <- read.xls(xls = filename, sheet = sn,
                perl = "D:/apps/perl/perl/bin/perl.exe")
        tmp$sheet_name <- sn
        tmp$cell_line <- cell_line
        tmp
      }
        
        )
      
     }
    )

    
   
    FResults
}