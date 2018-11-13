


format_data3 <- function(inputdata) {
# messy and manual formatting of excel data
  # format individual data.frames
  # melt individual data.frames
  # rbind them
  require(reshape2)

file1sheet1 <- inputdata
cnf1s1 <- colnames(file1sheet1) 
indf1s1 <- which(grepl("rep1", as.character(file1sheet1[1,])))
# also try if rep2 (file 7 started at rep2 and went to rep4)/ another only had rep3
if(length(indf1s1) == 0){
  indf1s1 <- which(grepl("rep2", as.character(file1sheet1[1,])))
}
if(length(indf1s1) == 0){
  indf1s1 <- which(grepl("rep3", as.character(file1sheet1[1,])))
}




indcnf1s1 <- which(grepl("rep", as.character(file1sheet1[1,])))
if(length(indcnf1s1)!=0){
nRep1 <- length(unique(indcnf1s1)) %/% length(unique(indf1s1))

if( nRep1 != length(unique(indcnf1s1)) / length(unique(indf1s1)) ){
  stop("not all features identical replicate number?")
}
} else{
  nRep1 <- 1
}

ind2f1s1 <- rep(indf1s1, each = nRep1)

newColNames <-  paste(cnf1s1[ind2f1s1] ,  as.character(file1sheet1[1,])[indcnf1s1], sep = "__")

colnames(file1sheet1)[indcnf1s1] <- newColNames
colfrom2toheader <- which(grepl("(X.[0-9]{1,2})$", colnames(file1sheet1)) | colnames(file1sheet1) == "X")

colnames(file1sheet1)[colfrom2toheader] <- as.character(file1sheet1[1,])[colfrom2toheader]

ind <- which(sapply(file1sheet1, function(x)  sum(!is.na(x))) != 0)

file1sheet1 <- file1sheet1[, ind] 
if(!length(indf1s1) == 0 ) {
file1sheet1 <- file1sheet1[-1,  ]
}

file1sheet1_melt <- melt(data = file1sheet1, id.vars = c("treatment", "dose_uM", "timeID", "sheet_name", "cell_line")  )

if(!length(indf1s1) == 0 ) {
file1sheet1_melt$replID <- sapply(
  strsplit(as.character(file1sheet1_melt$variable), "__"), "[[", 2)

file1sheet1_melt$variable <- gsub("(__rep[0-9]{1})$", "", as.character(file1sheet1_melt$variable))
} else {
  file1sheet1_melt$replID <- "rep1"
}

file1sheet1_melt

  
}
