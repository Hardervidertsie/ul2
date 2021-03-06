makeDoseLevels_cytokine <- function(dataIn){
  dataIn$treatcytokine <- paste(dataIn$treatment, dataIn$cytokine)
  dataIn <- dataIn[ order(  dataIn[ , "treatcytokine"], dataIn[ , "dose_uM"]),   ]
  
  counts.d <- by(data = dataIn, INDICES = dataIn[, "treatcytokine"], function(x) d.levels = unique(x[, "dose_uM"]))
  if(sum(lapply(counts.d, length ) >1  ) > 0 ) {  # if any compound has more than 1 dose level:
    counts.d.l <- sapply(counts.d, as.list)
    counts.d.l <- melt(counts.d.l,  length)
    old.nrow <- nrow(dataIn)
    dataIn <- merge(dataIn, counts.d.l, by.x = c( "treatcytokine", "dose_uM"), by.y = c( "L1", "value"), sort = FALSE)
    if(nrow(dataIn) != old.nrow){
      stop("setting dose levels for density plots failed")
    }
    dataIn$L2 <- factor(dataIn$L2)
    colnames(dataIn)[ colnames(dataIn) == "L2"] <- "doseLevel"
  }
  return(dataIn)
}
