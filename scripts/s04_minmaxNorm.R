normFun <- function(dataIn, variable = NULL) {
 
  if(!is.null(variable)){
   dataIn <- dataIn[ dataIn$variable == variable, ]
  } else {
    variable <- unique(dataIn$variable)
  }
  minV <- aggregate(data = dataIn, value ~ plateID + cell_line, 
                    FUN = function(x) min(x, na.rm = TRUE))
  maxV <- aggregate(data = dataIn, value ~ plateID + cell_line, 
                    FUN = function(x) max(x, na.rm = TRUE))
  
  colnames(minV)[colnames(minV) == "value"] <- "minV"
  colnames(maxV)[colnames(maxV) == "value"] <- "maxV"
  
  dataIn <- merge(dataIn, minV, by = c("plateID", "cell_line"), all.x = TRUE)
  dataIn <- merge(dataIn, maxV, by = c("plateID", "cell_line"), all.x = TRUE)
  
  dataIn$norm_value <- (dataIn$value - dataIn$minV) / (dataIn$maxV - dataIn$minV)
  dataIn$minV <- NULL
  dataIn$maxV <- NULL
  dataIn$value <- NULL
  colnames(dataIn)[colnames(dataIn) == "norm_value"] <- "value"
  dataIn$variable <- paste0("mmn_", dataIn$variable)
  dataIn
}