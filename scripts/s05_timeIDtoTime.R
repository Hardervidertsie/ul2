
timeIDtoTimeFun <- function(dataIn, timeBetweenFrames,exposureDelay, plateID)
  
{
  
  timeBetweenFrames <- round(as.numeric(strftime(strptime(timeBetweenFrames, format = "%H:%M:%S"), "%H")) + 
                             1/60 * as.numeric(strftime(strptime(timeBetweenFrames, 
                                                                 format = "%H:%M:%S"), "%M")) +
                             1/3600 * as.numeric(strftime(strptime(timeBetweenFrames, 
                                                                   format = "%H:%M:%S"), "%S"))
                           , digit =3 )

exposureDelay <- round(as.integer(strftime(strptime(exposureDelay, format = "%H:%M"), "%H")) + 
                         1/60 * as.integer(strftime(strptime(exposureDelay, 
                                                             format = "%H:%M"), "%M")), digit =1 )

dataIn[dataIn$plateID %in% plateID, ]$timeID <- dataIn[dataIn$plateID %in% plateID, ]$timeID*timeBetweenFrames+exposureDelay - timeBetweenFrames
dataIn
}