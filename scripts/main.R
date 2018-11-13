rm(list=ls())
source("s00_lockAndload.R")
source("s02_run_format_functions.R")
source("s03_exploratory_plot.R")
source("s04_minmaxNorm.R")
source("s05_timeIDtoTime.R")
source("s06_makeDoseLevels.R")
options(stringsAsFactors = FALSE)

require(RColorBrewer )
require(ggplot2)
#install.packages("ggplot2")
# als nieuwe r sessie om nieuwe packages te installeren doe:
#trace(utils:::unpackPkgZip, edit=TRUE)
#regel 140 pas sys.sleep aan naar 2 seconden


#raw_data <- load_data(rootdir = "J:/Workgroups/FWN/LACDR/TOX/data steven wink/Unilever2/imaging_data_files",
#                      debug = FALSE)

#save(raw_data, file = "../tmp/raw_data.Rdata")

#rm(list=ls())
load("../tmp/raw_data.Rdata")

sort(unique(combined_data$treatment))

combined_data <- run_formats(raw_data)
head(combined_data)
dim(combined_data)

#write.table(combined_data, file = "../resultaten/combined_data2.txt", sep ="\t", col.names = NA)


#combined_data <- read.delim(file = "../resultaten/combined_data2.txt", sep ="\t", row.names = 1)

head(combined_data)


# data opschonen
head(combined_data)

unique(combined_data$variable)
#Integrated.Intensity.Cytoplasm -> Integrated.Intensity

combined_data$variable <- gsub("Integrated.Intensity.Cytoplasm", "Integrated.Intensity", combined_data$variable)
combined_data$variable <- gsub("Integrated.Intensity.Nuclei", "Integrated.Intensity", combined_data$variable)
combined_data$variable <- gsub("Integrated.Intenisty", "Integrated.Intensity", combined_data$variable)


combined_data$variable <- gsub("Count.Positives.mean...3x.sd", "countPositives_m3sd", combined_data$variable)
combined_data$variable <- gsub("Count.Positives.mean.3x.sd", "countPositives_m3sd", combined_data$variable)
combined_data$variable <- gsub("Count.Positives.mean...3x.sd", "countPositives_m3sd", combined_data$variable)



####Count.Positives.mean...3x.sd HMOX1
####Count.Positives.mean...3x.sd HSPA1B
####Count.Positives.mean...3x.sd Nrf2
####Count.Positives.mean...3x.sd SRXN1
####Count.Positives.mean.3x.sd A20
####Count.Positives.mean.3x.sd CHOP
####countPositives_m3sd BIP



#aanpassen variable naam
combined_data$variable <- paste(combined_data$variable, combined_data$cell_line)
unique(combined_data$variable)

#om data te verwijderen
#combined_data <- combined_data[ !(combined_data$variable == "first set SRXN1 rep1") , ]

#combined_data$treatment <- paste(combined_data$treatment, combined_data$sheet_name)
combined_data$plateID <- paste(combined_data$sheet_name, combined_data$cell_line, combined_data$replID)

head(combined_data)
unique(combined_data$sheet_name)
#om data te verwijderen
unique(combined_data$variable)
#combined_data <- combined_data[ !(combined_data$variable == "first set SRXN1 rep1") , ]
#combined_data <- combined_data[ !(combined_data$variable == "first set SRXN1 rep2") , ]
#combined_data <- combined_data[ !(combined_data$variable == "first set SRXN1 rep3") , ]
#combined_data <- combined_data[ !(combined_data$variable == "first set new concentrations SRXN1 rep3") , ]


sort(unique(combined_data$plateID))

####combined_data$variable <- gsub("Integrated.Intenisty", "Integrated.Intensity", combined_data$variable)


#combined_data <- combined_data[ !(combined_data$variable == "Integrated.Intensity" & 
                           #  combined_data$replID == "rep1") , ]

dim(combined_data)


# tijd aanpassen

combined_data$timeID <- as.numeric(combined_data$timeID)

  unique(combined_data[combined_data$plateID == "IL1 RelA rep3", c("timeID")])

  #first select data
  
  sort(unique(combined_data$variable))
  
  
  # Quercetin eruit
  
  combined_data <- combined_data[!combined_data$treatment == "Quercetin", ]
  
  unique(combined_data$treatment)
  combined_data$treatment <- gsub("Etacrynic acid" , "Ethacrynic acid", combined_data$treatment) 
  
  
  
  # selectie metingen
  #Ik wil in GraphPad time response figuren maken van de "Integrated.Intensity" voor de volgende reporters , ,  , , , , , , 
  #Voor de inflamation pathway reporters (A20, ICAM1 en ikbalpha)  dit graag van zowel de TNF exposure als de IL1 exposure en voor de No TNF or IL1 <name reporter> rep1 (hiervan is maar 1 replicate)
  
 
  combined_data$value <- as.numeric(combined_data$value)
  combined_data$dose_uM   <- as.numeric(combined_data$dose_uM)
  
  checkFun <- function(datain, showtreat = NULL, doselevel = FALSE) {
  dataout <- unique(datain[ , c("treatment", "dose_uM")] )
  dataout <- dataout[order(dataout$treatment, dataout$dose_uM), c("treatment", "dose_uM")]
  
  if(!doselevel) {
  if(!is.null(showtreat)){
    dataout[dataout$treatment %in% showtreat, ]
  }else{
    dataout
  }
  } else{
    dataout <- unique(datain[ , c("treatment", "dose_uM", "doseLevel")] )
    dataout <- dataout[order(dataout$treatment, dataout$dose_uM), c("treatment", "dose_uM", "doseLevel")]
    if(!is.null(showtreat)){
      dataout[dataout$treatment %in% showtreat, ]
    }else{
      dataout
    }
  }
  
  }
  
  checkFunCell <- function(datain, showtreat = NULL, cellLine = NULL) {
    dataout <- unique(datain[ , c("cell_line", "treatment", "dose_uM")] )
    dataout <- dataout[order(dataout$cell_line, dataout$treatment, dataout$dose_uM), c("cell_line", "treatment", "dose_uM")]
    dataout[ dataout$cell_line %in% cellLine & dataout$treatment %in% showtreat, ]
  }
  
  
  checkFun(datain = combined_data)
  checkFunCell(datain = combined_data, showtreat = "TBHQ", cellLine = "Nrf2")
  
  unique(combined_data[combined_data$cell_line == "Nrf2" & combined_data$treatment == "TBHQ",
                       c("cell_line", "treatment","variable", "dose_uM")])

  unique(combined_data[combined_data$cell_line == "RelA" , "variable"])
  combined_data$variable <- gsub("ratio.shNuclei.CytoN..Int.Int. RelA" , "ratio_shNuclei_CytoN", combined_data$variable)
  combined_data$variable <- gsub("ratio.shNuclei.CytoN RelA" , "ratio_shNuclei_CytoN", combined_data$variable)
  unique(combined_data$variable)
  sel_feats <- c("Integrated.Intensity A20", "Integrated.Intensity BIP", "Integrated.Intensity CHOP", 
                 "Integrated.Intensity HMOX1", "Integrated.Intensity HSPA1B", "Integrated.Intensity ICAM1",
                 "Integrated.Intensity IkBalpha", "Integrated.Intensity Nrf2", "Integrated.Intensity SRXN1",
                 "ratio_shNuclei_CytoN")
  
  sel_data <- combined_data[ combined_data$variable %in% sel_feats, ]
  head(sel_data)
  # check concentration levels for RelA
  unique(sel_data[sel_data$cell_line == "RelA", c("treatment", "dose_uM")])
  
  #write.table(unique(combined_data[combined_data$plateID == "IL1 RelA rep3", "timeID"]), file ="time6minEach.txt", sep = "\t")
  ##==##
  # create dose levels
  sel_data$dose_uM <- as.numeric(sel_data$dose_uM)
  head(sel_data)
  
  
#write.table(unique(sel_data[ , c("treatment", "dose_uM", "doseLevel")] ), file = "concentratie_selecties.txt", sep = "\t", col.names =NA)
  ##
  
  getwd()

  metadata <- read.table(file = '../metadata/concentratie_selecties_withcelldeath.txt', sep = "\t", header = TRUE, row.names = 1)   
   head(metadata)
   metadata$X.2 <- NULL
   metadata[metadata$treatment == "CDDO-Me", ]
   
   sel_data <- merge(sel_data, metadata, by = c("treatment", "dose_uM"), all.x = TRUE)
   unique(sel_data[sel_data$treatment == "CDDO-Me" ,c("treatment", "dose_uM")])
   
   #write.table(unique(subset(sel_data, X.1 == "afrondingsfout")[, c("treatment","dose_uM","doseLevel","X.1")]),
    #           file = "../metadata/afrondingsfouten.txt", sep ="\t")
   
   # merge based on dose levels. replace dose_uM. recalculate doseLevels
   fixed_dose <- read.table(file = "../metadata/afrondingsfouten.txt", sep ="\t", header = TRUE, row.names = 1)
   head(fixed_dose)
   colnames(fixed_dose)[2] <-"dose_uM"
   colnames(fixed_dose)[4] <- "soortFout"
   
   
   sel_data <- merge(sel_data, fixed_dose, by = c("treatment","dose_uM"), all.x = TRUE )
   
   sel_data$treatment <- gsub("Tert-Butylhydroquinone" , "TBHQ", sel_data$treatment) 
 head(sel_data)
 unique(sel_data[sel_data$treatment == "CDDO-Me" ,c("treatment", "dose_uM", "sheet_name", "X.1")])
 
 head(sel_data)
   #to takeoff_figure_acetaminophen.R
 
   
   sel_data <- sel_data[!sel_data$X.1 == "weg", ]
   
   
   sel_data<- sel_data[!sel_data$sheet_name %in% c("first set", "first set new concentrations"), ]
   
   
   yestmp <- unique(sel_data[, c("treatment", "include")])
   yestmp <- yestmp[yestmp$include %in% "yes",]
   yestmp # these compounds should remain.
   
   sel_data <- sel_data[ sel_data$treatment %in% yestmp$treatment, ]
   
   # check concentration levels for RelA
   unique(sel_data[sel_data$cell_line == "RelA", c("treatment", "dose_uM")])
   unique(sel_data[, c("treatment", "dose_uM")])
   # now manually fix round-off errors:
   tmp <- unique(sel_data[, c("treatment", "dose_uM",  "doseLevel.x",
                              "doseLevel.y", "X.1", "New.doselevel", "include")])
   
  
  
    fix_roundoff_error_fun <- function(x) {
     sel_data <- x
   sel_data[ sel_data$treatment == "Andrographolide" & sel_data$dose_uM == 0.320, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Andrographolide" & sel_data$dose_uM == 0.320, "dose_uM"] <- 0.316
   sel_data[ sel_data$treatment == "Andrographolide" & sel_data$dose_uM == 0.680, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Andrographolide" & sel_data$dose_uM == 0.680, "dose_uM"] <- 0.681 
   sel_data[ sel_data$treatment == "Andrographolide" & sel_data$dose_uM == 1.470, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Andrographolide" & sel_data$dose_uM == 1.470, "dose_uM"] <- 1.467 
   sel_data[ sel_data$treatment == "Andrographolide" & sel_data$dose_uM == 14.670, "include"] <- "yes" 
   sel_data[ sel_data$treatment == "Andrographolide" & sel_data$dose_uM == 14.670, "dose_uM"] <- 14.667
   sel_data[ sel_data$treatment == "Andrographolide" & sel_data$dose_uM == 6.810, "include"] <- "yes" 
   sel_data[ sel_data$treatment == "Andrographolide" & sel_data$dose_uM == 6.810, "dose_uM"] <- 6.808
   sel_data[ sel_data$treatment == "Andrographolide" & sel_data$dose_uM == 46.20, "include"] <- "yes" 
   sel_data[ sel_data$treatment == "Andrographolide" & sel_data$dose_uM == 46.20, "dose_uM"] <-   46.416
   sel_data[ sel_data$treatment == "Andrographolide" & sel_data$dose_uM == 146.67, "include"] <- "yes" 
   sel_data[ sel_data$treatment == "Andrographolide" & sel_data$dose_uM == 146.67, "dose_uM"] <-   146.674
    
   
   
   
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 0.020, "include"] <- "yes"
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 0.020, "dose_uM"] <- 0.018
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 0.030, "include"] <- "yes"
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 0.030, "dose_uM"] <- 0.0316
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 0.032, "include"] <- "yes" # stond er nog niet
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 0.032, "dose_uM"] <- 0.0316
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 0.060 , "include"] <- "yes"
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 0.060 , "dose_uM"] <- 0.056
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 0.180 , "include"] <- "yes"
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 0.180 , "dose_uM"] <- 0.178 
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 1.780 , "include"] <- "yes"
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 1.780 , "dose_uM"] <- 1.777 
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 0.320 , "include"] <- "yes"
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 0.320 , "dose_uM"] <- 0.316
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 0.560 , "include"] <- "yes"
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 0.560 , "dose_uM"] <- 0.562
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 1.000 , "include"] <- "yes"
   sel_data[ sel_data$treatment == "CDDO-Me" & sel_data$dose_uM == 1.000 , "dose_uM"] <- 0.999
   
   
   sel_data[ sel_data$treatment == "DEM" & sel_data$dose_uM == 2.150, "include"] <- "yes"
   sel_data[ sel_data$treatment == "DEM" & sel_data$dose_uM == 2.150, "dose_uM"] <- 2.154
   sel_data[ sel_data$treatment == "DEM" & sel_data$dose_uM == 4.640, "include"] <- "yes"
   sel_data[ sel_data$treatment == "DEM" & sel_data$dose_uM == 4.640, "dose_uM"] <- 4.642
   sel_data[ sel_data$treatment == "DEM" & sel_data$dose_uM == 21.540, "include"] <- "yes"
   sel_data[ sel_data$treatment == "DEM" & sel_data$dose_uM == 21.540, "dose_uM"] <- 21.544
   sel_data[ sel_data$treatment == "DEM" & sel_data$dose_uM == 46.420, "include"] <- "yes"
   sel_data[ sel_data$treatment == "DEM" & sel_data$dose_uM == 46.420, "dose_uM"] <- 46.416
   sel_data[ sel_data$treatment == "DEM" & sel_data$dose_uM == 215.440, "include"] <- "yes"
   sel_data[ sel_data$treatment == "DEM" & sel_data$dose_uM == 215.440, "dose_uM"] <- 215.443
   sel_data[ sel_data$treatment == "DEM" & sel_data$dose_uM == 464.16, "include"] <- "yes"
   sel_data[ sel_data$treatment == "DEM" & sel_data$dose_uM == 464.16, "dose_uM"] <- 464.159
   
   
   sel_data[ sel_data$treatment == "Diclofenac" & sel_data$dose_uM == 2.150, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Diclofenac" & sel_data$dose_uM == 2.150, "dose_uM"] <- 2.154
   sel_data[ sel_data$treatment == "Diclofenac" & sel_data$dose_uM == 4.640, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Diclofenac" & sel_data$dose_uM == 4.640, "dose_uM"] <- 4.642
   sel_data[ sel_data$treatment == "Diclofenac" & sel_data$dose_uM == 21.540, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Diclofenac" & sel_data$dose_uM == 21.540, "dose_uM"] <- 21.544
   sel_data[ sel_data$treatment == "Diclofenac" & sel_data$dose_uM == 46.420, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Diclofenac" & sel_data$dose_uM == 46.420, "dose_uM"] <- 46.416
   sel_data[ sel_data$treatment == "Diclofenac" & sel_data$dose_uM == 215.440, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Diclofenac" & sel_data$dose_uM == 215.440, "dose_uM"] <- 215.443
   sel_data[ sel_data$treatment == "Diclofenac" & sel_data$dose_uM == 464.160, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Diclofenac" & sel_data$dose_uM == 464.160, "dose_uM"] <- 464.159
   
   sel_data[ sel_data$treatment == "Ethacrynic acid" & sel_data$dose_uM == 1.780, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Ethacrynic acid" & sel_data$dose_uM == 1.780, "dose_uM"] <- 1.777
   sel_data[ sel_data$treatment == "Ethacrynic acid" & sel_data$dose_uM == 5.620, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Ethacrynic acid" & sel_data$dose_uM == 5.620, "dose_uM"] <- 5.619
   sel_data[ sel_data$treatment == "Ethacrynic acid" & sel_data$dose_uM == 9.990, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Ethacrynic acid" & sel_data$dose_uM == 9.990, "dose_uM"] <- 9.993
   sel_data[ sel_data$treatment == "Ethacrynic acid" & sel_data$dose_uM == 56.190, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Ethacrynic acid" & sel_data$dose_uM == 56.190, "dose_uM"] <- 56.194
   sel_data[ sel_data$treatment == "Ethacrynic acid" & sel_data$dose_uM == 99.930, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Ethacrynic acid" & sel_data$dose_uM == 99.930, "dose_uM"] <- 99.928
   
   sel_data[ sel_data$treatment == "Sulforaphane" & sel_data$dose_uM == 0.220, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Sulforaphane" & sel_data$dose_uM == 0.220, "dose_uM"] <- 0.215
   sel_data[ sel_data$treatment == "Sulforaphane" & sel_data$dose_uM == 0.460, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Sulforaphane" & sel_data$dose_uM == 0.460, "dose_uM"] <- 0.464
   sel_data[ sel_data$treatment == "Sulforaphane" & sel_data$dose_uM == 2.150, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Sulforaphane" & sel_data$dose_uM == 2.150, "dose_uM"] <- 2.154
   sel_data[ sel_data$treatment == "Sulforaphane" & sel_data$dose_uM == 4.640, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Sulforaphane" & sel_data$dose_uM == 4.640, "dose_uM"] <- 4.642
   sel_data[ sel_data$treatment == "Sulforaphane" & sel_data$dose_uM == 21.540, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Sulforaphane" & sel_data$dose_uM == 21.540, "dose_uM"] <- 21.544
   sel_data[ sel_data$treatment == "Sulforaphane" & sel_data$dose_uM == 46.42, "include"] <- "yes"
   sel_data[ sel_data$treatment == "Sulforaphane" & sel_data$dose_uM == 46.42, "dose_uM"] <- 46.416
   
   
   
   sel_data[ sel_data$treatment == "TBHQ" & sel_data$dose_uM == 0.320, "include"] <- "yes"
   sel_data[ sel_data$treatment == "TBHQ" & sel_data$dose_uM == 0.320, "dose_uM"] <- 0.316
   sel_data[ sel_data$treatment == "TBHQ" & sel_data$dose_uM == 0.680, "include"] <- "yes"
   sel_data[ sel_data$treatment == "TBHQ" & sel_data$dose_uM == 0.680, "dose_uM"] <- 0.681
   sel_data[ sel_data$treatment == "TBHQ" & sel_data$dose_uM == 1.470, "include"] <- "yes"
   sel_data[ sel_data$treatment == "TBHQ" & sel_data$dose_uM == 1.470, "dose_uM"] <- 1.467
   sel_data[ sel_data$treatment == "TBHQ" & sel_data$dose_uM == 6.810 , "include"] <- "yes"
   sel_data[ sel_data$treatment == "TBHQ" & sel_data$dose_uM == 6.810 , "dose_uM"] <- 6.808
   sel_data[ sel_data$treatment == "TBHQ" & sel_data$dose_uM == 14.670 , "include"] <- "yes"
   sel_data[ sel_data$treatment == "TBHQ" & sel_data$dose_uM == 14.670 , "dose_uM"] <- 14.667
   sel_data[ sel_data$treatment == "TBHQ" & sel_data$dose_uM == 146.670 , "include"] <- "yes"
   sel_data[ sel_data$treatment == "TBHQ" & sel_data$dose_uM == 146.670 , "dose_uM"] <- 146.674
   
   sel_data[ sel_data$treatment == "Troglitazone" & sel_data$dose_uM == 0.320 , "include"] <- "yes"
   sel_data[ sel_data$treatment == "Troglitazone" & sel_data$dose_uM == 0.320 , "dose_uM"] <- 0.316
   sel_data[ sel_data$treatment == "Troglitazone" & sel_data$dose_uM == 0.680 , "include"] <- "yes"
   sel_data[ sel_data$treatment == "Troglitazone" & sel_data$dose_uM == 0.680 , "dose_uM"] <- 0.681
   sel_data[ sel_data$treatment == "Troglitazone" & sel_data$dose_uM == 6.810 , "include"] <- "yes"
   sel_data[ sel_data$treatment == "Troglitazone" & sel_data$dose_uM == 6.810 , "dose_uM"] <- 6.808
   sel_data[ sel_data$treatment == "Troglitazone" & sel_data$dose_uM == 14.670 , "include"] <- "yes"
   sel_data[ sel_data$treatment == "Troglitazone" & sel_data$dose_uM == 14.670 , "dose_uM"] <- 14.667
   sel_data[ sel_data$treatment == "Troglitazone" & sel_data$dose_uM == 146.670 , "include"] <- "yes"
   sel_data[ sel_data$treatment == "Troglitazone" & sel_data$dose_uM == 146.670 , "dose_uM"] <- 146.674
   sel_data
   }
   
   #test_data <- fix_roundoff_error_fun(sel_data)
   #test_data <- test_data[test_data$include == "yes", ]
   
   
   #unique(sel_data[sel_data$cell_line == "HMOX1" & sel_data$treatment == "Andrographolide", c("treatment", "dose_uM",  "cell_line")   ]) 
   #unique(test_data[test_data$cell_line == "HMOX1" & test_data$treatment == "Andrographolide", c("treatment", "dose_uM",  "cell_line")   ]) 
   
   
   sel_data <- fix_roundoff_error_fun(sel_data)
   
   sel_data <- sel_data[sel_data$include == "yes", ]
   
   
    sel_data$doseLevel <- NULL
   
   
   sel_data <- makeDoseFacFun(sel_data)
   
   #unique(sel_data[sel_data$cell_line == "RelA" , c("treatment", "dose_uM", "doseLevel")] )
   unique(sel_data[ , c("treatment", "dose_uM", "doseLevel")] )
   # plaat normalizatie
  sel_data$value <- as.numeric(sel_data$value)
   
  head(sel_data[sel_data$cell_line == "RelA", ])
  #tmp <- subset(sel_data, cell_line == "RelA")
  #sel_data <- subset(sel_data, cell_line != "RelA")
  
  head(tmp)
  
  ggplot(data = subset(tmp, cell_line == "RelA"), aes(x = timeID, y = value)) + 
    geom_point(aes(color = sheet_name, shape = replID))+
    facet_grid(treatment ~ doseLevel) + theme_classic()
  
  
  # relA should not be plate normalized? (plates with/ without IL1? ) Is this also the case for ICAM1/ A20 etc?
  # or do all plates have same positive controls?
  
  sel_data_list <- split(sel_data, f = sel_data$variable)
 
    sel_data_list <- lapply(sel_data_list, normFun)
    sel_data_fixed <- do.call('rbind', sel_data_list) 
  tail(sel_data_fixed)
  
 # sel_data_fixed <- rbind(sel_data_fixed, tmp)
  unique(sel_data_fixed[, "variable"])
    dim(sel_data_fixed)
    sel_data_fixed <- aggregate(value ~ plateID + cell_line + treatment + dose_uM + timeID + sheet_name + variable + replID + doseLevel,
                                    data = sel_data_fixed, FUN = function(x) mean(x, na.rm = TRUE))
  
    
    
    #head(sel_data_fixed)
 # unique(sel_data[, c("treatment", "dose_uM", "doseLevel")])
  # format data for graphpad
  
  sel_data_fixed$treatment <- paste(sel_data_fixed$treatment, sel_data_fixed$sheet_name, sep ="_")
 # unique(sel_data_fixed$treatment)
#  sel_data_fixed[grepl("CDDO" ,sel_data_fixed$treatment) ,]
  
  
all.feats <- unique(sel_data_fixed$variable)

  

 for( i in seq_along(all.feats)){
  sel_feat <- sel_data_fixed[sel_data_fixed$variable == all.feats[i], ]
  all.treats <- unique(sel_feat$treatment)
  for( j in seq_along(all.treats)){
  sel_feat_treat <- sel_feat[sel_feat$treatment == all.treats[j], ]
  
  head(sel_feat_treat)

    sel_feat_casted <- dcast(  timeID ~ doseLevel + replID + sheet_name, value.var = "value"  , data = sel_feat_treat)
filename <- paste( unique(sel_feat_treat$cell_line ), unique(sel_feat_treat$treatment )
                   ,unique(sel_feat_treat$variable ), sep = "__")

 write.table(sel_feat_casted, file = paste('../resultaten/graphpad format/withdeadconditions/', filename, ".txt", sep =""), sep = "\t", row.names = FALSE)
  
 
  }
 }

unique(sel_data[ sel_data$cell_line == "Nrf2" & sel_data$treatment == "TBHQ" , "doseLevel"])

# Nrf2__TBHQ_Data__mmn_Integrated.Intensity Nrf2 4, 7 & 8.
# bij HMOX ook.
## == processing of HepG2 data:


unique( combined_data$cell_line)
unique(combined_data[ combined_data$cell_line == "WT cell death", "variable"])

sel_feats_wt <- c("ratio AreaPI/AreaNuclei WT cell death", "ratio AreaAnV/AreaNuclei WT cell death")

sel_data_wt <- combined_data[ combined_data$variable %in% sel_feats_wt, ]



head(sel_data_wt)
table(unique(sel_data_wt[, c("plateID", "timeID")]))

sel_data_wt <- sel_data_wt[ sel_data_wt$timeID == "22",]

# create dose levels
sel_data_wt$dose_uM <- as.numeric(sel_data_wt$dose_uM)

sel_data_wt <- makeDoseFacFun(sel_data_wt)

unique(sel_data_wt[ , c("treatment", "dose_uM", "doseLevel")] )

# plaat normalizatie
sel_data_wt$value <- as.numeric(sel_data_wt$value)

sel_data_wt_fixed <- sel_data_wt
dim(sel_data_wt_fixed)

sel_data_wt_fixed <- aggregate(value ~ plateID + cell_line + treatment + dose_uM + timeID + sheet_name + variable + replID + doseLevel,
                            data = sel_data_wt_fixed, FUN = function(x) mean(x, na.rm = TRUE))



# format data for graphpad


# first plot data, then format for graphpad.
ggplot(data = sel_data_wt_fixed, aes(x = doseLevel, y= value, color = replID )) + facet_grid(variable + sheet_name ~ treatment) +
  geom_point() + geom_smooth(aes(group = replID), method = 'loess') + theme_bw()

#sel_data_wt_fixed$treatment <- paste(sel_data_wt_fixed$treatment, sel_data_wt_fixed$sheet_name, sep ="_")

all.feats <- unique(sel_data_wt_fixed$variable)

i=1;j=1

for( i in seq_along(all.feats)){
  sel_feat <- sel_data_wt_fixed[sel_data_wt_fixed$variable == all.feats[i], ]
  all.treats <- unique(sel_feat$treatment)
  for( j in seq_along(all.treats)){
    sel_feat_treat <- sel_feat[sel_feat$treatment == all.treats[j], ]
    
    
    
    sel_feat_casted <- dcast(  dose_uM ~ replID + sheet_name, value.var = "value"  , data = sel_feat_treat)
    filename <- paste( unique(sel_feat_treat$treatment )
                       ,unique(sel_feat_treat$variable ), sep = "__")
    filename <- gsub("/", ".div.", filename)
    write.table(sel_feat_casted, file = paste('../resultaten/graphpad format/HepG2WT/', filename, ".txt", sep =""), sep = "\t", row.names = FALSE)
    
    
  }
}

# dose levels voor plotten
  #sel_data$treatment <- as.character(sel_data$treatment)
  sel_data <- makeDoseFacFun(dataIn = sel_data)
head(sel_data)
  
  unique(sel_data[, c("treatment", "replID", "dose_uM", "doseLevel")])


# data normalisatie
  # swap to s07_heatmap.R
unique(sel_data$variable)

norm_data <- normFun(dataIn = sel_data, variable = "") # min max: kan allen als op plaat basis neg en pos controles hebt


head(sel_data)

exploratory_plot(norm_data)
exploratory_plot_dose(norm_data)

exploratory_plot(sel_data)
exploratory_plot_dose(sel_data)
