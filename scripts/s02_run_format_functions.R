



run_formats <- function(raw_data) {

source("s01_format_data.R")

length(raw_data) # 11
sapply(raw_data, length) # 3 1 1 1 1 3 2 1 2 3 2


head(raw_data[[1]][[2]])
head(raw_data[[1]][[3]])

data1 <- format_data3(raw_data[[1]][[1]])
data2 <- format_data3(raw_data[[1]][[2]])
data3 <- format_data3(raw_data[[1]][[3]])
head(data3)

##==

head(raw_data[[2]][[1]])

data4 <- format_data3(raw_data[[2]][[1]])
head(data4)

##==

head(raw_data[[3]][[1]])

data5 <- format_data3(raw_data[[3]][[1]])
head(data5)

##==

head(raw_data[[4]][[1]])

data6 <- format_data3(raw_data[[4]][[1]])
head(data6)

##==

head(raw_data[[5]][[1]])

data7 <- format_data3(raw_data[[5]][[1]])
head(data7)

##==

# kolom naam in raw_data[[6]][[2]] miste: PI rep3. In excel file toegevoegd. lage waarden vergeleken met andere 2 reps.
head(raw_data[[6]][[2]])

data8 <- format_data3(raw_data[[6]][[1]])
data9 <- format_data3(raw_data[[6]][[2]])
data10 <- format_data3(raw_data[[6]][[3]])

head(data10)

##==
head(raw_data[[7]][[2]])
# added missing repl4 name column 10 of file 7 sheet 1
data11 <- format_data3(raw_data[[7]][[1]])
data12 <- format_data3(raw_data[[7]][[2]])


head(data12)


##==
head(raw_data[[8]][[1]])

data13 <- format_data3(raw_data[[8]][[1]])

head(data13)


##==
head(raw_data[[9]][[2]])

data14 <- format_data3(raw_data[[9]][[1]])
data15 <- format_data3(raw_data[[9]][[2]])

head(data15)


##==
head(raw_data[[10]][[2]])

data16 <- format_data3(raw_data[[10]][[1]])
data17 <- format_data3(raw_data[[10]][[2]])
data18 <- format_data3(raw_data[[10]][[3]])

head(data18)



##==
head(raw_data[[11]][[2]])

# WT data headers is different. Changing in script.

head(raw_data[[11]][[1]])

raw_data[[11]][[1]]$DMEM[-c(1,2)] <- paste(raw_data[[11]][[1]]$DMEM[-c(1,2)] , "DMEM", sep = "_")
raw_data[[11]][[1]]$TNF[-c(1,2)] <- paste(raw_data[[11]][[1]]$TNF[-c(1,2)] , "TNF", sep = "_")


colnames(raw_data[[11]][[1]]) <- raw_data[[11]][[1]][1, ]
raw_data[[11]][[1]] <- raw_data[[11]][[1]][ -1, ]

colnames(raw_data[[11]][[1]])[colnames(raw_data[[11]][[1]]) == ""] <- paste0("X", 
                                                                             c("" ,paste0(
                                                                               "." ,seq(
                                                                                 1, (
                                                                                   sum(
                                                                                   colnames(
                                                                                     raw_data[[11]][[1]]
                                                                                     ) == ""
                                                                                   )-1
                                                                                )  
                                                                              )
                                                                            )
                                                                          )
                                                                        )

colnames(raw_data[[11]][[1]])[ 33:34] <- c("sheet_name", "cell_line")

head(raw_data[[11]][[1]][, c(1:16, 33,34)])
head(raw_data[[11]][[1]][, 17:34])

data19 <- format_data3(raw_data[[11]][[1]][, c(1:16, 33,34)])

data20 <- format_data3(raw_data[[11]][[1]][, 17:34])


##===##

head(raw_data[[11]][[2]])

raw_data[[11]][[2]]$DMEM[-c(1,2)] <- paste(raw_data[[11]][[2]]$IL1[-c(1,2)] , "DMEM", sep = "_")
raw_data[[11]][[2]]$IL1[-c(1,2)] <- paste(raw_data[[11]][[2]]$IL1[-c(1,2)] , "IL1", sep = "_")
bkuptmp <- raw_data[[11]][[2]]
head(bkuptmp)
colnames(raw_data[[11]][[2]]) <- raw_data[[11]][[2]][1, ]
raw_data[[11]][[2]] <- raw_data[[11]][[2]][ -1, ]

colnames(raw_data[[11]][[2]])[colnames(raw_data[[11]][[2]]) == ""] <- paste0("X", 
                                                                             c("" ,paste0(
                                                                               "." ,seq(
                                                                                 1, (
                                                                                   sum(
                                                                                     colnames(
                                                                                       raw_data[[11]][[2]]
                                                                                     ) == ""
                                                                                   )-1
                                                                                 )  
                                                                               )
                                                                             )
                                                                           )
                                                                        ) 


head(raw_data[[11]][[2]])

colnames(raw_data[[11]][[2]])[ 27:28] <- c("sheet_name", "cell_line")

head(raw_data[[11]][[2]][, c(1:14, 27,28)])

head(raw_data[[11]][[2]][, 15:28])

data21 <- format_data3(raw_data[[11]][[2]][, c(1:14, 27,28)])

data22 <- format_data3(raw_data[[11]][[2]][, 15:28])


combined_data <- rbind(data1,data2,data3, data4,data5,data6,data7,data8,data9,data10,data11,data12,
                       data13,data14,data15,data16,data17,data18,data19,data20,data21,data22)


combined_data

}
#file1sheet1 <-raw_data[[11]][[1]]

