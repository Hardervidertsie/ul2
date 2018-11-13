
# script that 
#load('../tmp/12_02_2018.RData')
ls()
tail(sel_data_fixed)

unique(sel_data_fixed$variable)

tmp <- unique(sel_data_fixed[, c("treatment", "dose_uM")])
tmp[order(tmp$treatment, tmp$dose_uM),]

tmp_wt <- unique(sel_data_wt_fixed[, c("treatment", "dose_uM")])
tmp_wt[order(tmp_wt$treatment, tmp_wt$dose_uM),]
unique(tmp_wt[, c("treatment", "dose_uM")])
#1 select same compounds for wt as reporter data
#2 fix wt concentration round off errors, 
#3 create a compound-concentration aggregate table to select correct data set from wt
#4 add wt to reporter data
head(sel_data_wt_fixed)
head(sel_data_fixed)
unique(sel_data_fixed$treatment)
unique(sel_data_fixed[sel_data_fixed$sheet_name %in% "TNF", "cell_line"])

cytokineReporters <- c("A20", "ICAM1", "IkBalpha", "RelA")

unique(sel_data_fixed[sel_data_fixed$treatment == "Sulforaphane_Data", "cell_line"])

# split treatment by underscores into treatment and cytokine
# add reporter level TNF, IL1 and DMEM to the cytokineReporters
tmp_reporter <- sel_data_fixed
tmp_wt <- sel_data_wt_fixed

tmp_reporter$cytokine <- sapply(sapply(tmp_reporter$treatment, strsplit, "_"), "[[", 2)
tmp_reporter$treatment <- sapply(sapply(tmp_reporter$treatment, strsplit, "_"), "[[", 1)
unique(tmp_reporter$treatment)
unique(tmp_reporter$cytokine)

tmp_wt$cytokine <- sapply(sapply(tmp_wt$treatment, strsplit, "_"), "[[", 2)
tmp_wt$treatment <- sapply(sapply(tmp_wt$treatment, strsplit, "_"), "[[", 1)
unique(tmp_wt$treatment)
unique(tmp_wt$cytokine)

tmp_reporter$cytokine[ !tmp_reporter$cytokine %in% c("TNF", "IL1")] <- "DMEM"
# now fix round off errors
# then select correct treatment & doses


tmp_wt$include <- "dummy"
tmp_reporter <- fix_roundoff_error_fun(tmp_reporter)



tmp_wt <- fix_roundoff_error_fun(tmp_wt)

chck_wt <-unique(tmp_wt[, c("treatment", "dose_uM")])
chck_wt[order(chck_wt$treatment, chck_wt$dose_uM), ]

chck_rp <-unique(tmp_reporter[, c("treatment", "dose_uM")])
chck_rp[order(chck_rp$treatment, chck_rp$dose_uM), ]

tmp_wt$mergecol <- paste0(tmp_wt$treatment, tmp_wt$dose_uM)
tmp_reporter$mergecol <- paste0(tmp_reporter$treatment, tmp_reporter$dose_uM)

sel_ind <- which(tmp_wt$mergecol %in% tmp_reporter$mergecol)

tmp_wt <- tmp_wt[sel_ind, ]
cbind(chck_rp[order(chck_rp$treatment, chck_rp$dose_uM), ],
      chck_wt[order(chck_wt$treatment, chck_wt$dose_uM), ])
colnames(tmp_reporter)
colnames(tmp_wt)
tmp_reporter$doseLevel <- NULL
tmp_wt$doseLevel <- NULL
source('s08_makeDoseLevels_cytokine.R')


tmp_reporter<-makeDoseLevels_cytokine(tmp_reporter)
tmp_wt<-makeDoseLevels_cytokine(tmp_wt)

unique(tmp_reporter[, c("treatment", "dose_uM", "doseLevel")])

head(tmp_reporter)
head(tmp_reporter)

heatmap_data <- rbind(tmp_reporter, tmp_wt)
head(heatmap_data)

unique(heatmap_data$treatcytokine)


# now create a square table, with rows as unique treatments and columns as features + doseLevel
# first select max per time course, then calculate average over replicates
# what about 

heatmap_data[ heatmap_data$variable == "ratio AreaAnV/AreaNuclei WT cell death", "cell_line"] <- "HepG2 WT"
heatmap_data[ heatmap_data$variable == "ratio AreaPI/AreaNuclei WT cell death", "cell_line"] <- "HepG2 WT"
heatmap_data[ heatmap_data$variable == "ratio AreaAnV/AreaNuclei WT cell death", "variable"] <- "Apoptosis WT"
heatmap_data[ heatmap_data$variable == "ratio AreaPI/AreaNuclei WT cell death", "variable"] <- "Necrosis WT"


heatmap_data$feature <- paste(heatmap_data$cell_line, heatmap_data$cytokine, sep = "_")
unique(heatmap_data$feature)

# fix DMSO (only 1 dose level)



# based on heatmap result: might need some control based normalization (e.g. ICAM1 treatment < control)
require(dplyr)
head(heatmap_data)
#HSPA1B eruit
unique(heatmap_data$cell_line)

heatmap_data <- heatmap_data[heatmap_data$cell_line != "HSPA1B", ]


heatmap_data$variable

# save data with time
all.feats <- unique(heatmap_data$variable)
# hier gebleven: 13-04-2018

for( i in seq_along(all.feats)){
  sel_feat <- heatmap_data[heatmap_data$variable == all.feats[i], ]
  all.treats <- unique(sel_feat$treatment)
  for( j in seq_along(all.treats)){
    sel_feat_treat <- sel_feat[sel_feat$treatment == all.treats[j], ]
    
    sel_feat_casted <- dcast(  timeID ~ doseLevel + replID + sheet_name + feature, value.var = "value"  , data = sel_feat_treat)
    filename <- gsub("/", "div", paste( unique(sel_feat_treat$cell_line ), unique(sel_feat_treat$treatment )
                       ,unique(sel_feat_treat$variable ), sep = "__"))
    
    write.table(sel_feat_casted, file = paste('../resultaten/graphpad format/deadcellsRemoved/', filename, "_.txt", sep =""), sep = "\t", row.names = FALSE)
    
    
  }
}

savereps <- heatmap_data %>% group_by(cell_line, treatment, doseLevel, dose_uM, variable, replID, cytokine, feature) %>%
  summarise(value = max(value))
tail(savereps)
write.table(file = '../resultaten/graphpad format/datawithReps.txt', savereps, sep = '\t', col.names = NA)
heatmap_data_time <- heatmap_data
head(heatmap_data_time)



# plot time courses, check of autofl of celdood invloed heeft. Hernormalizeer als mogelijk.

## eerst twee time heatmaps maken. Dan celdood data er weer bij doen. Dan de laatste concentratie-gegroepeerde HC maken. (heatmap_data)
unique(heatmap_data_time$variable)
##==##_HC1 create time course heatmap data of Count.Positives.2x.mean voor  Nrf2, Hmox1 & SRXN1, 
head(heatmap_data_time)
dir("../resultaten/figuren")
pdf(file = "../tmp/norm_validation.pdf", height =40, width = 32)
ggplot(data = subset(heatmap_data_time, grepl("mmn", variable)), aes(x = timeID, y = value)) + geom_point(aes(color = doseLevel)) + facet_grid(treatcytokine ~ variable, scales = "free")
dev.off()


getwd()

unique(heatmap_data_time$cell_line)
unique(heatmap_data_time$variable)
head(heatmap_data_time)
require(pheatmap)

HC1 <- heatmap_data_time %>% 
  filter(cell_line %in% c("Nrf2", "HMOX1", "SRXN1")) %>%
  group_by(cell_line, dose_uM, treatment, timeID, variable, doseLevel) %>%
  summarise(mean(value))

head(HC1)
unique(HC1$cell_line)
HC1$cell_line <- factor(HC1$cell_line, levels = c("Nrf2", "HMOX1", "SRXN1"))
HC1_wide <- dcast(treatment + doseLevel ~ cell_line + timeID, data = HC1, value.var = 'mean(value)')
head(HC1_wide)
lapply(HC1_wide,  function(x) sum(is.na(x)))
dim(HC1_wide)






myColors = list()
myColors$concentrationLevel <- brewer.pal(9, "BuGn" )
names(myColors$concentrationLevel) <- 1:9
# order columns

colnames(HC1_wide)

# only cluster row based treatment (leave dose-range)
myColors$cell_line <- c("purple", "light green", "light blue")

myColors$time <- colorRampPalette(brewer.pal(9, "Oranges"))(24)

  names(myColors$time) <- 1:24

names(myColors$cell_line) <- c("Nrf2", "HMOX1", "SRXN1")

concentrationtAnot<-data.frame(concentrationLevel = HC1_wide$doseLevel)

# color anotations for columns. 
colnames(HC1_wide)
annCol <- data.frame(cell_line=c( rep("Nrf2", 22 ), rep("HMOX1", 23), 
                          rep("SRXN1", 24)),
               time = c(1:22, 1:23, 1:24)
)

#colBreaks <-c(seq(-0.01, 0.3,length.out=15), seq(0.301,1.01, length.out = 14))

max(HC1_wide[, -c(1,2)])

#colBreaks <-c(  -seq(1.14, 0.411, length.out = 8), -seq( 0.4, 0.01, length.out=16)  , 
#                c(seq(0.01, 0.4,length.out=16), seq(0.411,1.01, length.out = 8)) )

my.colors <-brewer.pal(9, "GnBu")

pie(rep(1, length(my.colors)), labels = sprintf("%d (%s)", seq_along(my.colors), 
                                                my.colors), col = my.colors)
# cluster mean of treatments
head(HC1_wide)
aggrtreats <- HC1_wide %>% 
  group_by(treatment) %>% 
  select(-(treatment:doseLevel)) %>%
  summarise_all( mean)
aggrtreats <- data.frame(aggrtreats)
rownames(aggrtreats) <- aggrtreats$treatment
dist_hm_data <- dist(aggrtreats[,-1], method = c("euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski")[1])
clust_hm_data<-hclust(dist_hm_data, method = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")[1])

plot(clust_hm_data, hang = -1)
?plot.hclust
tmp <- table(HC1_wide$treatment)[clust_hm_data$order]
indl = alist()
for(i in seq_along(tmp)) {
  indl[[i]] <- which( HC1_wide$treatment %in% names(tmp)[i])
  
}
ind_rows<- unlist(indl)
HC1_wide$treatment[ind_rows]


rownames(HC1_wide) <- paste(HC1_wide$treatment, HC1_wide$doseLevel, sep = "_")

length(ind_rows)
nrow(HC1_wide)
?pheatmap

rownames(concentrationtAnot) <- rownames(HC1_wide)
rownames(annCol) <- colnames(HC1_wide)[-c(1,2)]


pdf(file = "../resultaten/figuren/ox_stress_time_heatmap.pdf", height =20, width = 16)
pheatmap(as.matrix(HC1_wide[ ind_rows , -c(1,2) ]), cluster_cols = FALSE, cluster_rows = FALSE, 
color = my.colors, border_color = "grey90",
fontsize = 12, width = 10, height=10, legend = TRUE,
 annotation_row = concentrationtAnot,
  annotation_col = annCol, annotation_colors = myColors   #,txt =hm_data_data[ indOrder ,colIndex]
)
dev.off()

##==##_HC2 create time course heatmap data of mmn_integrint voor  ICAM1 en A20. En alleen na TNF exposure
table(heatmap_data_time$cytokine)

HC2 <- heatmap_data_time %>% 
 # filter(cytokine == "TNF") %>%
  filter(cell_line %in% c("A20", "ICAM1")) %>%
  group_by(cell_line, cytokine, dose_uM, treatment, timeID, variable, doseLevel) %>%
  summarise(mean(value))

head(HC2)
unique(HC2$cell_line)
unique(HC2$cytokine)
HC2$cell_line <- factor(HC2$cell_line, levels = c("A20", "ICAM1"))
HC2$cytokine <- factor(HC2$cytokine, levels = c("DMEM", "TNF", "IL1"))
HC2_wide <- dcast(treatment + doseLevel ~ cell_line + cytokine + timeID, data = HC2, value.var = 'mean(value)')
head(HC2_wide)
lapply(HC2_wide,  function(x) sum(is.na(x)))
dim(HC2_wide)






myColors = list()
myColors$concentrationLevel <- colorRampPalette(brewer.pal(9, "BuGn" ))(10)
names(myColors$concentrationLevel) <- 1:10
# order columns

colnames(HC2_wide)

# only cluster row based treatment (leave dose-range)
myColors$cell_line <- c("purple", "light green")

myColors$cytokine <- brewer.pal(3, "Greys")

myColors$time <- colorRampPalette(brewer.pal(9, "Oranges"))(24)

names(myColors$time) <- 1:24
names(myColors$cytokine) <- c("DMEM", "TNF", "IL1")

names(myColors$cell_line) <- c("A20", "ICAM1")

concentrationtAnot<-data.frame(concentrationLevel = HC2_wide$doseLevel)

# color anotations for columns. 
colnames(HC2_wide)
dim(HC2_wide)
24*3*2 + 2
annCol <- data.frame(cell_line=c( rep("A20", 72 ), rep("ICAM1", 72)),
                     cytokine = rep(c( rep(rep("DMEM", 24), 1), rep(rep("TNF", 24), 1), rep(rep("IL1", 24), 1)), 2),
                     time = rep(1:24, 6)
)

#colBreaks <-c(seq(-0.01, 0.3,length.out=15), seq(0.301,1.01, length.out = 14))

max(HC2_wide[, -c(1,2)])

#colBreaks <-c(  -seq(1.14, 0.411, length.out = 8), -seq( 0.4, 0.01, length.out=16)  , 
#                c(seq(0.01, 0.4,length.out=16), seq(0.411,1.01, length.out = 8)) )

my.colors <-colorRampPalette(brewer.pal(9, "PuOr"))(99)

pie(rep(1, length(my.colors)), labels = sprintf("%d (%s)", seq_along(my.colors), 
                                                my.colors), col = my.colors)
# cluster mean of treatments
head(HC2_wide)
aggrtreats <- HC2_wide %>% 
  group_by(treatment) %>% 
  select(-(treatment:doseLevel)) %>%
  summarise_all( mean)
aggrtreats <- data.frame(aggrtreats)
rownames(aggrtreats) <- aggrtreats$treatment
dist_hm_data <- dist(aggrtreats[,-1], method = c("euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski")[1])
clust_hm_data<-hclust(dist_hm_data, method = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")[1])

plot(clust_hm_data, hang = -1)
?plot.hclust
tmp <- table(HC2_wide$treatment)[clust_hm_data$order]
indl = alist()
for(i in seq_along(tmp)) {
  indl[[i]] <- which( HC2_wide$treatment %in% names(tmp)[i])
  
}
ind_rows<- unlist(indl)
HC2_wide$treatment[ind_rows]

rownames(HC2_wide) <- paste(HC2_wide$treatment, HC2_wide$doseLevel, sep = "_")

length(ind_rows)
nrow(HC2_wide)

rownames(concentrationtAnot) <- rownames(HC2_wide)
rownames(annCol) <- colnames(HC2_wide)[-c(1,2)]

dim(HC2_wide)
length(HC2_wide[HC2_wide$treatment == "DMSO", ])
head(HC2_wide)
HC2_wide_dmso <- as.matrix(HC2_wide[, -c(1,2)]) - t(matrix(as.numeric(HC2_wide[HC2_wide$treatment == "DMSO", -c(1,2)]), ncol = 61, nrow = 144 ))

my_breaks <- seq(from = min(HC2_wide_dmso), to = abs(min(HC2_wide_dmso)), length.out = 100)
length(my.colors)



pdf(file = "../resultaten/figuren/infl_stress_time_heatmap_DMSOsubtracted_allcytokines.pdf", height =20, width = 32)
pheatmap(as.matrix(HC2_wide_dmso[ ind_rows ,  ]), cluster_cols = FALSE, cluster_rows = FALSE, 
         color = my.colors, breaks = my_breaks, border_color = NA,
         fontsize = 12, width = 10, height=10, legend = TRUE,
         annotation_row = concentrationtAnot,
         annotation_col = annCol, annotation_colors = myColors   #,txt =hm_data_data[ indOrder ,colIndex]
)
dev.off()


##==##_HC3  van de volgende reporters: Hmox1,  SRXN1, A20, ICAM1, BIP, CHOP, , + necrosis informatie van de wt : "ratio AreaPI/AreaNuclei WT cell death"+ "apoptoseinformatie v.d. wt: "ratio AreaAnV/AreaNuclei WT cell death"
# hier moet wel celdood data erbij, dus terug naar main.R en van metadata wisselen.
head(heatmap_data)
table(heatmap_data$timeID)
unique(heatmap_data$cell_line)
head(heatmap_data)
heatmap_data[heatmap_data$variable == "Apoptosis WT", "cell_line"] <- "HepG2 WT apoptosis"
heatmap_data[heatmap_data$variable == "Necrosis WT", "cell_line"] <- "HepG2 WT necrosis"
heatmap_data$feature <- paste(heatmap_data$cell_line, heatmap_data$cytokine)
head(heatmap_data)

unique(heatmap_data[heatmap_data$cell_line =="HepG2 WT necrosis",  "timeID"])

sum(is.na(heatmap_data_w))

#write.table(heatmap_data_w, file = '../resultaten/heatmap_tabel_doselevelRijen.txt', sep ="\t")

# hier verder 20042018: 2 heatmapjes maken ivm verschil in concentraties inflammatory reporters
ind_cols <- c("IkBalpha TNF", "IkBalpha IL1", "RelA TNF", "RelA IL1")
# seperate out the inflam and non-inflam dataset.

heatmap_data_rela_ikba <- heatmap_data[heatmap_data$cell_line %in% c("IkBalpha", "RelA", "HepG2 WT apoptosis", "HepG2 WT necrosis") , ]
heatmap_data_notInflam <- heatmap_data[ !heatmap_data$cell_line %in% c("IkBalpha", "RelA"), ]
head(heatmap_data_notInflam)
table(heatmap_data_notInflam$timeID)
#group_by(cell_line, treatment, doseLevel, variable, replID, cytokine, feature)
heatmap_data_notInflam <- heatmap_data_notInflam %>%
  filter(timeID %in%22) %>% group_by(cell_line, treatment, doseLevel, variable, cytokine, feature) %>%
  dplyr::summarise(meanValue = mean(value))


heatmap_data_notInflam_w <- dcast( treatment + doseLevel ~  feature, data = heatmap_data_notInflam)
head(heatmap_data_notInflam_w)



#heatmap_data <- heatmap_data[!heatmap_data$cytokine %in% c("IL1"), ]
head(heatmap_data)

unique(heatmap_data$variable)
unique(heatmap_data$feature)
unique(heatmap_data$variable)


heatmap_data_notInflam_w[which(is.na(heatmap_data_notInflam_w), arr.ind = T),, drop = F] 
write.table(heatmap_data_notInflam_w, file = '../resultaten/heatmap_tabel_doselevelRijen_noinflam.txt', sep ="\t")
# er mist data (8 datapunten, ), uitzoeken waar dit vandaan komt.

#heatmap_data_w[9, 11 ] <- 0.0315
#heatmap_data_w[62, 19 ] <- 0.0805

head(heatmap_data)


head(heatmap_data_w)

require(RColorBrewer)
unique(heatmap_data$doseLevel)
myColors = list()

myColors$concentrationLevel <- colorRampPalette(brewer.pal(9, "BuGn" ))(10)
names(myColors$concentrationLevel) <- 1:10
# order columns

# show no difference in cell death dependent on TNF in suppl figure
suppl_CD <- heatmap_data_notInflam_w[heatmap_data_notInflam_w$doseLevel == 10 , c("treatment" ,"HepG2 WT apoptosis DMEM", "HepG2 WT apoptosis TNF", "HepG2 WT necrosis DMEM", "HepG2 WT necrosis TNF")]

write.table(suppl_CD, file = "J:/Workgroups/FWN/LACDR/TOX/data steven wink/Unilever2/resultaten/small_tabel_CD_noIL1.txt", sep ="\t", col.names = TRUE, row.names = FALSE)

head(heatmap_data_notInflam_w)
colnames(heatmap_data_notInflam_w)

heatmap_data_notInflam_w <- heatmap_data_notInflam_w[ , c(1, 2, 18, 14, 19, 7, 6, 3, 5, 4, 15, 17, 16, 8, 10, 9, 11, 13, 12 )]





# only cluster row based treatment (leave dose-range)
myColors$cell_line <- c("light blue", "light green", "purple", "black")

names(myColors$cell_line) <- c("Oxidative-stress", "ER-stress", "Inflammatory-stress", "Cell Death")
myColors$cytokine <- brewer.pal(3, "Greys")
names(myColors$cytokine) <- c("DMEM", "TNF", "IL1")


concentrationtAnot<-data.frame(concentrationLevel = heatmap_data_notInflam_w$doseLevel)

# color anotations for columns. 
colnames(heatmap_data_notInflam_w)
annCol <- data.frame(feature=c( rep("Oxidative-stress",3 ), rep("ER-stress",2), 
                          rep("Inflammatory-stress", 6), rep("Cell Death", 6))
                 ,
               cytokine = c(rep("DMEM", 5), rep(c("DMEM", "TNF", "IL1" ), 4))
               )

#colBreaks <-c(seq(-0.01, 0.3,length.out=15), seq(0.301,1.01, length.out = 14))

min(heatmap_data_notInflam_w[, -c(1,2)])
head(heatmap_data_notInflam_w)
#colBreaks <-c(  -seq(1.14, 0.411, length.out = 8), -seq( 0.4, 0.01, length.out=16)  , 
#                c(seq(0.01, 0.4,length.out=16), seq(0.411,1.01, length.out = 8)) )

# normalize TNF data
heatmap_data_w_tnf_norm <- heatmap_data_notInflam_w

A20_TNF <- heatmap_data_w_tnf_norm  %>% 
  filter(treatment == "DMSO") %>% 
  select(`A20 TNF`)
heatmap_data_w_tnf_norm <- heatmap_data_w_tnf_norm %>% mutate( `A20 DMEM` = `A20 DMEM` - as.numeric(A20_TNF),
                                            `A20 TNF` = `A20 TNF` - as.numeric(A20_TNF),
                                            `A20 IL1` = `A20 IL1` - as.numeric(A20_TNF))

ICAM1_TNF <- heatmap_data_w_tnf_norm  %>% 
  filter(treatment == "DMSO") %>% 
  select(`ICAM1 TNF`)
heatmap_data_w_tnf_norm <- heatmap_data_w_tnf_norm %>% mutate( `ICAM1 DMEM`= `ICAM1 DMEM` - as.numeric(ICAM1_TNF),
                                                               `ICAM1 TNF` = `ICAM1 TNF` - as.numeric(ICAM1_TNF),
                                                               `ICAM1 IL1` = `ICAM1 IL1` - as.numeric(ICAM1_TNF))
head(heatmap_data_w_tnf_norm)



nrow(heatmap_data_w_tnf_norm)
head(heatmap_data_w_tnf_norm)


my.colors <-colorRampPalette(brewer.pal(9, "PuOr"))(99)
#my.colors <- c(rev(brewer.pal(9, "YlGnBu")) , brewer.pal(9, "YlOrRd"))


pie(rep(1, length(my.colors)), labels = sprintf("%d (%s)", seq_along(my.colors), 
                                                my.colors), col = my.colors)
colnames(heatmap_data_w_tnf_norm)
# cluster mean of treatments
aggrtreats <- heatmap_data_w_tnf_norm %>% 
  group_by(treatment) %>% 
  select(-(treatment:doseLevel)) %>%
  summarise_all( mean)
aggrtreats <- as.data.frame(aggrtreats)
rownames(aggrtreats) <- aggrtreats$treatment

dist_hm_data <- dist(aggrtreats[,-1], method = c("euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski")[1])
clust_hm_data<-hclust(dist_hm_data, method = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")[1])
table(heatmap_data_w$treatment)
plot(clust_hm_data, hang = -1)
tmp <- table(heatmap_data_w_tnf_norm$treatment)[clust_hm_data$order]
indl = alist()
for(i in seq_along(tmp)) {
indl[[i]] <- which( heatmap_data_w_tnf_norm$treatment %in% names(tmp)[i])

}
ind_rows<- unlist(indl)
heatmap_data_w_tnf_norm$treatment[ind_rows]
table(heatmap_data_w_tnf_norm$treatment)

length(ind_rows)
nrow(heatmap_data_w_tnf_norm)

library(dendextend) # for reorder dendogram (if needed)

head(heatmap_data_w_tnf_norm)


#cutree(dd.reorder,h=7) # it now works on a dendrogram#
# It is like using:
#dend1 <-dendextend:::cutree.dendrogram(dd.reorder,h=7)

# dend1 <- color_branches(dd, k = 4)
# dend1 <- color_labels(dend1, k = 4)
# plot(dend1)

# pdf("../generated/results/heatmap/orig_hm_rowClusterDendogram_manh_wardD_noreo.pdf", width = 24, height = 10)
# plot(dend1)
# dev.off()



head(heatmap_data_w_tnf_norm)

rownames(heatmap_data_w_tnf_norm) <- paste(heatmap_data_w_tnf_norm$treatment, heatmap_data_w_tnf_norm$doseLevel, sep ="_")

rownames(concentrationtAnot) <- rownames(heatmap_data_w_tnf_norm)
rownames(annCol) <- colnames(heatmap_data_w_tnf_norm)[-c(1,2)]

my_breaks <- seq(from = min(heatmap_data_w_tnf_norm[, -c(1,2)]), 
                 to = (max(heatmap_data_w_tnf_norm[, -c(1,2)])), length.out = 100)



pdf("../resultaten/figuren/large_heatmap_eucl_wardD1.pdf", width= 12, height = 16)
pheatmap(as.matrix(heatmap_data_w_tnf_norm[ ind_rows , -c(1,2) ]), cluster_cols = FALSE, cluster_rows = FALSE, 
         color = my.colors, border_color = "grey90",
         fontsize = 12, width = 10, height=10, legend = TRUE,
         annotation_row = concentrationtAnot,
         annotation_col = annCol, annotation_colors = myColors 
         )
dev.off()

##===### HC4
# repeat but now for ikba and rela (less tp)
heatmap_data_rela_ikba <- heatmap_data[heatmap_data$cell_line %in% c("IkBalpha", "RelA", "HepG2 WT apoptosis", "HepG2 WT necrosis") , ]

head(heatmap_data_rela_ikba)
table(heatmap_data_rela_ikba$timeID)
#

heatmap_data_rela_ikba <- heatmap_data_rela_ikba %>% 
  group_by(cell_line, treatment, doseLevel, variable, replID, cytokine, feature) %>%
  dplyr::summarise(max_value = max(value)) %>% group_by(cell_line, treatment, doseLevel, variable, cytokine, feature) %>%
  dplyr::summarise(meanValue = mean(max_value))


heatmap_data_rela_ikba_w <- dcast( treatment + doseLevel ~  feature, data = heatmap_data_rela_ikba)
head(heatmap_data_rela_ikba)



#heatmap_data <- heatmap_data[!heatmap_data$cytokine %in% c("IL1"), ]
head(heatmap_data_rela_ikba)

unique(heatmap_data_rela_ikba$variable)
unique(heatmap_data_rela_ikba$feature)
unique(heatmap_data_rela_ikba$variable)

heatmap_data_rela_ikba_w <- na.omit(heatmap_data_rela_ikba_w)


heatmap_data_rela_ikba_w[which(is.na(heatmap_data_rela_ikba_w), arr.ind = T),, drop = F] 
write.table(heatmap_data_rela_ikba_w, file = '../resultaten/heatmap_tabel_doselevelRijen_inflam.txt', sep ="\t")
# er mist data (8 datapunten, ), uitzoeken waar dit vandaan komt.

#heatmap_data_w[9, 11 ] <- 0.0315
#heatmap_data_w[62, 19 ] <- 0.0805

head(heatmap_data)


head(heatmap_data_w)

require(RColorBrewer)
unique(heatmap_data$doseLevel)
myColors = list()

(unique(heatmap_data_rela_ikba_w$doseLevel))
heatmap_data_rela_ikba_w$doseLevel<- factor(heatmap_data_rela_ikba_w$doseLevel)

myColors$concentrationLevel <- colorRampPalette(brewer.pal(9, "BuGn" ))(10)

names(myColors$concentrationLevel) <- 1:10
# order columns

# show no difference in cell death dependent on TNF in suppl figure
suppl_CD <- heatmap_data_notInflam_w[heatmap_data_notInflam_w$doseLevel == 10 , c("treatment" ,"HepG2 WT apoptosis DMEM", "HepG2 WT apoptosis TNF", "HepG2 WT necrosis DMEM", "HepG2 WT necrosis TNF")]

write.table(suppl_CD, file = "J:/Workgroups/FWN/LACDR/TOX/data steven wink/Unilever2/resultaten/small_tabel_CD_noIL1.txt", sep ="\t", col.names = TRUE, row.names = FALSE)

head(heatmap_data_rela_ikba_w)
colnames(heatmap_data_rela_ikba_w)

heatmap_data_rela_ikba_w <- heatmap_data_rela_ikba_w[ , c(1, 2, 10, 9, 12, 11, 3, 5, 4, 6, 8, 7  )]


# only cluster row based treatment (leave dose-range)
myColors$cell_line <- c( "purple", "black")

names(myColors$cell_line) <- c("Inflammatory-stress", "Cell Death")
myColors$cytokine <- brewer.pal(3, "Greys")
names(myColors$cytokine) <- c("DMEM", "TNF", "IL1")


concentrationtAnot<-data.frame(concentrationLevel = heatmap_data_rela_ikba_w$doseLevel)

# color anotations for columns. 
colnames(heatmap_data_rela_ikba_w)
annCol <- data.frame(feature=c( rep("Inflammatory-stress", 4), 
                                rep("Cell Death", 6))
                     ,
                     cytokine = c(rep(c("TNF", "IL1"), 2), rep(c("DMEM", "TNF", "IL1" ), 2))
)

#colBreaks <-c(seq(-0.01, 0.3,length.out=15), seq(0.301,1.01, length.out = 14))

min(heatmap_data_rela_ikba_w[, -c(1,2)])
max(heatmap_data_rela_ikba_w[, -c(1,2)])


#colBreaks <-c(  -seq(1.14, 0.411, length.out = 8), -seq( 0.4, 0.01, length.out=16)  , 
#                c(seq(0.01, 0.4,length.out=16), seq(0.411,1.01, length.out = 8)) )

# normalize TNF data

my.colors <-colorRampPalette(brewer.pal(9, "GnBu"))(99)
#my.colors <- c(rev(brewer.pal(9, "YlGnBu")) , brewer.pal(9, "YlOrRd"))


pie(rep(1, length(my.colors)), labels = sprintf("%d (%s)", seq_along(my.colors), 
                                                my.colors), col = my.colors)
colnames(heatmap_data_rela_ikba_w)
# cluster mean of treatments
aggrtreats <- heatmap_data_rela_ikba_w %>% 
  group_by(treatment) %>% 
  select(-(treatment:doseLevel)) %>%
  summarise_all( mean)
aggrtreats <- as.data.frame(aggrtreats)
rownames(aggrtreats) <- aggrtreats$treatment

dist_hm_data <- dist(aggrtreats[,-1], method = c("euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski")[1])
clust_hm_data<-hclust(dist_hm_data, method = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")[1])
table(heatmap_data_w$treatment)
plot(clust_hm_data, hang = -1)
tmp <- table(heatmap_data_rela_ikba_w$treatment)[clust_hm_data$order]
indl = alist()
for(i in seq_along(tmp)) {
  indl[[i]] <- which( heatmap_data_rela_ikba_w$treatment %in% names(tmp)[i])
  
}
ind_rows<- unlist(indl)
heatmap_data_rela_ikba_w$treatment[ind_rows]
table(heatmap_data_rela_ikba_w$treatment)

length(ind_rows)
nrow(heatmap_data_rela_ikba_w)

library(dendextend) # for reorder dendogram (if needed)

head(heatmap_data_rela_ikba_w)


#cutree(dd.reorder,h=7) # it now works on a dendrogram#
# It is like using:
#dend1 <-dendextend:::cutree.dendrogram(dd.reorder,h=7)

# dend1 <- color_branches(dd, k = 4)
# dend1 <- color_labels(dend1, k = 4)
# plot(dend1)

# pdf("../generated/results/heatmap/orig_hm_rowClusterDendogram_manh_wardD_noreo.pdf", width = 24, height = 10)
# plot(dend1)
# dev.off()



head(heatmap_data_rela_ikba_w)

rownames(heatmap_data_rela_ikba_w) <- paste(heatmap_data_rela_ikba_w$treatment, heatmap_data_rela_ikba_w$doseLevel, sep ="_")

rownames(concentrationtAnot) <- rownames(heatmap_data_rela_ikba_w)
rownames(annCol) <- colnames(heatmap_data_rela_ikba_w)[-c(1,2)]

my_breaks <- seq(from = min(heatmap_data_rela_ikba_w[, -c(1,2)]), 
                 to = (max(heatmap_data_rela_ikba_w[, -c(1,2)])), length.out = 100)



pdf("../resultaten/figuren/small_heatmap_eucl_wardD1.pdf", width= 12, height = 16)
pheatmap(as.matrix(heatmap_data_rela_ikba_w[ ind_rows , -c(1,2) ]), cluster_cols = FALSE, cluster_rows = FALSE, 
         color = my.colors, border_color = "grey90",
         fontsize = 12, width = 10, height=10, legend = TRUE,
         annotation_row = concentrationtAnot,
         annotation_col = annCol, annotation_colors = myColors 
)
dev.off()
