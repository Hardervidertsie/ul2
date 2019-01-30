

## Data analysis pipeline BioSpyder datasets ##

# 0 - Data assembly, normalization, QCs

# Marije Niemeijer
# 06/07/2018

# Version 1

## Steps:
# Loading of data
# Library size calculation and filtering
# Normalization of raw counts using DESEq2 package
# Quality control plots (library size and count distribution, replicate correlation, variance across means)
# Exluding specific samples for further analysis (NoCells and MAQC controls, and others if necessary)
# DEGs determination: Log2FC and p-value calculation using DESEq2 package
# Filtering of DEGs (based on p-value/Log2FC)
# Preparation of txt files for IPA as input files
# Assessing gene expression for specific genes
# Principal component analysis (PCA plots, top25 genes in PCs, contribution of variables)
# Hierarchical clustering of samples and genes in heatmap
# Saving data

#############################################

## Settings ##
rm(list=ls(all=TRUE))
options(stringsAsFactors = FALSE)
options(scipen=999)

## Libraries ##

library("gridExtra")
library("stringr")
library("ggplot2")
library("pheatmap")
library("reshape2")
library("RColorBrewer")
library("plyr")
library("dplyr")
library("tidyr")
library("colorspace")
library("scales")
library("data.table")

library("compare")
library("DESeq2")
library("readxl")
library("PoiClaClu")
library("hexbin")
library("ggalt")

#biocLite()
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggbiplot")

## INPUT USER ##


MetaName <- "UniL2_metadata.txt" #use template and save as txt file
CountThreshold <- 100000 #Exclude samples lower than this number for library size
VariableQCs <- "SAMPLE_ID" #Variable to check distribution Library size and counts (choose CELL_ID, TREATMENT, SAMPLE_ID, EXP_ID, TIMEPOINT or REPLICATE)
RemoveSamples <- c("MAQC", "NoCells") #Give list of SAMPLE_IDs for samples which you would like to remove for further analysis (for example MAQC, NoCells samples, outliers)
ControlSample <- 'HepG2_HepG2-WT_medium_conc0_TP8' # Sample name (meanID, same for replicates) of control for log2FC and p-value calculations
                                                       # meanID: 'EXP_ID'_'CELL_ID'_'TREATMENT'_conc'CONCENTRATION'_TP'TIMEPOINT

Threshold_padj <- 0.1 #Filtering of genes with p-adj lower than threshold
Threshold_log2FC <- 0.5 #Filtering of genes with absolute log2FC higher than threshold (genes lower than -0.5, and genes higher than 0.5)
Filtering <- "padj" #Fill in filtering method ("padj", "log2FC" or "pajd&log2FC")
CheckGenes <- c("^NFE2L2_", "^SRXN1_", "^HMOX1_", "^NQO1_", "^NFKBIA_", "^ICAM1_", "^TNFAIP3_", "^BIRC3_") #Check expression of individual genes among samples (use following template: "^GenesymbolX_")
                                                                    #You can list as many genes as you like, however graph may become too full
Geneset <- "UPR_TOP50genes_TempoSeq.txt" #Fill in file name of specific gene list of interest (GeneSymbols, txt format)


# Working directories #

CountTableDir <- "J:/Workgroups/FWN/LACDR/TOX/Data Bas/Post-Doc/Biospyder/Unilever2/raw counts" #Folder containing only raw count table files
inputDir <-  'J:/Workgroups/FWN/LACDR/TOX/Data Bas/Post-Doc/Biospyder/Unilever2/metadata/' #Store metadata file and raw count files here
outputDir_GR <- 'J:/Workgroups/FWN/LACDR/TOX/Data Bas/Post-Doc/Biospyder/Unilever2/Graphs/'
GeneListDir <- 'J:/Workgroups/FWN/LACDR/TOX/Data Bas/Post-Doc/Biospyder/Unilever2/Genelists/' #store specific gene lists as txt files here
output_IPA <- 'J:/Workgroups/FWN/LACDR/TOX/Data Bas/Post-Doc/Biospyder/Unilever2/IPA_input/' #Directory for txt files to use as IPA input 

## Functions ##



LoadCountData <- function(FilePath) {
  setwd(FilePath)
  CountData <- list()
  
  for(i in 1:length(dir(FilePath))){
    x <- dir(FilePath)[i]
    if(grepl("xlsx", x)){
      tmp  <- read_excel(x, sheet = 1)
    } else {
      tmp  <- read.csv(x, sep = "\t")
    }
    
    CountData[i] <- list(as.data.frame(tmp))
  }  
  
  probeNames <- CountData[[1]][[1]]
  
  for(i in 1:length(CountData)) {
    rownames(CountData[[i]]) <- probeNames
    if (!isTRUE(all.equal(rownames(CountData[[i]]),CountData[[i]][[1]]))) { 
      mismatches <- paste(which(rownames(CountData[[i]]) != CountData[[i]][[1]]), collapse = ",")
      print(paste0("ProbeIDs do not match between count tables: ", i)); print(mismatches)
    }
    
    CountData[[i]] <- CountData[[i]][-1]
  }
  CountData <- do.call(cbind, CountData)
  CountData
}

getProbeName <- function(gene, data = log2Norm) return(grep(toupper(gene), rownames(data), value = TRUE))

PathwayProbesCounts <- function(Genelist = "OX_genes.txt", CountTable = df, probes = TRUE, matrix = FALSE, GenelistType = "Txtfile"){
  if(GenelistType == "Txtfile"){
    Genes <- read.table(paste0(GeneListDir, Genelist), sep = "\t", header = TRUE)	
    GenesOfInterest <- Genes[,1]
  } else {
    GenesOfInterest <- Genelist
  }
  
  if (probes == FALSE){
    GenesOfInterest  <- unique(unlist(sapply(GenesOfInterest, getProbeName, data = log2Norm)))
  }
  if (matrix == TRUE){
    df_sub <- CountTable[which(grepl(paste0(GenesOfInterest, collapse = "|"), row.names(CountTable))),]
  } else {
    df_sub <- CountTable[which(grepl(paste0(GenesOfInterest, collapse = "|"), CountTable$Probe_ID)),]
  }
  return(df_sub)
}

runDESeq2_DEGs <- function(rawCounts, meta, 
                           meanID_control = 'UniL2_HepG2_DMSO NO_concC0_TP8', #Sample name (meanID) of control which you want to use to calculate fold changes (not unique for different replicates)
                           meanID_column = "meanID") {
  
  rslt <- c()
  
  for(c in 1:length(unique(meta[,c(meanID_column)]))){
    
    sample <- unique(meta[,c(meanID_column)])[c]
    
    samples <- c(rownames(meta[meta[,c(meanID_column)] == meanID_control | 
                                 meta[,c(meanID_column)] == sample, ]))
    
    count_tmp <- rawCounts[, samples]
    meta_tmp <- meta[samples, ]
    meta_tmp$CONTRAST <- factor(meta_tmp[,c(meanID_column)] != meanID_control)
    
    if(length(unique(meta_tmp$CONTRAST))==2){
      dds_tmp <- DESeqDataSetFromMatrix(countData = count_tmp, colData = meta_tmp, design = ~ CONTRAST)
      dds_tmp <- estimateSizeFactors(dds_tmp)
      dds_tmp <- DESeq(dds_tmp)
      
      rslt_tmp <- as.data.frame(results(dds_tmp))
      rslt_tmp[,c(meanID_column)] <- sample
      rslt_tmp$Probe_ID <- rownames(rslt_tmp)
      rownames(rslt_tmp) <- c()
      
      rslt <- rbind(rslt, rslt_tmp)
    }
    
    if(c %% 10==0){
      print(c)
    }
    
  }
  return(rslt)
} #Column name of the sample names which are not unique for different replicates

PCADF <- function(DF = log2Norm){
  pca <- prcomp(t(DF))
  return(pca)
}

PCAplots <- function(pca = df_pca, metaDF = meta_Filtered, colourVar = "TREATMENT", #Fill in variable for which you want different colors in plot
                     Clusts = FALSE, NrClusts = 3, 
                     FileID = "Allsamples", #Unique name to save plot
                     width = 10, height = 8, #Size of plot
                     sizeP = 0.2, #size of points in plot
                     Dataclass = "log2", #Fill in log2 or log2FC
                     a = 0.7){ #Transparency of dots in plot
  
  if(Dataclass == "log2"){
    Vari <- c("SAMPLE_ID", "CELL_ID", "TREATMENT", "CONCENTRATION", "TIMEPOINT", "REPLICATE", "EXP_ID", "LIB_SIZE")
    ID <- "SAMPLE_ID"
  }
  
  if(Dataclass == "log2FC"){
    Vari <- c("meanID", "CELL_ID", "TREATMENT", "CONCENTRATION", "TIMEPOINT", "EXP_ID")
    ID <- "meanID"
  }

  pca_eigs <- pca$sdev^2
  pca_proportion <- pca_eigs / sum(pca_eigs)
  pca_proportion <- data.frame("variable" = c("PC1", "PC2", "PC3", "PC4", "PC5"), "proportion" = pca_proportion[1:5])
  pca_clust <- melt(kmeans(pca$x[,1:5], centers = NrClusts)$cluster)
  pca_clust[,ID] <- rownames(pca_clust)
  colnames(pca_clust)[1] <- "clust"
  
  pca_df <- data.frame(pca$x[which(grepl(paste(metaDF[,ID], collapse = "|"), rownames(pca$x))), 1:5])
  pca_df <- data.frame(metaDF[rownames(pca_df),], pca_df)
  pca_Long <- melt(pca_df, id.vars = Vari, measure.vars = grep("PC", colnames(pca_df), value = TRUE))
  pca_Long <- left_join(pca_Long, pca_proportion)
  pca_Long$PlotLabels <- paste0(pca_Long$variable, " (Variance: ", 100*round(pca_Long$proportion, digits = 3), "%)")
  pca_Long <- left_join(pca_Long, pca_clust)
  
  pca_GR <- do.call(rbind, lapply(1:5, function(i) {
    pca_Long$x <-   pca_Long$value[pca_Long$variable == paste0('PC', i)]
    pca_Long$xPC <- factor(paste0('PC', i), ordered = TRUE)
    pca_Long$yPC <- factor(pca_Long$variable, ordered = TRUE)
    pca_Long
  })) 
  
  levels(pca_GR$xPC) <- unique(pca_GR$PlotLabels) 
  
  pca_GR$colourVar <- pca_GR[,which(names(pca_GR)==colourVar)]
  
  if(Clusts == TRUE){
    cl <- "Clusts"
    pdf(paste0(outputDir_GR, "PCA_", FileID, "_", colourVar, "_", cl,".pdf"), width = width, height = height)
    print(ggplot(pca_GR, aes(x = x, y = value)) +
            geom_encircle(aes(group = factor(clust), fill = factor(clust)), s_shape = 0.7, expand = 0, color = "#666666", alpha = a) +    
            geom_point(aes(colour = colourVar), alpha = a, size = sizeP) +
            facet_grid(yPC ~ xPC, scales = 'free') +
            theme_bw() +
            theme(panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank()) +
            labs(x = "", 
                 y = "", 
                 title = paste0("PCA_", FileID, "_", colourVar, "_", cl),
                 colour = paste0(colourVar)))
    dev.off()
    
  } else { 
    cl <- "NoClusts"
    pdf(paste0(outputDir_GR, "PCA_", FileID, "_", colourVar, "_", cl,".pdf"), width = width, height = height)
    print(ggplot(pca_GR, aes(x = x, y = value)) +
            geom_point(aes(colour = colourVar), alpha = a, size = sizeP) +
            facet_grid(yPC ~ xPC, scales = 'free') +
            theme_bw() +
            theme(panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank()) +
            labs(x = "", 
                 y = "", 
                 title = paste0("PCA_", FileID, "_", colourVar, "_", cl),
                 colour = paste0(colourVar)))
    dev.off()
  }
}

PCAcontrVar <- function(pcaDF = df_pca){
  tmp <- data.frame(PC1Genes = names(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,1])))[1:25],1]),
                    PC1rotation = unname(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,1])))[1:25],1]),
                    PC2Genes = names(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,2])))[1:25],2]),
                    PC2rotation = unname(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,2])))[1:25],2]),
                    PC3Genes = names(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,3])))[1:25],3]),
                    PC3rotation = unname(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,3])))[1:25],3]),
                    PC4Genes = names(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,4])))[1:25],4]),
                    PC4rotation = unname(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,4])))[1:25],4]),
                    PC5Genes = names(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,5])))[1:25],5]),
                    PC5rotation = unname(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,5])))[1:25],5]))
  return(tmp)
}

pcaCoV <- function(pcaDF = df_pca, metaDF = meta_Filtered, filename = "Allsamples"){
  
  pcaDF <- data.frame(metaDF[rownames(pcaDF$x),], pcaDF$x[, 1:5])
  pcaDF[, which(grepl("character|logical|factor", sapply(pcaDF, class)))] <- apply(pcaDF[,which(grepl("character|logical|factor", sapply(pcaDF, class)))], 2, function(x) as.numeric(as.factor(x)))
  pearsonMatrix <- na.omit(cor(pcaDF, method = 'pearson')[grep('PC', colnames(pcaDF), value = TRUE, invert = TRUE), 
                                                          grep('PC', colnames(pcaDF), value = TRUE)])
  
  pheatmap(pearsonMatrix, main = paste0('Pearson correlation -', filename),
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           breaks = c(seq((-1 * 1), 0, length.out = ceiling(100/2) + 1), 
                      seq(1/100, 1, length.out = floor(100/2))),
           cellwidth = 25, 
           cellheight = 15,
           width = 7,
           file = paste0(outputDir_GR, "Heatmap_PCACoV_", filename, ".pdf"))
}

Heatmap <- function(Data = log2Norm, Meta = meta_Filtered, 
                    FC = FALSE,
                    heightGene = 0.2,
                    widthGenes = 1,
                    FileName = "Allsamples",
                    GeneName = FALSE) {
  
    Meta <- Meta[colnames(Data),]
    
    if(FC == FALSE){
      breaksList <- seq(min(Data, na.rm = TRUE)*1.1, max(Data, na.rm = TRUE)*1.1, by = 0.01)
      meta_Filtered_ann <- Meta[, c("EXP_ID", "CELL_ID", "REPLICATE", "TREATMENT", "CONCENTRATION", "TIMEPOINT")]
      hm_color <- rev(heat_hcl(length(breaksList), c = c(80, 30), l = c(30, 90), power = c(1/5, 1.3)))
    } else {
      meta_Filtered_ann <- Meta[, c("EXP_ID", "CELL_ID", "TREATMENT", "CONCENTRATION", "TIMEPOINT")]
      Data[is.na(Data)] <- 0
      breaksList <- c(seq(min(Data), -0.5, length.out=50), seq(-0.5, 0.5,length.out=50)[-1],
                      seq(0.5, max(Data),length.out=50)[-1])
      hm_color <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(150))
    }
  
    meta_Filtered_ann[] <- lapply(meta_Filtered_ann, factor)
    df_hm <- Data[,rownames(meta_Filtered_ann)]
    
    df_hm <- df_hm[apply(df_hm, MARGIN = 1, FUN = function(x) sd(x, na.rm = TRUE) != 0),]
    df_hm <- df_hm[, apply(df_hm, MARGIN = 2, FUN = function(x) sd(x, na.rm = TRUE) != 0)]
    df_hm <- df_hm[, apply(df_hm, MARGIN = 2, FUN = function(x) all(is.finite(x)))]
      
    pheatmap(as.matrix(df_hm),
                      color = hm_color,
                      breaks = breaksList,
                      clustering_distance_rows = "correlation",
                      clustering_distance_cols = "correlation",
                      clustering_method = "complete", 
                      annotation_col = meta_Filtered_ann,
                      display_numbers = FALSE, 
                      fontsize_row = 6,
                      fontsize_col = 6,
                      fontsize_number = 3, 
                      border_color = NA, 
                      cluster_cols = TRUE,
                      show_colnames = TRUE,
                      show_rownames = GeneName,
                      cellwidth = widthGenes,
                      cellheight = heightGene,
                      treeheight_row = 25,
                      treeheight_col = 30,
                      main = FileName,
                      filename = paste0(outputDir_GR, paste("Heatmap_", FileName, ".pdf")))
}


## Assembly of Data ##

# Load data #
meta <- read.delim(paste0(inputDir, MetaName), stringsAsFactors = FALSE)


rownames(meta) <- meta[,1]
count <- LoadCountData(CountTableDir)


meta$meanID <- paste0(meta$EXP_ID, "_",
                      meta$CELL_ID, "_",
                      meta$TREATMENT, "_conc",
                      meta$CONCENTRATION,"_TP", 
                      meta$TIMEPOINT)

meta$SAMPLENAME <- paste0(meta$EXP_ID, "_",
                          meta$CELL_ID, "_",
                          meta$TREATMENT, "_TP", 
                          meta$TIMEPOINT, "_conc",
                          meta$CONCENTRATION, "_rep",
                          meta$REPLICATE)

if(!all(rownames(meta) == colnames(count))){
  warning("No identical sample names (meta & counts)")
} 

ind <- (rownames(meta) %in%
          colnames(count))


count <- count[, rownames(meta)]

## Library size filter ##

meta$LIB_SIZE <- colSums(count, na.rm = TRUE)
summary(meta$LIB_SIZE)
meta_LowLibSize <- meta[which(meta$LIB_SIZE < CountThreshold), ]
CELLID_LowLibSize <- as.data.frame(table(meta_LowLibSize $CELL_ID))
meta_Filtered <- meta[which(meta$LIB_SIZE > CountThreshold), ] 
count_Filtered <- count[, rownames(meta_Filtered)]

## Check if CountTable contains NAs ##

count_NA <- count_Filtered[rowSums(is.na(count_Filtered)) > 0, colSums(is.na(count_Filtered)) > 0]
if(nrow(count_NA) == 0){
  print("No NAs in counts table")
} else {
  warning("NAs present in counts table")
  print(table(count_NA)); print(dim(count_NA))
}

count_Filtered <- na.omit(count_Filtered)

## Base normalization for normalized countdata using DESeq2 ##
unique(gsub( "(_[0-9]{1,30})$", "" ,rownames(count1)))

dds <- DESeqDataSetFromMatrix(countData = count_Filtered, 
                              colData = meta_Filtered, 
                              design = ~ 1)


dds <- estimateSizeFactors(dds)
Norm <- as.data.frame(counts(dds, normalized = TRUE)) #Normalized by size factors
log2Norm <- log2(Norm + 1) #Log2 normalization of size factor corrected counts

## Quality control of data ##

# Dataframe for plots

count_long <- melt(as.matrix(count))
colnames(count_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
count_long <- left_join(count_long, meta[,c("SAMPLE_ID", "SAMPLENAME", "CELL_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID")])

Norm_long <- melt(as.matrix(Norm))
colnames(Norm_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
Norm_long <- left_join(Norm_long, meta_Filtered[,c("SAMPLE_ID", "SAMPLENAME", "CELL_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID")])

log2Norm_long <- melt(as.matrix(log2Norm))
colnames(log2Norm_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
log2Norm_long <- left_join(log2Norm_long, meta_Filtered[,c("SAMPLE_ID", "SAMPLENAME", "CELL_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID")])

# Library size distribution #

meta$VariableQCs = meta[,which(names(meta)==VariableQCs)]
pdf(paste0(outputDir_GR, "LibSize_", VariableQCs, ".pdf"), width = 24, height = 8)
print(ggplot(meta, aes(x = reorder(VariableQCs, LIB_SIZE, FUN = median), y = LIB_SIZE)) +
        geom_boxplot(size = 0.3, outlier.size = 0.5) +
        scale_y_log10(limits = c(1, max(meta$LIB_SIZE))) +
        theme_classic() +
        theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
              axis.title.x = element_text(size=12, vjust = 0.25),
              axis.title.y = element_text(size=12, vjust = 1),
              axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
              axis.text.y = element_text(size=12)) +
        ggtitle("Library size distribution") + ylab('Library size') + xlab(VariableQCs))
dev.off()

# Counts & normalized counts distribution amoung CELL IDs #

count_long$VariableQCs = count_long[,which(names(count_long)==VariableQCs)]
pdf(paste0(outputDir_GR, "Distribution-Counts_", VariableQCs, ".pdf"), width = 24, height = 8)
print(ggplot(count_long, aes(x = reorder(VariableQCs, COUNTS, FUN = median), y = COUNTS+1)) +
        geom_boxplot(size = 0.3, outlier.size = 0.5) +
        scale_y_log10(limits = c(1, max(count_long$COUNTS))) +
        theme_classic() +
        theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
              axis.title.x = element_text(size=12, vjust = 0.25),
              axis.title.y = element_text(size=12, vjust = 1),
              axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
              axis.text.y = element_text(size=12)) +
        ggtitle("Distribution counts") + ylab('counts') + xlab(VariableQCs))
dev.off()

log2Norm_long$VariableQCs = log2Norm_long[,which(names(log2Norm_long)==VariableQCs)]
pdf(paste0(outputDir_GR, "Distribution-CountsNorm_", VariableQCs, ".pdf"), width = 24, height = 8)
print(ggplot(log2Norm_long, aes(x = reorder(VariableQCs, COUNTS, FUN = median), y = COUNTS)) +
        geom_boxplot(size = 0.3, outlier.size = 0.5) +
        theme_classic() +
        theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
              axis.title.x = element_text(size=12, vjust = 0.25),
              axis.title.y = element_text(size=12, vjust = 1),
              axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
              axis.text.y = element_text(size=12)) +
        ggtitle("Distribution Normalized counts") + ylab('Normalized counts') + xlab(VariableQCs))
dev.off()

# SDs across means after different normalization methods

notAllZero <- (rowSums(counts(dds))>0)

require("vsn")

pdf(paste0(outputDir_GR, "SDvsMean_log2Norm.pdf"), width = 24, height = 8)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))
dev.off()

# Replicate correlation #

repCors <- sapply(unique(meta_Filtered[, "meanID"]), function(expID) {
  replicateIDs <- rownames(meta_Filtered)[which(meta_Filtered[, "meanID"] == expID)]
  if(length(replicateIDs)>1){
    experimentMean <- rowMeans(log2Norm[, replicateIDs])
    apply(log2Norm[, replicateIDs], 2, function(.expID) cor(.expID, experimentMean))
  }
})

repCors <- data.frame(meanID = as.factor(rep(names(repCors), do.call(rbind, lapply(repCors, function(x) length(x))))),
                      SAMPLE_ID  = unlist(sapply(repCors, names)),
                      pearsonR     = unlist(repCors))
repCors$clr <- repCors$pearsonR > 0.95
repCors <- left_join(repCors, meta_Filtered)

repCors$VariableQCs = repCors[,which(names(repCors)==VariableQCs)]
pdf(paste0(outputDir_GR, "RepCor_PearsonR_", VariableQCs, ".pdf"), width = 24, height = 8)
ggplot(repCors, aes(x = reorder(VariableQCs, pearsonR, FUN = min), y = pearsonR)) +
  geom_hline(yintercept = 0.95, colour = 'grey75', size = 0.75) +
  geom_point(aes(color = clr), alpha = 0.7, size = 0.8) +
  theme_classic() +
  theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
        axis.title.x = element_text(size=12, vjust = 0.25),
        axis.title.y = element_text(size=12, vjust = 1),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=12),
        legend.position = 'none') +
  ggtitle("Replicate correlation") + ylab('PearsonR') + xlab(VariableQCs)
dev.off()

repCors_lowPR <- filter(repCors, clr == FALSE)

## Excluding specific samples for further analysis

meta_Filtered <- meta_Filtered[which(!grepl(paste0(RemoveSamples, collapse = "|"), meta_Filtered$SAMPLE_ID)),]
count_Filtered <- count_Filtered[, rownames(meta_Filtered)]
log2Norm <- log2Norm[, rownames(meta_Filtered)]
Norm <- Norm[, rownames(meta_Filtered)]

## Log2FC & pvalues with DESEq2 - DEGs ##

unique(colnames(count_Filtered))
unique((meta_Filtered$meanID))[grepl("DMSO", unique(meta_Filtered$meanID))]


ncol(count_Filtered)
nrow(meta_Filtered)

# select relevant subsets
#
ControlSample <- "UniL2_HepG2_DMSO NO_concC0_TP8"
ind_meta_hepg2_no <- grepl("HepG2",meta_Filtered$SAMPLE_ID) & grepl("noTNF", meta_Filtered$SAMPLE_ID)
ind_count_hepg2_no <- grepl("HepG2",colnames(count_Filtered)) & grepl("noTNF", colnames(count_Filtered))

identical(ind_meta_hepg2_no, ind_count_hepg2_no)

# 
ControlSample <- "UniL2_HepG2_DMSO YES_concC0_TP8"
ind_meta_hepg2_yes <- grepl("HepG2",meta_Filtered$SAMPLE_ID) & grepl("_TNF_", meta_Filtered$SAMPLE_ID)
ind_count_hepg2_yes <- grepl("HepG2",colnames(count_Filtered)) & grepl("_TNF_", colnames(count_Filtered))
identical(ind_meta_hepg2_yes, ind_count_hepg2_yes)

#
ControlSample <- "UniL2_hipscHEPS_DMSO NO_concC0_TP8"
ind_meta_hipscHEPS_no <- grepl("hipscHEPS",meta_Filtered$SAMPLE_ID) & grepl("noTNF", meta_Filtered$SAMPLE_ID)
ind_count_hipscHEPS_no <- grepl("hipscHEPS",colnames(count_Filtered)) & grepl("noTNF", colnames(count_Filtered))
identical(ind_meta_hipscHEPS_no, ind_count_hipscHEPS_no)

# 
ControlSample <- "UniL2_hipscHEPS_DMSO YES_concC0_TP8"
ind_meta_hipscHEPS_yes <- grepl("hipscHEPS",meta_Filtered$SAMPLE_ID) & grepl("_TNF_", meta_Filtered$SAMPLE_ID)
ind_count_hipscHEPS_yes <- grepl("hipscHEPS",colnames(count_Filtered)) & grepl("_TNF_", colnames(count_Filtered))
identical(ind_meta_hepg2_yes, ind_count_hepg2_yes)

count_Filtered[1:5,1:5]
meta_Filtered[1:5, 1:5]



df_DEGs <- runDESeq2_DEGs(count_Filtered[,ind_count_hipscHEPS_yes], 
                          meta_Filtered[ind_meta_hipscHEPS_yes,], 
                          meanID_control = "UniL2_hipscHEPS_DMSO YES_concC0_TP8", 
                          meanID_column = "meanID")

df_DEGs <- unique(left_join(df_DEGs, meta_Filtered[,c("CELL_ID", "TREATMENT", "CONCENTRATION", "TIMEPOINT", "EXP_ID", "meanID")], by = "meanID"))
df_DEGs <- separate(df_DEGs, c("Probe_ID"), into = c("GeneSymbol", "ProbeNr"), sep = "_", remove = FALSE)

df_nrDEGs <- ddply(df_DEGs, .(meanID, CELL_ID, TREATMENT, TIMEPOINT, EXP_ID), summarize, 
                   NrDEGs = sum(padj < Threshold_padj , na.rm = TRUE),
                   NrDEGsFC = sum(padj < Threshold_padj & abs(log2FoldChange) > Threshold_log2FC , na.rm = TRUE))

if(Filtering == "padj"){
  df_DEGs_filtered <- subset(df_DEGs, padj < Threshold_padj)
} else {
  if(Filtering == "padj&log2FC"){
    df_DEGs_filtered <- subset(df_DEGs, padj < Threshold_padj & abs(log2FoldChange) > Threshold_log2FC)
  } else {
    if(Filtering == "log2FC"){
      df_DEGs_filtered <- subset(df_DEGs, abs(log2FoldChange) > Threshold_log2FC)
    } else {
      warning("Unrecognized input for Filtering variable. Check what has been filled in for Filtering at 'Input user' section")
    }
  }
}

df_DEGs_FCmatrix <- dcast(df_DEGs, Probe_ID ~ meanID, value.var = "log2FoldChange")
rownames(df_DEGs_FCmatrix) <- df_DEGs_FCmatrix[,1]
df_DEGs_FCmatrix <- df_DEGs_FCmatrix[,-1]
meta_Filtered_meanID <- unique(meta_Filtered[, c("CELL_ID", "TREATMENT", "CONCENTRATION", "TIMEPOINT", "EXP_ID", "meanID")])
rownames(meta_Filtered_meanID) <- meta_Filtered_meanID$meanID

## Separate text data frames for each sample for Ingenuity Pathway Analysis ##

for(i in 1:length(unique(df_DEGs$meanID))){
  df_DEGs_sub <- df_DEGs[which(df_DEGs$meanID == unique(df_DEGs$meanID)[i]), c("GeneSymbol", "log2FoldChange", "padj")]
  write.table(df_DEGs_sub, file = paste0(output_IPA, unique(df_DEGs$meanID)[i], ".txt"), sep = "\t", row.names = FALSE)
}

## Check gene expression of individual genes 

# Log2FC
df_DEGs_GenesManual <- PathwayProbesCounts(Genelist = CheckGenes, CountTable = df_DEGs, matrix = FALSE, probes = FALSE, GenelistType = "list")
pdf(paste0(outputDir_GR, "GenesManual_DESeq2Log2FC.pdf"), width = 15, height = 15)
ggplot(df_DEGs_GenesManual, aes(x = meanID, y = log2FoldChange)) + 
  geom_bar(aes(fill = padj < Threshold_padj), stat="identity", position=position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), 
                color = "black",  position=position_dodge(0.9), width = 0.5) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
        axis.title.x = element_text(size=12, vjust = 0.25),
        axis.title.y = element_text(size=12, vjust = 1),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) +
  facet_grid(Probe_ID ~ CELL_ID, scales = "free") +
  ggtitle("Log2FC") + ylab('Log2FC (vs control)') + xlab('Conditions')
dev.off()

#Log2 normalized counts

log2Norm_long_GenesManual <- PathwayProbesCounts(Genelist = CheckGenes, CountTable = log2Norm_long, matrix = FALSE, probes = FALSE, GenelistType = "list")
log2Norm_long_GenesManual_sum <- setDT(log2Norm_long_GenesManual)[,.(meanLog2Count = mean(COUNTS), sdLog2Count = sd(COUNTS)), by = .(TREATMENT, CELL_ID, TIMEPOINT, Probe_ID, EXP_ID, meanID, CONCENTRATION)]

pdf(paste0(outputDir_GR, "GenesManual_Log2count.pdf"), width = 15, height = 15)
ggplot(log2Norm_long_GenesManual_sum, aes(x = meanID, y = meanLog2Count)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black", fill = "steelblue") +
  geom_errorbar(aes(ymin = meanLog2Count, ymax = meanLog2Count + sdLog2Count), 
                color = "black",  position=position_dodge(0.9), width = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20)) +
  theme_bw() +
  theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
        axis.title.x = element_text(size=12, vjust = 0.25),
        axis.title.y = element_text(size=12, vjust = 1),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) +
  facet_grid(Probe_ID ~ CELL_ID) +
  ggtitle("Log2 normalized counts") + ylab('Log2 normalized counts') + xlab('Conditions')
dev.off()

## Subsetting of data for specific geneset
Geneset <- CheckGenes


log2Norm_Genes <- PathwayProbesCounts(Genelist = Geneset, CountTable = log2Norm, matrix = TRUE, probes = FALSE, GenelistType = "no")

df_DEGs_FCmatrix_Genes <- PathwayProbesCounts(Genelist = Geneset, CountTable = df_DEGs_FCmatrix, matrix = TRUE, probes = FALSE, GenelistType = "no")

## Principle component analysis ##

df_log2Norm_pca <- PCADF(log2Norm)
PCAplots(df_log2Norm_pca, meta_Filtered, colourVar = "TREATMENT", Clusts = FALSE, FileID = "Log2-AllSamples-AllGenes", width = 10, height = 8, sizeP = 0.2, Dataclass = "log2", a = 0.7)
PCAplots(df_log2Norm_pca, meta_Filtered, colourVar = "CELL_ID", Clusts = FALSE, FileID = "Log2-AllSamples-AllGenes", width = 10, height = 8, sizeP = 0.2, Dataclass = "log2", a = 0.7)
pca_ALL_top25_log2Norm <- PCAcontrVar(df_log2Norm_pca)
pcaCoV(df_log2Norm_pca, meta_Filtered, filename = "Log2-AllSamples-AllGenes")

df_log2Norm_Genes_pca <- PCADF(log2Norm_Genes)
PCAplots(df_log2Norm_Genes_pca, meta_Filtered, colourVar = "TREATMENT", Clusts = FALSE, FileID = "Log2-AllSamples-Geneset", width = 10, height = 8, sizeP = 0.2, Dataclass = "log2", a = 0.7)
PCAplots(df_log2Norm_Genes_pca, meta_Filtered, colourVar = "CELL_ID", Clusts = FALSE, FileID = "Log2-AllSamples-Geneset", width = 10, height = 8, sizeP = 0.2, Dataclass = "log2", a = 0.7)
pca_ALL_top25_log2Norm_Genes <- PCAcontrVar(df_log2Norm_Genes_pca)
pcaCoV(df_log2Norm_Genes_pca, meta_Filtered, filename = "Log2-AllSamples-Geneset")

df_log2FC_pca <- PCADF(na.omit(df_DEGs_FCmatrix))
PCAplots(df_log2FC_pca, meta_Filtered_meanID, colourVar = "TREATMENT", Clusts = FALSE, FileID = "Log2FC-AllSamples-AllGenes", width = 10, height = 8, sizeP = 0.2, Dataclass = "log2FC", a = 0.7)
PCAplots(df_log2FC_pca, meta_Filtered_meanID, colourVar = "CELL_ID", Clusts = FALSE, FileID = "Log2FC-AllSamples-AllGenes", width = 10, height = 8, sizeP = 0.2, Dataclass = "log2FC", a = 0.7)
pca_ALL_top25genes_log2FC <- PCAcontrVar(df_log2FC_pca)
pcaCoV(df_log2FC_pca, meta_Filtered_meanID, filename = "Log2FC-AllSamples-AllGenes")

df_log2FC_Genes_pca <- PCADF(na.omit(df_DEGs_FCmatrix_Genes))
PCAplots(df_log2FC_Genes_pca, meta_Filtered_meanID, colourVar = "TREATMENT", Clusts = FALSE, FileID = "Log2FC-AllSamples-Geneset", width = 10, height = 8, sizeP = 0.2, Dataclass = "log2FC", a = 0.7)
PCAplots(df_log2FC_Genes_pca, meta_Filtered_meanID, colourVar = "CELL_ID", Clusts = FALSE, FileID = "Log2FC-AllSamples-Geneset", width = 10, height = 8, sizeP = 0.2, Dataclass = "log2FC", a = 0.7)
pca_ALL_top25genes_log2FC_Genes <- PCAcontrVar(df_log2FC_Genes_pca)
pcaCoV(df_log2FC_Genes_pca, meta_Filtered_meanID, filename = "Log2FC-AllSamples-Geneset")

## Heatmap ##

Heatmap(Data = log2Norm, Meta = meta_Filtered, FC = FALSE, heightGene = 0.2, widthGene = 7, FileName = "Log2-AllSamples-AllGenes")
Heatmap(Data = log2Norm_Genes, Meta = meta_Filtered, FC = FALSE, heightGene = 7, widthGene = 7, FileName = "Log2-AllSamples-Geneset", GeneName = TRUE)

Heatmap(Data = df_DEGs_FCmatrix, Meta = meta_Filtered_meanID, FC = TRUE, heightGene = 0.2, widthGene = 7, FileName = "Log2FC-AllSamples-AllGenes")
Heatmap(Data = df_DEGs_FCmatrix_Genes, Meta = meta_Filtered_meanID, FC = TRUE, heightGene = 7, widthGene = 7, FileName = "Log2FC-AllSamples-Geneset", GeneName = TRUE)

## Saving data ##	

setwd(inputDir)
write.table(meta_Filtered, file = 'meta_Filtered.txt', sep = '\t')
write.table(meta_Filtered_meanID, file = 'meta_Filtered_meanID.txt', sep = '\t')
write.table(meta, file = 'meta.txt', sep = '\t')
write.table(count, file = 'count.txt', sep = '\t')
write.table(count_Filtered, file = 'count_Filtered.txt', sep = '\t')
write.table(Norm, file = 'count_Filtered_LibSizenorm.txt', sep = '\t')
write.table(log2Norm, file = 'count_Filtered_log2LibSizenorm.txt', sep = '\t')

write.table(df_DEGs, file = paste0('DEGs_vs_', ControlSample, ".txt"), sep = "\t")

write.table(df_DEGs_filtered, file = paste0('DEGsfiltered_vs_', ControlSample, ".txt"), sep = "\t")
write.table(df_nrDEGs, file = paste0('CountDEGsfiltered_vs_', ControlSample, ".txt"), sep = "\t")

write.table(pca_ALL_top25_log2Norm, file = 'pca_ALLgenes_log2Norm_top25.txt', sep = '\t')
write.table(pca_ALL_top25_log2Norm_Genes, file = 'pca_Geneset_log2Norm_top25.txt', sep = '\t')
write.table(pca_ALL_top25genes_log2FC, file = 'pca_ALLgenes_log2FC_top25.txt', sep = '\t')
write.table(pca_ALL_top25genes_log2FC_Genes, file = 'pca_Geneset_log2FC_top25.txt', sep = '\t')

save.image(paste0(inputDir, "0_BioSpyder analysis pipeline_General_DataAssembly_v1.RData"))


