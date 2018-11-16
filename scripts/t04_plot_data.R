

head(formatted_data)



col_names <- c(paste("siNrf2", 1:24), paste("siRelA", 1:24), paste("siCtrl1", 1:24))
ind <- match(col_names, colnames(formatted_data[, -1]))

formatted_data <- as.data.frame(formatted_data)
rownames(formatted_data) <- formatted_data$row_id
formatted_data$row_id <- NULL
formatted_data <- formatted_data[,  ind ]

tmpcol <-brewer.pal(9, "PuOr")
colfunc <- colorRampPalette(c("white", tmpcol[9]))
my_color <- colfunc(100)

pie(rep(1, length(my_color)), labels = sprintf("%d (%s)", seq_along(my_color), 
                                               my_color), col = my_color)

head(formatted_data)
pdf(file = "../results_kd/heatmap_no_scaling.pdf", height = 4, width = 12)
pheatmap(formatted_data, cluster_rows = FALSE, cluster_cols = FALSE, color =  my.colors)
dev.off()

# scale per reporter

formatted_data_block_mmn <-  formatted_data
minmaxFun <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
formatted_data_block_mmn[1:4, ]  <- minmaxFun(formatted_data_block_mmn[1:4, ])
formatted_data_block_mmn[5:8, ]  <- minmaxFun(formatted_data_block_mmn[5:8, ])
formatted_data_block_mmn[9:12, ]  <- minmaxFun(formatted_data_block_mmn[9:12, ])
formatted_data_block_mmn[13:16, ]  <- minmaxFun(formatted_data_block_mmn[13:16, ])


# anotation
head(formatted_data_block_mmn)
colnames(formatted_data_block_mmn)

col_anot <- data.frame(siRNA = gsub(" [0-9]{1,2}$", "", colnames(formatted_data_block_mmn)),
                        time_index = str_match(colnames(formatted_data_block_mmn), "[0-9]{1,2}$")
)

rownames(col_anot) <- colnames(formatted_data_block_mmn)
col_anot$time_index <- as.numeric(col_anot$time_index)

row_anot <- rownames(formatted_data_block_mmn)
  
row_anot <- data.frame(cell_line = gsub(
                         " ", "", str_match(
                          row_anot, "^[A-Za-z0-2]{3,5} ")
                          ), 
                       TNF = gsub( " ", "", gsub(
                         " [1-4]{1}", "",str_match(
                         row_anot, " [DEMTNF]{3,4} [1-4]{1}$")
                         )),
                       concentration = str_match(
                         row_anot, "[1-4]{1}$"
                       )
  )
  
row_anot$concentration <- as.numeric(row_anot$concentration)
  
rownames(row_anot) <- rownames(formatted_data_block_mmn)

anot_colors = list(cell_line = c("blue", "blue", "orange", "orange"),
                   TNF = c("grey", "black"),
                   siRNA = c("blue", "black", "green")
                   
                   )
names(anot_colors$cell_line) <- c("Nrf2", "SRXN1", "A20", "ICAM1")
names(anot_colors$TNF) <- c("DMEM", "TNF")
names(anot_colors$siRNA) <- c("siNrf2", "siRelA", "siCtrl1")




pdf(file = "../results_kd/heatmap_reporter_mmn.pdf", height = 6, width = 12)
pheatmap(formatted_data_block_mmn, cluster_rows = FALSE, cluster_cols = FALSE, color =  my_color,
         annotation_row = row_anot, annotation_col = col_anot, annotation_colors = anot_colors)
dev.off()



