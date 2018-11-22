


plot_data <- function(input_data, pdf_name, scale_reporters = FALSE) {

col_names <- c(paste("siNrf2", 1:24), paste("siRelA", 1:24), paste("siCtrl1", 1:24))
ind <- match(col_names, colnames(input_data[, -1]))

input_data <- as.data.frame(input_data)
rownames(input_data) <- input_data$row_id
input_data$row_id <- NULL

input_data <- input_data[ , ind ]

tmpcol <-brewer.pal(9, "PuOr")
colfunc <- colorRampPalette(c("white", tmpcol[9]))
my_color <- colfunc(100)

#pie(rep(1, length(my_color)), labels = sprintf("%d (%s)", seq_along(my_color), 
#                                               my_color), col = my_color)

# head(input_data)
# pdf(file = "../results_kd/heatmap_no_scaling.pdf", height = 4, width = 12)
# pheatmap(input_data, cluster_rows = FALSE, cluster_cols = FALSE, color =  my.colors)
# dev.off()

# scale per reporter

if(scale_reporters) {

minmaxFun <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
input_data[1:4, ]  <- minmaxFun(input_data[1:4, ])
input_data[5:8, ]  <- minmaxFun(input_data[5:8, ])
input_data[9:12, ]  <- minmaxFun(input_data[9:12, ])
input_data[13:16, ]  <- minmaxFun(input_data[13:16, ])
}

# anotation
head(input_data)
colnames(input_data)

col_anot <- data.frame(siRNA = gsub(" [0-9]{1,2}$", "", colnames(input_data)),
                        time_index = paste0("t", str_match(colnames(input_data), "[0-9]{1,2}$"))
)

rownames(col_anot) <- colnames(input_data)


row_anot <- rownames(input_data)
  
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
  
rownames(row_anot) <- rownames(input_data)


colfunc2 <- colorRampPalette(c("orange", "red"))
time_colors <- colfunc2(24)

anot_colors = list(cell_line = c("#69359c", "#fc8eac", "Blue", "Sky Blue"),
                   TNF = c("grey", "black"),
                   siRNA = c("#fff44f", "Tan", "#826644"),
                             time_index = time_colors
                   
                   )
names(anot_colors$cell_line) <- c("Nrf2", "SRXN1", "A20", "ICAM1")
names(anot_colors$TNF) <- c("DMEM", "TNF")
names(anot_colors$siRNA) <- c("siNrf2", "siRelA", "siCtrl1")
names(anot_colors$time_index) <- paste0("t", 1:24)

# Heatmap reporter MMN
# Parameter:	Color:
#   Time index	gradient met orange red
# siNRF2		#69359c
# siRelA		Tan
# siCtrl1		#826644
# 
# Nrf2		#69359c
# SRXN1		#fc8eac
# A20		Blue
# ICAM1 		Sky Blue




pdf(file = paste0("../results_kd/", pdf_name,".pdf"), height = 6, width = 12)
pheatmap(input_data, cluster_rows = FALSE, cluster_cols = FALSE, color =  my_color,
         annotation_row = row_anot, annotation_col = col_anot, annotation_colors = anot_colors)
dev.off()

}

