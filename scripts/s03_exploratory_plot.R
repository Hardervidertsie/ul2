


exploratory_plot <- function(inputdata){
  
  
  
  plotname <- unique(inputdata$variable)
  
  pdf(file = paste0("../output/", plotname, ".pdf"), height = 80, width = 10)
  
  p <- ggplot(data = inputdata, aes(x = as.numeric(timeID), y = as.numeric(value), color = as.numeric(doseLevel) )) +
    geom_point(aes(shape = replID)) +
    facet_grid(treatment ~ cell_line, scales = "free") + # "free_x" "free_y" 
    scale_color_gradient(low="blue", high="red")
  print(p)
    
  dev.off()
  
}



exploratory_plot_dose <- function(inputdata){
  
  if(length(unique(inputdata$cell_line)) !=1 ){
    stop("deze functie mag maar voor 1 cell lijn")
  }
 
  
  plotname <- unique(inputdata$variable)
  
  pdf(file = paste0("../output/", plotname, "dose.pdf"), height = 40, width = 30)
  
  p <- ggplot(data = inputdata, aes(x = as.numeric(timeID), y = as.numeric(value), color = plateID )) +
    geom_point(aes(shape = replID)) +
    facet_grid(treatment ~ doseLevel) #+   # "free_x" "free_y" , scales = "free"
    #scale_color_gradient(low="blue", high="red")
  
  print(p)
  
  dev.off()
  
}

