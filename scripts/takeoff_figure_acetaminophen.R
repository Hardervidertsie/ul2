



acetaminophen <- sel_data[sel_data$treatment %in% c("Acetaminophen", "DMSO", "DMEM") & sel_data$cell_line == "SRXN1",]
head(acetaminophen)

unique(acetaminophen[ acetaminophen$treatment== "Acetaminophen", "variable" ])

require(ggplot2)
require(dplyr)


acetaminophen <- acetaminophen %>% filter( !(treatment == "DMSO" & dose_uM == 0.63)) %>%
  filter(variable == "Integrated.Intensity SRXN1") %>%
  mutate(dose_uM = if_else(treatment == "DMSO",   0,  dose_uM  )) %>%
  mutate(treatment = "Acetaminophen")






dir('../')
pdf('../resultaten/figuren/takeoff.pdf', width = 12, height =50)
ggplot(data = acetaminophen, aes( x = timeID, y = value , color = replID )) + facet_wrap( ~dose_uM, ncol = 1 ) +
  geom_point() + geom_smooth() + theme_bw()

dev.off()
head(acetaminophen)
