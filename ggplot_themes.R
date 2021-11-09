library(ggplot2)

palette5.customized <- c("#615da3", "#5bc1da", "#f88600", "#a1c436", "#0383bc") 
palette8 <- c("#424242", "#FF0000", "#43a047", "#0033CC", "#4fc3f7", "#f8bbd0", "#ff8f00", "#ffeb3b")
linetypes <- c("solid", "longdash", "dotted", "dotdash", "dashed", "twodash")

theme_CECLMortgageStudy <- function(gg.size.scale = 1){
  # Provides custom theme which was used for JHA reports at first, then for CECL Mortgage book written by J.Breeden
  #
  # gg.size.scale: Scale for all text in pdf files
  
  # Size in ggplot is given in mm and thus should be devided on the following "magic number" 
  mn <- ggplot2:::.pt
  
  theme(axis.title.x = element_text(family = "sans", colour = "#333333", size = 12*gg.size.scale, face = "bold",
                                    margin = margin(15*gg.size.scale,0,0,0)),
        axis.title.y = element_text(family = "sans", colour = "#333333", size = 12*gg.size.scale, face = "bold",
                                    margin = margin(0,15*gg.size.scale,0,0)),
        axis.text.x = element_text(family = "sans", colour = "#333333", size = 12*gg.size.scale),
        axis.text.y = element_text(family = "sans", colour = "#333333", size = 12*gg.size.scale),
        legend.title = element_text(family = "sans", colour = "#333333", size = 12*gg.size.scale),
        legend.text = element_text(family = "sans", colour = "#333333", size = 12*gg.size.scale),
        plot.title = element_text(family = "sans", colour = "#333333", face = "bold", size = 14, hjust = 0.5),
        plot.margin = unit(c(0.5, 1, 0.2, 0.5)*gg.size.scale, "cm"),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", size = 0.5, color = "#333333"),
        panel.grid.major.y = element_line(colour = "#d6d6d6", size = 1/mn),  
        panel.grid.minor.y = element_line(colour = "#d6d6d6", size = 1/mn),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.key = element_rect(fill = "white")) 
  
}