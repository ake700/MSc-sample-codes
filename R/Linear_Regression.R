library(ggpubr)
library(ggplot2)
theme_set(theme_pubr())

shime_scale <- read.csv("shime_scale.csv", header=TRUE)
names(shime_scale)[1] <- "Time"
colnames(shime_scale)
colnames(shime_scale)[15] <- paste0('C. sakazakii')
head(shime_scale)

#Example with one metabolite
acetate <- ggscatter(shime_scale2, x = "Acetate", y = "C. sakazakii",
          color = "black", palette = c("#00AFBB", "#E7B800", "#FC4E07"), size = 3,
          # Points color, shape and size
          add = "reg.line",  # Add regression line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
           xlab =FALSE, ylab = 'C. sakazakii (log CFU/mL)',xlim=c(37, 50),
          ylim=c(0,6),
          # change axes labels and scale
          xticks.by=5
          # change tick marks on x/y axis
) + 
   stat_cor(
      aes(label = paste(..rr.label.., sep = "~`,`~")),
      p.accuracy = 0.05, r.accuracy = 0.01, hjust=-9, vjust=3, size=6
      ) + 
   stat_regline_equation(label.x = 42, label.y = 2)

#Combine linear regression plots as one figure
cowplot::plot_grid(acetate, ethanol + theme(axis.title.y=element_blank()),
                   acetone, glycine + theme(axis.title.y=element_blank()),
                   nrow=2, ncol=2, align='hv',vjust=1, hjust=-1.4,
                   labels='AUTO')
