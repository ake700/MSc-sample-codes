library(GGally)
library(Hmisc)
library (ltm) #for function rcor.test
library(ComplexHeatmap)
library(ggplot2)

#Example code for one SHIME compartment

#Read in dataset and convert to matrix
S1_spear <- read.csv("S1_Spear_allb.csv", header=TRUE)
names(S1_spear)[1] <- "Time"
S1x_spear <- as.matrix(S1_spear[,2:16], 'numeric')
names(S1x_spear)[10] <- "Propylene glycol"

#Conduct Spearman's correlation analysis
S1x_spear_corr <- rcorr(S1x_spear, type="spearman")
df.S1x_spear_r <- data.frame(S1x_spear_corr$r) #generates the rho and other statistical values
df.S1x_spear_n <- data.frame(S1x_spear_corr$n)
df.S1x_spear_p <- data.frame(S1x_spear_corr$P)
#Optional adjust p-values by Benjamini-Hochberg before converting into data frame/CSV
s1_adj <- p.adjust(as.matrix(df.S1x_spear_p), method = 'BH')

#Optional, create CSVs
write.csv(df.S1x_spear_r, 'S1_spear_r')
write.csv(df.S1x_spear_n, 'S1_spear_n')
write.csv(df.S1x_spear_p, 'S1_spear_p') #or using s1_adj 

#read in CSV again if needed
shime1_spear <- read.csv("S1_spear_r.csv", header=TRUE)
names(shime1_spear)[1] <- "All"

#Heatmap preparation
shime1x_spear <- as.matrix(shime1_spear[,2:16], 'numeric')
rownames(shime1x_spear) <- paste0(shime1_spear$All)
rownames(shime1x_spear)
head(shime1x_spear)

#Generate heatmap from Spearman's values
#read in p-values
s1_p <- read.csv('S1_spear_p.csv', header=TRUE)
s1_spear <- Heatmap(shime1x_spear, circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1), 
                                                                        c('#d7191c', '#fdae61',
                                                                          '#FFFFFF', '#abd9e9', 
                                                                          '#2c7bb6')),
                                column_names_gp=grid::gpar(fontsize=20,
                                                           col=c(rep('#4daf4a', 8), rep('#984ea3',7))),
                                row_names_gp=grid::gpar(fontsize=20, 
                                                        col=c(rep('#4daf4a', 8), rep('#984ea3',7))),
                                heatmap_legend_param=list(title="Spearman's rho", 
                                                          direction='horizontal',
                                                          at = c(-1, -0.5, 0, 0.5, 1)),
                   cell_fun = function(j, i, x, y, w, h, fill){
                     if(s1_p[i, j] < 0.05 & s1_p[i,j] > 0.01) {
                       grid.text('*', x, y, gp = gpar(fontsize=20))
                       } else {
                         if(s1_p[i, j] < 0.01 & s1_p[i,j] > 0.001){
                       grid.text('**', x, y, gp = gpar(fontsize=20))}
                         else {
                           if(s1_p[i,j] < 0.001){
                             grid.text('***', x, y, gp = gpar(fontsize=20))}
                             }
                     }
                   })
s1_corr <- draw(s1_spear, heatmap_legend_side='top')
