#Cluster analysis of metabolites on a heatmap
library(NbClust)
library(factoextra)
library(circlize)
library(ComplexHeatmap)

#Read in metabolite data
shime1b_pca <- read.csv("shime1_pca_meta.csv", header=TRUE)
names(shime1b_pca)[1] <- "Time"
names(shime1b_pca)[10] <- "Propylene glycol"
shime1bmat_pca <- as.matrix(shime1b_pca[,2:14], 'numeric')
rownames(shime1bmat_pca) <- paste0(shime1b_pca$Time)
g <- shime1bmat_pca
g2 <- scale(g, scale=TRUE)
g3 <- t(g2)

#To determine the best number of clusters
res.nbclust2 <- NbClust(g, distance = 'euclidean',
          min.nc = 2, max.nc = 10, 
          method = "complete", index ="all")
fviz_nbclust(res.nbclust, ggtheme = theme_minimal())

#Set up heat map
col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
lgd <- heatmap_legend_param=list(title='Concentration (mm)', 
    at = c(-3, -2, -1, 0, 1, 2, 3))

g_heat <- Heatmap(g3, cluster_columns=FALSE, row_km=3, 
                   column_title = 'Time (d)', column_title_side='bottom', 
                   column_names_gp=grid::gpar(fontsize=20),
                   row_names_gp=grid::gpar(fontsize=20),
                   heatmap_legend_param=list(title='Normalized change in concentration', 
                                             direction='horizontal',
                                             at = c(-3, -2, -1, 0, 1, 2, 3)))
S1_cluster <- draw(g_heat, heatmap_legend_side='top')

#Principal Component Analysis
g_scale <- prcomp(g, scale=TRUE, center=TRUE)
clusters <- c(rep("SCFA", times=3), rep("Organic acid", times=3), 
            rep("Alcohols", times=3), rep('Amino acid', times=2), rep('Ketone', times=1), rep('Lipid', times=1))
group <- c(rep("S1T1: C. sakazakii", times=8), rep('S1T2: MRS media', times=6))

pca_S1 <- fviz_pca_biplot(g_scale, repel=TRUE, pointsize=5, pointshape=21, col.var=factor(clusters),
          arrowsize=0.6, labelsize=4,geom='point', fill.ind=group, invisible='quali',
          legend.title=list(fill='Treatment',color='Categories')) +
          ggpubr::fill_palette('jco') + ggpubr::color_palette("lancet") 
