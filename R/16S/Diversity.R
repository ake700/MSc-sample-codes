library(phyloseq)
library(ggplot2)
library(ggpubr)
library(vegan)
library(Rcpp)
library(dplyr)
library(Rcolorbrewer)
library(tidyr)

#Alpha diversity
sample_data(physeq_shime)$Treatment <- as(sample_data(physeq_shime)$Treatment, 'character')
sample_data(physeq_shime)$SHIME <- as(sample_data(physeq_shime)$SHIME, 'character')

#Box plot of Shannon and Simpsons by Treatment
palette <- c("#B0F2E7", "#166AD0", "#F89EE9", "#DA0000", "#C6C3D3", "#23202C")

#Shannon
sh_shime <- plot_richness(physeq_shime, x='SHIME', measures='Shannon') + 
            theme_bw() + xlab('SHIME') + scale_y_continuous(limits=c(2.0,3.2)) + 
            theme(legend.position='none') + 
            geom_boxplot(lwd=0.7, alpha=0.7, fill='light green') + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank())

#Simpsons
si_shime <- plot_richness(physeq_shime, x='SHIME', measures='Simpson') + theme_bw() + 
            xlab('SHIME') + scale_y_continuous(limits=c(0.75,1)) + 
            theme(legend.position='none') + 
            geom_boxplot(lwd=0.7, alpha=0.7, fill='light green') + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank())

#Combine all indices into one ggplot
ggarrange(sh_shime, si_shime,  
          labels=c('A','B'), 
          align='hv') 

#Optional: export data frame containing number of standard alpha diversity estimates
rich <- estimate_richness(physeq_shime, measures=c('Simpson','Shannon'))

#Evaluate significant differences in species diversity
pairwise.wilcox.test(rich$Shannon, sample_data(physeq_shime3)$SHIME, p.adjust.method = 'BH')
pairwise.wilcox.test(rich$Simpson, sample_data(physeq_shime3)$SHIME, p.adjust.method = 'BH')

#Beta diversity
wunif_dist <- phyloseq::distance(physeq_shime, method='unifrac', weighted=F)
ordination <- ordinate(physeq_shime, method="PCoA", distance=wunif_dist)
bray <- distance(physeq_shime, method='bray', type='samples')
ordinate <- ordinate(physeq_shime, "PCoA", distance=bray)

palette <- c("#B0F2E7", "#166AD0", "#F89EE9", "#DA0000", "#C6C3D3", "#23202C")

#Unweighted Unifrac
unif <- plot_ordination(physeq_shime, ordination, color="System") + theme(aspect.ratio=1) +
  geom_point(size=5) + scale_color_manual(values = palette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size=15), legend.text=element_text(size=18))
unif <- unif + ggtitle(paste("PCoA using Unweighted Unifrac",sep=''))
unif_system <- unif + stat_ellipse(geom = "polygon", type="norm", 
                                   level = 0.8, alpha=0.05, aes(fill=System)) + 
  scale_fill_manual(values=palette) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) + xlab('PC1 [24.8%]') + ylab('PC2 [17.7%]')

#Bray-Curtis Distance
bray1 <- plot_ordination(physeq_shime, ordinate, color='System') + theme(aspect.ratio=1) +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              text = element_text(size=15), legend.text=element_text(size=18))+ geom_point(size=5)+
         scale_color_manual(values = palette)
bray1 <- bray1 + ggtitle(paste("PCoA using Bray-Curtis Distance Matrix",sep='')) 
bray_system <- bray1 + stat_ellipse(geom = "polygon", type="norm", level = 0.8,
                                           alpha=0.05, aes(fill=System)) + 
                                    scale_fill_manual(values=palette) + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                                                           panel.grid.minor = element_blank(),
                                                                                           panel.background = element_blank()) +
                                    xlab('PC1 [57.4%]') + ylab('PC2 [17.3%]')

#Combine Unif and Bray figures
ggarrange(unif_system, bray_system, legend = FALSE,  
          labels=c('A','B'), 
          align='hv') 

#PERMANOVA statistics
adonis(formula = wunif_dist ~ sample_data(physeq_shime)$System, permutations = 10000)
adonis(formula = bray ~ sample_data(physeq_shime)$System, permutations = 10000)
#Compile P-values and adjust for Benjamini-Hochberg using
p.adjust(as.matrix(filename), method='BH')

#Relative abundance of genera and phyla 
S1_RA1 <- read.csv("S1_RA1.csv", header=TRUE)
S1_RA1_Phyla <- read.csv('S1_RA1_Phyla.csv', header=TRUE)

Shime1_RA1 <- S1_RA1 %>%
              group_by(Timepoint, Taxa) %>%
              summarise(n = sum(Relative.Abundance)) %>%
              mutate (percentage = n / sum(n))

Shime1_RA1_Phyla <- S1_RA1_Phyla %>%
                    group_by(Timepoint, Taxa) %>%
                    summarise(n = sum(Relative.Abundance)) %>%
                    mutate (percentage = n / sum(n))

#Set color palette
colourCount <- length(unique(Shime1_RA1$Taxa))
nb.cols <- 15
getPalette <- colorRampPalette(brewer.pal(9, 'Set1'))(nb.cols)
Phyla_palette <- c('#ca0020', '#f4a582', '#f7f7f7', '#92c5de', '#0571b0')

Shime1_RA1 <- Shime1_RA1 %>% ungroup %>%
              complete(Timepoint, Taxa, fill = list(n = 0, percentage = 0))
Shime1_RA1_Phyla <- Shime1_RA1_Phyla %>% ungroup %>%
                    complete(Timepoint, Taxa, fill = list(n = 0, percentage = 0))
  
S1_16S_Genus <- ggplot(Shime1_RA1, aes(x = factor(Timepoint), y = percentage, fill = Taxa, group = Taxa)) +
                geom_area(position = "fill", colour = "black", size = .5, alpha = .7) +
                scale_y_continuous(name="Relative Abundance", expand=c(0,0)) +
                scale_x_discrete(expand=c(0,0), 
                                 limits=c('-10', '-3', '0', "1", "2", "3",'7','9','10','11','12','16')) + 
                scale_fill_manual(values = getPalette) + 
                theme(legend.position='bottom') + theme_bw() + geom_vline(xintercept=3, linetype='dotted') +
                geom_vline(xintercept=8, linetype='solid') + theme(axis.title.x=element_blank()
                                                                   ,text=element_text(size=13.5))  +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank()) 

S1_16S_Phyla <- ggplot(Shime1_RA1_Phyla, aes(x = factor(Timepoint), y = percentage, fill = Taxa, group = Taxa)) +
                geom_area(position = "fill", colour = "black", size = .5, alpha = .7) +
                scale_y_continuous(name="Relative Abundance", expand=c(0,0)) +
                scale_x_discrete(name="Timepoint (d)", expand=c(0,0), 
                                 limits=c('-10', '-3', '0', "1", "2", "3",'7','9','10','11','12','16')) + 
                scale_fill_manual(values = Phyla_palette) + 
                theme(legend.position='bottom') + theme_bw() + geom_vline(xintercept=3, linetype='dotted') +
                geom_vline(xintercept=8, linetype='solid') + theme(axis.title.x=element_blank()
                                                                   ,text=element_text(size=13.5))  +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank()) 
