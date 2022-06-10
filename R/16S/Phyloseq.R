#Read in OTU tables
otu <- read.table(file="shime_otu_table.txt", header=TRUE)

#Read in taxonomy table
tax <-read.table(file='shime_taxonomy.tsv', sep='\t', header=TRUE)

#Merge files
merged_shime <- merge(otu, tax, by.x=c('OTUID'), by.y=c('OTUID'))

#Output merged .txt file
write.table(merged_shime, file='combined_otu_tax', sep='\t', col.names=TRUE, row.names=FALSE)

#Read in updated OTU table and convert to a matrix
otu_shime_table <- read.csv("shime_otu2.csv", sep=',', row.names=1)
otu_shime <- as.matrix(otu_shime_table)

#Read in taxonomy and convert to a matrix
taxonomy_shime <- read.csv('shime_tax.csv', sep=',', row.names=1)
tax_shime <- as.matrix(taxonomy_shime)

#Read in metadata and phylogenetic tree
metadata_shime <- read.csv('shime_meta2.csv', row.names=1)
phy_tree_shime <- read_tree('shime_tree.nwk')

#Import as phyloseq objects
OTU_SHIME <- otu_table(otu_shime, taxa_are_rows = TRUE)
TAX_SHIME <- tax_table(tax_shime)
META_SHIME <- sample_data(metadata_shime)

#Check if OTU names are consistent across objects
taxa_names(OTU_SHIME)
taxa_names(TAX_SHIME)
taxa_names(phy_tree_shime)

#Do the files have the same sample names?
sample_names(OTU_SHIME)
sample_names(META_SHIME)

#Merge to one phyloseq object
physeq_shime <- phyloseq(OTU_SHIME, TAX_SHIME, META_SHIME, phy_tree_shime)
