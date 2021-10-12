#load libraries
library(ggplot2)
library(vegan)
library(phyloseq)

otulakegenera <- read.csv(file="TimeSeries_lake_top10_genera_OTUs.csv", row.names = 1)
taxlakegenera <- read.csv(file="TimeSeries_lake_top10_genera_taxa.csv", row.names = 1)
otu1mat<-as.matrix(otulakegenera)
tax1mat<-as.matrix(taxlakegenera)
OTU1 = otu_table(otu1mat, taxa_are_rows = TRUE)
TAX1 = tax_table(tax1mat)
lakegenera = phyloseq(OTU1, TAX1)
lakegenera
mapfile = "Sample_metadata.csv"
# Import sample metadata
map <- read.csv(mapfile)
map <- sample_data(map)
# Assign rownames to be Sample ID's
rownames(map) <- map$Sample.name

LakeGen <- merge_phyloseq(lakegenera, map)
LakeGen

theme_set(theme_bw())
#dev.new()
plot_bar(LakeGen, "time", fill="Species", facet_grid=~depth)
ggsave("lake.pdf")
#------------------------------
otumcgenera <- read.csv(file="MixedCultures_lake_anO2_top10_genera_OTUs.csv", row.names = 1)
taxmcgenera <- read.csv(file="MixedCultures_lake_anO2_top10_genera_taxa.csv", row.names = 1)
otu2mat<-as.matrix(otumcgenera)
tax2mat<-as.matrix(taxmcgenera)
OTU2 = otu_table(otu2mat, taxa_are_rows = TRUE)
TAX2 = tax_table(tax2mat)
mcgenera = phyloseq(OTU2, TAX2)
mcgenera
mapfile = "Sample_metadata.csv"
# Import sample metadata
map <- read.csv(mapfile)
map <- sample_data(map)
# Assign rownames to be Sample ID's
rownames(map) <- map$Sample.name

MCGen <- merge_phyloseq(mcgenera, map)
MCGen

theme_set(theme_bw())
#dev.new()
plot_bar(MCGen,"Culture_name", fill="Species")
ggsave("modelCommunities.pdf")

