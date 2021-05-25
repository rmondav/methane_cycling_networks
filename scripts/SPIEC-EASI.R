## install package
#library(devtools)
#install_github("zdk123/SpiecEasi", lib="/proj/snic2020-16-196/private/")
library(SpiecEasi, lib.loc='/proj/snic2020-16-196/private/')
library(phyloseq)

## set processing parameters and add seed for reproducibility
paramargs <- list(rep.num=1000, seed=261, ncores=4)

## read in biom files and parse taxonomy information
MC_anO2_lake_clean_psob<-import_biom(
  "./OTU_table_anO2_lake_cultures_2000_min17S.ovr36count.biom", 
  parseFunction=parse_taxonomy_default)
## make OTU and taxon tables
MC_anO2_lake_clean_otus=otu_table(MC_anO2_lake_clean_psob, taxa_are_rows)
MC_anO2_lake_clean_taxa=tax_table(MC_anO2_lake_clean_psob)
## save tables for later use
write.csv(MC_anO2_lake_clean_otus, file = "./MC_anO2_lake_clean_otus.csv", row.names = FALSE)
write.csv(MC_anO2_lake_clean_taxa, file = "./MC_anO2_lake_clean_taxa.csv", row.names = FALSE)
## run spiec.easi 
MC_anO2_lake_spieceasi <-spiec.easi(
  MC_anO2_lake_clean_otus, method='mb', lambda.min.ratio=0.01, 
  nlambda=50, sel.criterion='stars', pulsar.select=TRUE, pulsar.params=paramargs)
saveRDS(MC_anO2_lake_spieceasi, "../results/MC_anO2_lake_spieceasi.rds")
#MC_anO2_lake_spieceasi <- readRDS("../results/MC_anO2_lake_spieceasi.rds")


MC_anO2_pond_clean_psob<-import_biom(
  "./OTU_table_anO2_pond_cultures_2000_min11S.ovr30count.biom", 
  parseFunction=parse_taxonomy_default)
MC_anO2_pond_clean_otus=otu_table(MC_anO2_pond_clean_psob, taxa_are_rows)
MC_anO2_pond_clean_taxa=tax_table(MC_anO2_pond_clean_psob)
write.csv(MC_anO2_pond_clean_otus, file = "./MC_anO2_pond_clean_otus.csv", row.names = FALSE)
write.csv(MC_anO2_pond_clean_taxa, file = "./MC_anO2_pond_clean_taxa.csv", row.names = FALSE)
## run spiec.easi
MC_anO2_pond_spieceasi <-spiec.easi(
  MC_anO2_pond_clean_otus, method='mb', lambda.min.ratio=0.05, 
  nlambda=50, sel.criterion='stars', pulsar.select=TRUE, pulsar.params=paramargs)
## and check lamda selection
#getOptInd(MC_anO2_pond_spieceasi)
#sum(getRefit(MC_anO2_pond_spieceasi))/2
#MC_anO2_pond_spieceasi$select$stars$summary
#getStability(MC_anO2_pond_spieceasi)
saveRDS(MC_anO2_pond_spieceasi, "../results/MC_anO2_pond_spieceasi.rds")

MC_O2_lake_clean_psob<-import_biom(
  "OTU_table_O2_lake_cultures_2000_min11S.ovr30count.biom", 
  parseFunction=parse_taxonomy_default)
MC_O2_lake_clean_otus=otu_table(MC_O2_lake_clean_psob, taxa_are_rows)
MC_O2_lake_clean_taxa=tax_table(MC_O2_lake_clean_psob)
write.csv(MC_O2_lake_clean_otus, file = "./MC_O2_lake_clean_otus.csv", row.names = FALSE)
write.csv(MC_O2_lake_clean_taxa, file = "./MC_O2_lake_clean_taxa.csv", row.names = FALSE)
MC_O2_lake_spieceasi <-spiec.easi(
  MC_O2_lake_clean_otus, method='mb', lambda.min.ratio=0.01, 
  nlambda=50, sel.criterion='stars', pulsar.select=TRUE, pulsar.params=paramargs)
saveRDS(MC_O2_lake_spieceasi, "../results/MC_O2_lake_spieceasi.rds")

TS_lake_2000_clean_psob<-import_biom(
  "OTU_table_TS_lake_2000_min18S.ovr38count.biom", 
  parseFunction=parse_taxonomy_default)
TS_lake_2000_clean_otus=otu_table(TS_lake_2000_clean_psob, taxa_are_rows)
TS_lake_2000_clean_taxa=tax_table(TS_lake_clean_psob)
write.csv(TS_lake_2000_clean_otus, file = "./TS_lake_2000_clean_otus.csv", row.names = FALSE)
write.csv(TS_lake_2000_clean_taxa, file = "./TS_lake_2000_clean_taxa.csv", row.names = FALSE)
TS_lake_2000_spieceasi <-spiec.easi(
  TS_lake_2000_clean_otus, method='mb', lambda.min.ratio=0.01, nlambda=50, 
  sel.criterion='stars', pulsar.select=TRUE, pulsar.params=paramargs)
saveRDS(TS_lake_2000_spieceasi, "../results/TS_lake_2000_spieceasi.rds")

