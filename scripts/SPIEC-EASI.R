## install package
#library(devtools)
#install_github("zdk123/SpiecEasi", lib="/proj/snic2020-16-196/private/")
library(SpiecEasi, lib.loc='/proj/snic2020-16-196/private/')
library(phyloseq)
library(igraph)

## read in biom files and parse taxonomy information
MC_anO2_lake1_clean_psob<-import_biom(
  "./OTU_table_anO2_lake1_cultures_2000_min9S.ovr28count.biom", 
  parseFunction=parse_taxonomy_default)
MC_anO2_lake2_clean_psob<-import_biom(
  "OTU_table_anO2_lake2_cultures_2000_min7S.ovr26count.biom",
  parseFunction=parse_taxonomy_default)
TS_lake_2000_clean_psob<-import_biom(
  "OTU_table_TS_lake_2000_min18S.ovr38count.biom",
  parseFunction=parse_taxonomy_default)

## set spieceasi processing parameters and add seed for reproducibility
paramargs <- list(rep.num=1000, seed=261, ncores=4)

tn = list("MC_anO2_lake1","TS_lake_2000","MC_anO2_lake2")

## set up for loop to iterate over list "tn"
for (tmnt in tn) {
  ## make file and variable names
  infile=paste0(tmnt,"_clean_psob")
  otus=paste0(tmnt,"_clean_otus.csv")
  taxs=paste0(tmnt,"_clean_taxa")
  out_OTU=paste0("processed_tables/",tmnt,"_clean_otus.csv")
  out_TAX=paste0("processed_tables/",tmnt,"_clean_taxa.csv")
  out_RDS=paste0("results/",tmnt,"_spieceasi.rds")
  ## make OTU and taxon tables
  otus<-otu_table(infile, taxa_are_rows)
  taxs<-tax_table(infile)
  ## save tables for later use
  write.csv(otus, file=out_OTU, row.names=T)
  write.csv(taxs, file=out_TAX, row.names=T)
  ## run spiec.easi 
  spieceasi<-spiec.easi(otus, method='mb', lambda.min.ratio=0.01, nlambda=50, 
    sel.criterion='stars', pulsar.select=TRUE, pulsar.params=paramargs)
  saveRDS(spieceasi, file=out_RDS)

  ## convert rds to csv
  edgevals=symBeta(getOptBeta(spieceasi), mode="maxabs")
  edgevalsm<-as.matrix(edgevals)
  #re-get OTU names as lost in matrix conversion
  rownames(edgevalsm) <- rownames(otus)
  colnames(edgevalsm) <- rownames(otus)
  ## convert matrix to adjacency graph
  edgevals.adj<-graph.adjacency(
    edgevalsm,mode="undirected",weighted=TRUE,diag=FALSE,add.colnames=NULL)
  Evals_df<-get.data.frame(edgevals.adj)
  ## change column name from "weight" to "SE weight"
  colnames(Evals_df)[colnames(Evals_df) == "weight"] <- "spieceasi.weight"
  ## make column with named edges
  Evals_df$edge<-paste(Evals_df$from,"to",Evals_df$to,sep="_")
  ## reorder columns so edge is first
  Spieceasi_edges <- Evals_df[,c(4,3,1,2)]
  ## save file as csv
  outfile=paste0("results/",tmnt,"_Spieceasi_edges.csv")
  write.csv(Spieceasi_edges, file=outfile, quote=F, row.names=F)
}


## if spieceasi doesnt find any edges can help to
## check lamda selection during testing
#getOptInd(MC_anO2_lake1_spieceasi)
#sum(getRefit(MC_anO2_lake1_spieceasi))/2
#MC_anO2_lake1_spieceasi$select$stars$summary
#getStability(MC_anO2_lake1_spieceasi)

