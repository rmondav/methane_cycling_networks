library(SpiecEasi, lib.loc='/proj/snic2020-16-196/private/')
library(phyloseq)
library(igraph)
tn = list("MC_anO2_lake1","MC_anO2_pond","MC_anO2_lake2","TS_lake_2000")
## set up for loop to iterate over list "tn"
for (tmnt in tn) {
  infile=paste0("processed_tables/",tmnt,"_clean_otus.csv")
  rdsfile=paste0("results/",tmnt,"_spieceasi.rds")
  clean_otus<-read.csv(file =infile, header=T, row.names=1)
  spieceasi_out <- readRDS(rdsfile)
  edgevals=symBeta(getOptBeta(spieceasi_out), mode="maxabs")
  edgevalsm<-as.matrix(edgevals)
  rownames(edgevalsm) <- rownames(clean_otus)
  colnames(edgevalsm) <- rownames(clean_otus)
  edgevals.adj<-graph.adjacency(
    edgevalsm,mode="undirected",weighted=TRUE,diag=FALSE,add.colnames=NULL)
  Evals_df<-get.data.frame(edgevals.adj)
  ## change column of r from "weight" to "Pcr"
  colnames(Evals_df)[colnames(Evals_df) == "weight"] <- "spieceasi.weight"
  ## make column with named edges
  Evals_df$edge<-paste(Evals_df$from,"to",Evals_df$to,sep="_")
  ## reorder columns so edge is first
  Spieceasi_edges <- Evals_df[,c(4,3,1,2)] 
  ## save file as csv
  outfile=paste0("results/",tmnt,"_Spieceasi_edges.csv")
  write.csv(Spieceasi_edges,file=outfile,quote=F)
}

