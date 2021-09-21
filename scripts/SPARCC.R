library(Matrix); library(igraph)
#library(SpiecEasi) #for local analyses
library(SpiecEasi, lib.loc='/proj/snic2020-16-196/private/')

## make list for treatment names
tn<-c("TS_lake_2000","MC_anO2_lake1","MC_anO2_lake2")
## set up for loop to iterate over list "tn"
for (tmnt in tn) { 
  #read in files
  infile=paste0("processed_tables/",tmnt,"_clean_otus.csv")
  clean_otus<-read.csv(file =infile, header=T, row.names=1)
  ## transpose
  clean_otust<-t(clean_otus)
  ## run sparcc
  sparcc.obj <- sparcc(clean_otust)
  ## replace names of OTUs as removed during processing
  rownames(sparcc.obj$Cor) <- colnames(clean_otust)
  colnames(sparcc.obj$Cor) <- colnames(clean_otust)
  rownames(sparcc.obj$Cov) <- colnames(clean_otust)
  colnames(sparcc.obj$Cov) <- colnames(clean_otust)
  ## make correlation results accessible, convert to matrix
  sparcc.corr.m<-as.matrix(sparcc.obj$Cor,dimnames(
    rownames(sparcc.obj$Cor),colnames(sparcc.obj$Cor)))
  ## convert to igraph object with named columns
  sparcc_corr.adj<-graph.adjacency(sparcc.corr.m,mode="undirected",
                  weighted=TRUE,diag=FALSE,add.colnames=NULL)
  ## convert to dataframe
  sparcc_corr.df<-get.data.frame(sparcc_corr.adj)
  ## change correlation column name
  colnames(sparcc_corr.df)[colnames(sparcc_corr.df)=="weight"]<-"sparcc.corr"
  ## calculate pseudo pvals using reboot method
  sparcc.p.out <- sparccboot(clean_otust, R = 1000, ncpus = 4)
  ## convert results to list
  sparcc.P <- pval.sparccboot(sparcc.p.out, sided = "both")
  ## convert list to dataframe
  sparcc.P_df<-as.data.frame(sparcc.P)
  ## add pseudo pvals to corr matrix and rename
  sparcc_corr.df$sparcc.ppvals<-sparcc.P_df$pvals
  ## remove edges with (arbitrary) correlation weaker than |0.3|
  sparcc_corr.rdcd<-
    subset(sparcc_corr.df, sparcc.corr>0.3 | sparcc.corr<(-0.3))
  ## remove those with high-ish probability of being false positives
  sparcc_corr.clean<-subset(sparcc_corr.rdcd, sparcc.ppvals<0.05)
  ## make column with named edges
  sparcc_corr.clean$edge<-
    paste(sparcc_corr.clean$from,"to",sparcc_corr.clean$to,sep="_")
  ## reorder columns so edge is first
  sparcc_edges <- sparcc_corr.clean[,c(5,3,4,1,2)] 
  ## save file as csv
  outfile=paste0("results/",tmnt,"_sparcc_edges.csv")
  ## save file as csv
  write.csv(sparcc_edges, file=outfile, quote=F, row.names=F)
}

