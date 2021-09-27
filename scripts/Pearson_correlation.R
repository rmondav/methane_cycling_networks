library(psych)
library(igraph)

#make list for treatment names
tn <-c("TS_lake_2000","MC_anO2_lake1","MC_anO2_lake2")
## set up for loop to iterate over list "tn"
for (tmnt in tn) {
  #read in files
  infile=paste0("processed_tables/",tmnt,"_clean_otus.csv")
  clean_otus<-read.csv(file =infile, header=T, row.names=1)
  clean_otust<-t(clean_otus)
  ## pearson correlation with adjustment
  P<-corr.test(clean_otust, method="pearson", adjust="BH", alpha=0.05, ci=F)
  ## clean
  rm(clean_otus);rm(clean_otust)

  ## getting pval and pearsons rho adj as list from matrix
  ## get pearson correlation rho from matrix
  Pcorr<-as.matrix(P$r)
  ## make igraph obj of correlation coeff "r"
  ## weirdly add.colnames=NULL adds the OTU names 
  Pcorr.adj=graph.adjacency(Pcorr,mode="undirected",
                            weighted=T,diag=F,add.colnames=NULL)
  ## convert to dataframe
  Pcorr_df<-get.data.frame(Pcorr.adj)
  ## change column of r from "weight" to "Pcr"
  colnames(Pcorr_df)[colnames(Pcorr_df) == "weight"] <- "Pcr"

  ## get pval adj as list from matrix
  Ppval<-as.matrix(P$p)
  ## convert to igraph object
  Ppval.adj=graph.adjacency(Ppval,mode="undirected",
                            weighted=T,diag=F,add.colnames=NULL)
  ## convert to dataframe
  Ppval_df<-get.data.frame(Ppval.adj)
  ## clean
  rm(P); rm(Pcorr); rm(Pcorr.adj);rm(Ppval);rm(Ppval.adj)

  ## add pval to correlation dataframe and change column name to "P pval"
  Pcorr_df$Ppval<-Ppval_df$weight
  ## subset to only those with at least weak correlations
  Pcorr_rdcd<-subset(Pcorr_df, Pcr>0.5 | Pcr<(-0.2))
  ## keep those with low probability of being false positive
  Pcorr_clean<-subset(Pcorr_rdcd, Ppval<0.01)
  ## clean
  rm(Pcorr_df);rm(Ppval_df);rm(Pcorr_rdcd)
  ## make column with named edges
  Pcorr_clean$edge<-paste(Pcorr_clean$from,"to",Pcorr_clean$to,sep="_")
  ## reorder columns so edge is first
  Pcorr_edges <- Pcorr_clean[,c(5,3,4,1,2)] 
  ## save file as csv
  fname=paste0(tmnt,"_Pcorr_edges")
  outfile=paste0("results/",fname,".csv")
  #assign(fname,Pcorr_edges)
  ## save file as csv
  write.csv(Pcorr_edges, quote=F, file=outfile)
  ## clean
  rm(Pcorr_clean); rm(Pcorr_edges)
}

