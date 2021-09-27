library(dplyr)
## make sample list
tn <-c("TS_lake_2000","MC_anO2_lake1","MC_anO2_lake2")
## set up for loop to iterate over list "tn"
for (tmnt in tn) {
  ## make variables for the edge files
  infile1=paste0("results/",tmnt,"_Spieceasi_edges.csv")
  infile2=paste0("results/",tmnt,"_Pcorr_edges.csv")
  infile3=paste0("results/",tmnt,"_sparcc_edges.csv")
  ## read in files, make rownames from edge list, and remove numbered column
  if1<-read.csv(infile1,header=T); rownames(if1)<-as.character(if1[,1])
  if2<-read.csv(infile2,header=T); rownames(if2)<-as.character(if2[,1])
  if3<-read.csv(infile3,header=T); rownames(if3)<-as.character(if3[,1])
  ## make lists of edges for each analysis
  spcE<-rownames(if1); persE<-rownames(if2); sparccE<-rownames(if3)
  ## get list of first pair intersect
  spcE_persE<-intersect(spcE,persE)
  ## use list to subset
  spcE_df<-if1[spcE_persE,]
  persE_df<-if2[spcE_persE,]; persE_df[4:5]<-NULL #remove to and from columns
  ## add the columns together to make single dataframe
  spcE_persE_df<-cbind(spcE_df,persE_df); spcE_persE_df[1]<-NULL
  
  ## second pair intersect
  spcE_sparccE<-intersect(spcE,sparccE)
  spcE_df<-if1[spcE_sparccE,]
  sparccE_df<-if3[spcE_sparccE,]; sparccE_df[4:5]<-NULL
  spcE_sparccE_df<-cbind(spcE_df,sparccE_df); spcE_sparccE_df[1]<-NULL
  
  ## third pair intersect
  persE_sparccE<-intersect(persE,sparccE)
  persE_df<-if2[persE_sparccE,]
  sparccE_df<-if3[persE_sparccE,]; sparccE_df[4:5]<-NULL
  persE_sparccE_df<-cbind(persE_df,sparccE_df); persE_sparccE_df[1]<-NULL
  
  ## get unions of all paired intersects and outer join together 
  union.1.2_df<-full_join(spcE_persE_df,spcE_sparccE_df, keep=F, 
                          by=c("edge","from","to","spieceasi.weight"))
  union.1.2.3_df<-full_join(union.1.2_df,persE_sparccE_df, keep=F,
                            by=c("edge","from","to","Pcr","Ppval",
                                 "sparcc.corr","sparcc.ppvals"))
  ## reorder for ease of upload to cytoscape
  edge_union <- union.1.2.3_df[,c(4,2,3,1,5,6,7,8)]
  ## calculate the average of the correlation statistics (yes even though the are different!)
  edge_union$union.corr<-rowMeans(
    edge_union[,c("sparcc.corr","Pcr","spieceasi.weight")], na.rm=TRUE)
  ## save file
  out_union=paste0("results/",tmnt,"_edge_union.csv")
  write.csv(edge_union, file=out_union, quote=F, row.names=F)
}

