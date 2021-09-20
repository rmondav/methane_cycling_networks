library(phyloseq)
library(SpiecEasi, lib.loc='/proj/snic2020-16-196/private/')
library(igraph)
library(Matrix)
library(huge)

## making reduced tables for network calculations
objs<-list("MC_lake_1_anO2","MC_lake_2_anO2","MC_pond_anO2")
for (tmnt in objs) {
  psyobj<- readRDS(paste0("results/",tmnt,"_CH4cycling.rds"))
  table<-otu_table(psyobj, taxa_are_rows)
  name=paste0("results/",tmnt,"OTU_CH4cycling.csv")
  write.csv(table,file=name,quote=F); rm(table)
  psyobj_t<-filter_taxa(psyobj, function(x){sum(x >0) >=3}, prune = TRUE)
  psyobj_rdcd<-prune_taxa(taxa_sums(psyobj_t)>=0.1, psyobj_t)
  outpsyob<-paste0("results/",tmnt,"_rdcd_3S_1p_CH4cycling.rds")
  saveRDS(psyobj_rdcd, file=outpsyob)
  table<-otu_table(psyobj_rdcd, taxa_are_rows)
  name=paste0("results/",tmnt,"OTU_rdcd_3S_1p_CH4cycling.csv")
  write.csv(table,file=name,quote=F); rm(table); rm(psyobj)
}


## complete list of full and reduced tables
objs<-list("MC_lake_1_anO2","MC_lake_2_anO2","MC_pond_anO2","MC_lake_1_anO2_rdcd_3S_1p")

for (tmnt in objs) {
  ##read in phyloseq object
  psyobj<- readRDS(paste0("./results/",tmnt,"_CH4cycling.rds"))
  ## make OTU table
  table<-otu_table(psyobj, taxa_are_rows)
  
  ##SPARCC
  tablet<-t(table)
  ## run sparcc
  sparcc.obj <- sparcc(tablet)
  ## replace names of OTUs as removed during processing
  rownames(sparcc.obj$Cor) <- colnames(tablet)
  colnames(sparcc.obj$Cor) <- colnames(tablet)
  rownames(sparcc.obj$Cov) <- colnames(tablet)
  colnames(sparcc.obj$Cov) <- colnames(tablet)
  ## make correlation results accessible, convert to matrix
  sparcc.corr.m<-as.matrix(sparcc.obj$Cor,dimnames(
    rownames(sparcc.obj$Cor),colnames(sparcc.obj$Cor)))
  ## convert to igraph object with named columns
  sparcc_corr.adj<-graph.adjacency(sparcc.corr.m,mode="undirected",
                                   weighted=TRUE,diag=FALSE,add.colnames=NULL)
  ## convert to dataframe
  sparcc_corr.df<-get.data.frame(sparcc_corr.adj)
  rm(sparcc.corr); rm(sparcc.corr.m); rm(sparcc_corr.adj)
  ## change correlation column name
  colnames(sparcc_corr.df)[colnames(sparcc_corr.df)=="weight"]<-"sparcc.corr"
  ## calculate pseudo pvals using reboot method
  sparcc.p.out <- sparccboot(tablet, R =999, ncpus =4)
  ## convert results to list
  sparcc.P <- pval.sparccboot(sparcc.p.out, sided = "both")
  ## convert list to dataframe
  sparcc.P_df<-as.data.frame(sparcc.P)
  rm(sparcc.p.out); rm(sparcc.P)
  ## add pseudo pvals to corr matrix and rename
  sparcc_corr.df$sparcc.ppvals<-sparcc.P_df$pvals
  ## remove edges with (arbitrary) correlation weaker than |0.3|
  sparcc_corr.rdcd<-
    subset(sparcc_corr.df, sparcc.corr>0.3 | sparcc.corr<(-0.3))
  rm(sparcc.P_df); rm(sparcc_corr.df)
  ## remove those with high-ish probability of false positive
  sparcc_corr.clean<-subset(sparcc_corr.rdcd, sparcc.ppvals<0.05)
  ## make column with named edges
  sparcc_corr.clean$edge<-
    paste(sparcc_corr.clean$from,"to",sparcc_corr.clean$to,sep="_")
  ## reorder columns so edge is first
  sparcc_edges <- sparcc_corr.clean[,c(5,3,4,1,2)] 
  ## save file as csv
  outfile3=paste0("results/",tmnt,"_CH4cycling_sparcc_edges.csv")
  ## save file as csv
  write.csv(sparcc_edges, file=outfile3, quote=F)
  rm(sparcc_corr.rdcd); rm(sparcc_corr.clean)

  
  ## PEARSON correlation with adjustment
  P<-corr.test(tablet, method="pearson", adjust="BH", alpha=0.05, ci=F)
  ## getting pval and pearsons rho adj as list from matrix
  ## get pearson correlation rho from matrix
  Pcorr<-as.matrix(P$r)
  ## make igraph obj of correlation coeff "r"
  ## weirdly add.colnames=NULL adds the OTU names 
  Pcorr.adj=graph.adjacency(Pcorr,mode="undirected",
                            weighted=T,diag=F,add.colnames=NULL)
  ## convert to dataframe
  Pcorr_df<-get.data.frame(Pcorr.adj)
  rm(Pcorr); rm(Pcorr.adj)
  ## change column of r from "weight" to "Pcr"
  colnames(Pcorr_df)[colnames(Pcorr_df) == "weight"] <- "Pcr"
  ## get pval adj as list from matrix
  Ppval<-as.matrix(P$p)
  ## convert to igraph object
  Ppval.adj=graph.adjacency(
    Ppval,mode="undirected",weighted=T,diag=F,add.colnames=NULL)
  ## convert to dataframe
  Ppval_df<-get.data.frame(Ppval.adj)
  rm(Ppval); rm(Ppval.adj)
  ## weirdly rows with zero correlation (pval=1) are not printed 
  ## so check data if uneven row numbers and remove from pvals
  #Ppval_df<-Ppval_df[Ppval_df$weight!=1,]
  ## add pval to correlation dataframe and change column name to "P pval"
  Pcorr_df$Ppval<-Ppval_df$weight; rm(Ppval_df)
  ## subset to only those with at least weak correlations
  Pcorr_rdcd<-subset(Pcorr_df, Pcr>0.5 | Pcr<(-0.2))
  ## keep those with low probability of being false positive
  Pcorr_clean<-subset(Pcorr_rdcd, Ppval<0.01)
  rm(Pcorr_df); rm(Pcorr_rdcd)
  ## make column with named edges
  Pcorr_clean$edge<-paste(Pcorr_clean$from,"to",Pcorr_clean$to,sep="_")
  ## reorder columns so edge is first
  Pcorr_edges <- Pcorr_clean[,c(5,3,4,1,2)] 
  ## save file as csv
  outfile4=paste0("results/",tmnt,"_Pcorr_edges.csv")
  ## save file as csv
  write.csv(Pcorr_edges, quote=F, file=outfile4)
  rm(Pcorr_clean)
  
  ## SPIEC-EASI
  paramargs <- list(rep.num=999, seed=261, ncores=4)
  ## run spieceasi
  spieceasi_out<-spiec.easi(
    table, method='mb', lambda.min.ratio=0.01, nlambda=50,
    sel.criterion='stars', pulsar.select=TRUE, pulsar.params=paramargs)
  ## make file name
  outfile1=paste0("./results/",tmnt,"_CH4cycling_spieceasi.rds")
  ## save spieceasi output to file
  saveRDS(spieceasi_out, file=outfile1)
  ## convert to table
  edgevals<-symBeta(getOptBeta(spieceasi_out), mode="maxabs")
  edgevalsm<-as.matrix(edgevals)
  rownames(edgevalsm) <- rownames(table)
  colnames(edgevalsm) <- rownames(table)
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
  outfile2=paste0("results/",tmnt,"_CH4cycling_Spieceasi_edges.csv")
  write.csv(Spieceasi_edges,file=outfile2,quote=F); rm(spieceasi_out)


  
  ## UNION tables
  ## make variables for the edge files
  infile1=paste0("results/",tmnt,"_Spieceasi_edges.csv")
  infile2=paste0("results/",tmnt,"_Pcorr_edges.csv")
  infile3=paste0("results/",tmnt,"_sparcc_edges.csv")
  ## read in files, make rownames from edge list, and remove numbered column
  if1<-read.csv(infile1,header=T); rownames(if1)<-as.character(if1[,2]); if1[1]<-NULL
  if2<-read.csv(infile2,header=T); rownames(if2)<-as.character(if2[,2]); if2[1]<-NULL
  if3<-read.csv(infile3,header=T); rownames(if3)<-as.character(if3[,2]); if3[1]<-NULL
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
  ## add in column that defines pos and neg for edge mapping 
  edge_union$type<-(ifelse(edge_union$union.corr>0,"co","me"))
  ## save file
  out_union=paste0("results/",tmnt,"_edge_union.csv")
  write.csv(edge_union, file=out_union, quote=F, row.names=F)
}

