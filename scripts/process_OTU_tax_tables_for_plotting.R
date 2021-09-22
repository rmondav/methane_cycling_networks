## separate OTUtax tables and add in abundance info for network plotting
## make sample list
tn <-c("TS_lake_2000","MC_anO2_lake1","MC_anO2_lake2")
## set up for loop to iterate over list "tn"
for (tmnt in tn) {
  #read in files
  infile=paste0("processed_tables/",tmnt,"_clean_otus.csv")
  clean_otus<-read.csv(file =infile, header=T, row.names=1)
  ## make new df
  names<-rownames(clean_otus)
  ## make new df
  otu_df<-as.data.frame(names)
  ## get sum of r.a. abundances for each
  otu_df$sum.r.a<-rowSums(clean_otus)
  ## get prevalence for each
  otu_df$no.samples<-rowSums(clean_otus!=0)
  ## calculate mean r.a. abundance for each
  otu_df$av.r.a.abund<-(otu_df$sum.r.a*100/2000)/(otu_df$no.samples)
  ## calculate percentage of samples for each
  otu_df$pcnt.samples<-otu_df$no.samples*100/(ncol(clean_otus))
  ## save tables 
  out_OTU=paste0("results/",tmnt,"_OTU_table_ntwk.csv")
  write.csv(otu_df, file=out_OTU, quote=F, row.names=F)
}  


## reduce OTU table list that holds taxonomy and metabolic information
## to make it easier to work with in making figures etc
## read in full OTU table
full_OTU<-read.csv("processed_tables/full_taxonomy_table.csv", header=T, row.names = 1)
## get list of OTUs in the clean OTU tables
infile="processed_tables/MC_anO2_lake1_clean_otus.csv"
clean_otus<-read.csv(file =infile, header=T, row.names=1) 
names1<-rownames(clean_otus)

infile="processed_tables/MC_anO2_lake2_clean_otus.csv"
clean_otus<-read.csv(file =infile, header=T, row.names=1) 
names2<-rownames(clean_otus)

infile="processed_tables/TS_lake_2000_clean_otus.csv"
clean_otus<-read.csv(file =infile, header=T, row.names=1) 
names3<-rownames(clean_otus)

## make list
names123<-c(names1,names2,names3)
names123U<-unique(names123)
## subset full OTU table using list
ntwk_OTU_df<-full_OTU[names123U,]
## save table
write.csv(ntwk_OTU_df, file="results/ntwk_OTU_tax_metab.csv")

