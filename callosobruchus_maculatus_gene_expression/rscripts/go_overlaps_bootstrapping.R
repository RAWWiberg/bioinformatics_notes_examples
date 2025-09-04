#!/usr/bin/Rscript

Custom<-ViSEAGO::Custom2GO(file = here::here("rdata","C_mac_HiFi_customfile.18-03-2025.txt"))
# Create GENE2GO
myGENE2GO<-ViSEAGO::annotate(id="Cmac",object=Custom)


# The background gene set is defined as all transcripts in the genome.
# This is despite the fact that many did not have expression data.
Transcripts <- read.table(here::here("rdata","Transcript_coords.bed"),
                          header = FALSE, sep = "\t")
colnames(Transcripts) <- c('contig','start_T','end_T','transcript_ID','gene_ID')

# Load DE results
qlf_results_tables<-readRDS(file=here::here("rdata/adults","selection_lines_de_results.RData"))
load(file = here::here("rdata","go_results.RData"))

nboots<-100
#n_genes<-1000
boot_dat<-data.frame("boot"=seq(1,nboots),
                     "SAvL3_ovlp"=vector(length=nboots),
                     "L1vL3_ovlp"=vector(length=nboots),
                     "unique_L1_ovlp"=vector(length=nboots), # The overlap of unique L
                     "unique_SA_ovlp"=vector(length=nboots))
# The number of transcripts that should be sampled for each overlap comparison
# e.g. for the overlap with SAvL3 GO terms, I need to sample the same number of transcripts
# as are DE in the L1vL3 comparison.
set_n<-c(2842, 1550, 2076, 784)
names(set_n)<-c("SAvL3","L1vL3","unique_SA","unique_L1")
i<-1
#set<-"SAvL3"
set.seed(18032025)
for(i in 1:nboots){
  Backgrund_GENES <- Transcripts$transcript_ID
  for(set in names(set_n)){
    n_genes<-set_n[set]
    # Get a random set of transcripts
    # Because each of the contrasts tested the same number of transcripts
    # I simply use the table for L1vL3 here as a source for transcript names.
    sel_GENES<-qlf_results_tables[["L1vL3"]]$trans[
      sample(1:length(qlf_results_tables[["L1vL3"]]$trans),n_genes,replace=FALSE)]
    #sel_GENES<-Transcripts$transcript_ID[
    #sample(1:length(Transcripts$transcript_ID),n_genes,replace=FALSE)]
    
    # create topGOdata
    BP<-ViSEAGO::create_topGOdata(allGenes = Backgrund_GENES, geneSel = sel_GENES,
                                  gene2GO=myGENE2GO, 
                                  ont="BP",
                                  nodeSize=10)
    
    # ViSEAGO says there are "13350 feasible genes" that can be used.
    # How many of my "strict" genes are in these feasible sets
    #str(BP)
    feasible_genes<-BP@allGenes[BP@feasible]
    n<-length(sel_GENES[which(sel_GENES %in% feasible_genes)])
    tot<-length(sel_GENES)
    cat("boot ",i," set: ",set," -- Genes w/ GO: ",n,", ",tot,", ",n/tot,"\n")
    
    # perform TopGO test using clasic algorithm
    classic<-topGO::runTest(
      BP,
      algorithm ="classic",
      statistic = "fisher"
    )
    # merge results from topGO
    BP_sResults<-ViSEAGO::merge_enrich_terms(
      Input=list(
        p.value=c("BP","classic")
      )
    )
    # What proportion of the GO terms are overlapping with GO terms in:
    boot_dat[,paste(set,"_ovlp",sep="")][i]<-length(
      which(BP_sResults@data$GO.ID %in% go_results_table[[set]]$GO.ID))/length(BP_sResults@data$GO.ID)

  }
}

write.table(boot_dat,here::here("rdata",paste("go_boot_results_SAvL3-L1vL3_ovlap2.tab",sep="")),
            col.names=TRUE,row.names=FALSE,quote = FALSE,sep="\t")
