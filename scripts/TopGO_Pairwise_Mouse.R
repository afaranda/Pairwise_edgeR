##############################################################################
# File: TopGO_Pairwise_Mouse.R                                               #
# Purpose: Identifiy enriched ontology terms in a pairwise contrast between  #
#          samples from mouse experiments                                    #     
# Created: Feb 26, 2020                                                      #
# Author: Adam Faranda                                                       #
#                                                                            #
##############################################################################
library(topGO)
library(goseq)
library(org.Mm.eg.db)

fn<-"~/Documents/StrainContrasts/results/LF_Black6_vs_NMRI_DEG_Table.tsv"
idc='gene_id'
fdr='FDR'
lfc='logFC'

# Import data and build topGO object #########################################
deg<-read.table(fn, sep='\t', quote="", header=T, stringsAsFactors = F)
ag<-ifelse(deg[,fdr] < 0.05 & abs(deg[,lfc]) > 1, 1, 0)
names(ag)<-deg[,idc]
GOdata<-new(
  "topGOdata", ontology="BP", allGenes = ag,
  geneSelectionFun = function(allScore){
    return(allScore > 0)
  }, annot=annFUN.org,
  mapping="org.Mm.eg.db",
  ID="ensembl"
)
# Run Enrichment Tests #######################################################
res.F.C<-runTest(GOdata, algorithm = "classic", statistic = "fisher")
res.F.E<-runTest(GOdata, algorithm = "elim", statistic = "fisher")
res.F.W<-runTest(GOdata, algorithm = "weight", statistic = "fisher")

allRes<-GenTable(
  GOdata, classFisher=res.F.C, 
  elimFisher=res.F.E, weightFisher=res.F.W,
  orderBy="elimFisher", ranksOf="classFisher",
  topNodes=10, numChar=100
)





