##############################################################################
#                                                                            #
#  File: Annotate_Ensembl_Features.R                                         #
#  Author: Adam Faranda                                                      #
#  Created: June 18, 2020                                                    #
#  Purpose:  Add annotations to gene length table                            #
#                                                                            #
##############################################################################
library(edgeR)
library('AnnotationHub')
library('org.Mm.eg.db')
library(dplyr)
options(echo=T)

# Enter Working Directory and Load Raw Data
setwd('/Users/adam/Pairwise_edgeR')
wd<-getwd()
data_dir<-paste(wd,'data',sep='/')

# Import "Gene Length" Annotation to add to dgelist
lt<-read.table(
  'data/gene_coding_lengths.txt', 
  header=T, quote="", sep="\t", 
  stringsAsFactors = F
)


# Import Annotations, join on "Gene Length" Table
ah<-AnnotationHub()

# Run Query to find proper Annotation Set:
# AnnotationHub::query(ah, pattern=c("EnsDb", "Mus musculus", "98"))
edb<-ah[['AH75036']]
lt<-merge(
  lt, AnnotationDbi::select(
    edb, keys=lt$gene_id, 
    columns = c("SYMBOL", "DESCRIPTION", "GENEBIOTYPE", "SEQCOORDSYSTEM"), 
    keytype = "GENEID"
  ),
  by.x = 'gene_id', by.y='GENEID'
)

# Join Entrez ID's by gene symbol
lt<-merge(
  lt, AnnotationDbi::select(
    org.Mm.eg.db, keys=unique(lt$SYMBOL), 
    columns = c("ENTREZID"), 
    keytype = "SYMBOL"
  ), by='SYMBOL'
)
row.names(lt)<-lt$gene_id

# Write out data Gene Annotation Table:
fn<-paste(data_dir, "Gene_Annotations.csv", sep='/')
write.csv(lt, fn)




