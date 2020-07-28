##############################################################################
#                                                                            #
#  File: Pairwise_Tests.R                                                    #
#  Author: Adam Faranda                                                      #
#  Created: June 18, 2020                                                    #
#  Purpose:  Estimate magnitude and statistical significance of              #
#            differential expression using edgeR's exactTest and             #
#            Quasi-Likliehood methods.                                       #
#                                                                            #
##############################################################################
library(edgeR)
library(dplyr)
options(echo=T)

# Enter Working Directory and Load Raw Data
setwd('/Users/adam/Pairwise_edgeR')
source('scripts/BuildHTSeq_DataMatrix.R')
source('scripts/Excel_Write_Functions.R')
source('scripts/Overlap_Comparison_Functions.R')
wd<-getwd()
results<-paste(wd,'results',sep='/')
data_dir<-paste(wd,'data',sep='/')

# Generate table of data files / sample information
# After generating the file table, update it with 
# Experimental Grouping information in excel/oocalc/emacs
ft<-hc_getFileTable(
  wd = data_dir,
  dirList=c(data_dir),
  filename = "HTSeq_GeneCountFiles.csv"
)

# Assemble Data Sets
ds<-hc_loadFiles(ft)
ft<-hc_identifierConsistency(ds, ft, idCol=1)
row.names(ft)<-ft$sample



# Import Annotations, append to "Gene Length" Table
fn<-paste(data_dir,'Gene_Annotations.csv', sep='/')
if(file.exists(fn)){
  lt<-read.csv(fn, row.names = 1)
  print(TRUE)
} else {
  library('AnnotationHub')
  lt<-read.table(
    'data/gene_coding_lengths.txt', 
    header=T, quote="", sep="\t", 
    stringsAsFactors = F
  )
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
  
  lt<-merge(
    lt, AnnotationDbi::select(
      org.Mm.eg.db, keys=unique(lt$SYMBOL), 
      columns = c("ENTREZID"), 
      keytype = "SYMBOL"
    ), by='SYMBOL'
  )
  row.names(lt)<-lt$gene_id
  rm(ah, edb)
  detach(package:AnnotationHub, unload=T)
  detach(package:ensembldb, unload=T)
  detach(package:AnnotationFilter, unload=T)
  write.csv(lt, fn)
} 

# Assemble master "DGE List" Object. 
df=hc_buildDataFrame(ds, ft, return_matrix = F)
master<-DGEList( 
  counts = df,
  genes = lt[row.names(df),],
  samples = ft[colnames(df),],
)
  
# Define Groups As a list of the form: "Group1=c('S1', 'S2', 'S3')
# Each list element is the name of the Group and points to a
# Vector of sample labels, drop samples without a group. 
# A sample can only be assigned to one group.
f<-function(z, gr){
  (names(gr)[sapply(gr, function(x, y=z) y %in% x)])[1]
}
groups=list(
  WT=paste('WT',1:3,sep=''),
  P6=paste('P6',1:3,sep='')
)
dge<-master[,unlist(groups)]
dge$samples$group<-factor(
  sapply(dge$samples$sample, f, gr=groups),
  levels = names(groups)
)

# Add Average FPKMs to gene table in DGEList
for(g in names(groups)){
  cn<-paste(g, "Avg", sep="_")
  s<-groups[[g]]
  dge$genes[cn] <- apply(rpkm(dge)[,s], 1, mean)   
}

# Define A list of named contrasts; each element points to a vector with
# a pair of group labels. Positive fold changes will be associated
# with the second group listed. 
contrasts=list(
  Pax6cKO_vs_Wildtype=c('WT', 'P6')
)

# Function Calculates fold change from logFC
calcFC<-function(lfc){
	if(!is.na(lfc) & lfc > 0){
		return(2^lfc)
	}
	else if(!is.na(lfc) & lfc < 0){
	     return(0 - 2^abs(lfc))
	}
	else { return(0)}
}

# Prepare DGEList for Differential Expression Analysis
design<-model.matrix(~0+group, dge$samples)     # Define Experimental Design
colnames(design)<-gsub(
  'group', '', 
  colnames(design)
)
dge<-dge[filterByExpr(dge, design), ,keep.lib.sizes=F]  # drop low count genes
dge<-calcNormFactors(dge)                          # Calculate Scaling Factors
dge<-estimateDisp(dge, design, robust = T)               # Estimate Dispersion
fit<-glmQLFit(dge, design, robust = T)        # Fit QLF Model to design matrix

## Generate diagnostic plots. 
png("results/BCV_Plot.png")
plotBCV(dge)                                                 # BCV Plot
dev.off()

# png(
#   "results/PCA_Plot.png",   # PCA Plot
#   width=600, height = 400
# )
# plotPrinComp(
#   cpm(dge, log=T), ft=dge$samples,                  
#   idCol = 0, groupCol = 'group'
# )
# dev.off()

png(
  "results/MDS_Plot.png",  # MDS Plot
  width=600, height=500
)
par(mar=c(8.5,5,4.1,2.1))
plotMDS(
  dge, cex=1.5, cex.lab = 1.5,
  cex.axis = 1.5, cex.main = 2,
  pch = rep(15:18, each=3),
  main = "Experiment Title",
  col = rep(c('red', 'blue','black'), each=3)
)
legend(
  -3, -2.4, 
  legend = unique(dge$samples$group), 
  pch = 15:17, cex = 1.5,
  col=c('red', 'blue','black'), xpd=T
)
dev.off()

# Calculate Differential Expression based on 
for( c in names(contrasts)){
  cn<-sapply(
    colnames(design), 
    function(n){
      ifelse(
        n == contrasts[[c]][1], -1, 
        ifelse(n == contrasts[[c]][2], 1, 0
        )
      )
    }
  )
  
  gr1<-groups[[contrasts[[c]][1]]]
  gr2<-groups[[contrasts[[c]][2]]]
  print(paste("Contrast:",c))
  print(paste("Group 1:",paste(gr1, collapse=", ")))
  print(paste("Group 2:",paste(gr2, collapse=", ")))
  print(contrasts[[c]])
  print(cn)

  # Estimate Differential Expression

  deg.et<-as.data.frame(topTags(exactTest(dge, pair = contrasts[[c]]), n=Inf))
  deg.qt<-as.data.frame(topTags(glmQLFTest(fit,contrast=cn), n=Inf))
  
  fn<-paste(
    "results/Experiment_",c,"_Exact_Test_DEG.csv", sep="")
  write.csv(deg.et, fn)
  
  fn<-paste(
    "results/Experiment_",c,"_QLFTest_DEG.csv", sep="")
  write.csv(deg.qt, fn)
  
}
