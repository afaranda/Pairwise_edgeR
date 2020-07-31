##############################################################################
#                                                                            #
#  File: Pairwise_Tests_EDASeq.R                                             #
#  Author: Adam Faranda                                                      #
#  Created: July 27, 2020                                                    #
#  Purpose:  Estimate magnitude and statistical significance of              #
#            differential expression using edgeR's exactTest and             #
#            Quasi-Likliehood methods; Use CQN to normalize for              #
#            gene length OR GC Content in lieu of the TMM method             #
#                                                                            #
##############################################################################
library(edgeR)
library(dplyr)
library(cqn)
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
  lt<-read.csv(
    "data/EDASeq_Biomart_MouseGenes.csv",row.names = 1
  )
  lt$gene_id <- row.names(lt)
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
  lt<-lt %>% filter(!is.na(length) & !is.na(gc))
  write.csv(lt, fn)
} 

# Assemble master "DGE List" Object. 
df=hc_buildDataFrame(ds, ft, return_matrix = F)
master<-DGEList( 
  counts = df[row.names(lt),],
  genes = lt,
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

# Prepare DGEList for Differential Expression Analysis
design<-model.matrix(~0+group, dge$samples)     # Define Experimental Design
colnames(design)<-gsub(
  'group', '', 
  colnames(design)
)
dge<-dge[filterByExpr(dge, design), ,keep.lib.sizes=F]  # drop low count genes
dge.cqn<-dge                 # Store a copy for CQN Analysis

# Setup standard edgeR Analysis
dge<-calcNormFactors(dge)                          # Calculate Scaling Factors
dge<-estimateDisp(dge, design, robust = T)               # Estimate Dispersion
fit<-glmQLFit(dge, design, robust = T)        # Fit QLF Model to design matrix

## Generate diagnostic plots. 
png("results/BCV_Plot.png")
plotBCV(dge)                                                 # BCV Plot
dev.off()

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

# Run Standard Analyses
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

  # Estimate Differential Expression for standard model
  deg.et<-as.data.frame(topTags(exactTest(dge, pair = contrasts[[c]]), n=Inf))
  deg.qt<-as.data.frame(topTags(glmQLFTest(fit,contrast=cn), n=Inf))
  
  fn<-paste(
    "results/TMM_",c,"_Exact_Test_DEG.csv", sep="")
  write.csv(deg.et, fn)
  
  fn<-paste(
    "results/TMM_",c,"_QLFTest_DEG.csv", sep="")
  write.csv(deg.qt, fn)
}

## Run a Bootstrapped CQN Normalization to find
## the minimum stable fold change / FDR
iter=10
df.et<-data.frame()
df.qt<-data.frame()
for(i in 1:iter){
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
    cqn.obj<-cqn(
      counts = dge.cqn$counts,
      x = dge.cqn$genes$gc,
      lengths=dge.cqn$genes$length,
      sizeFactors = colSums(dge.cqn$counts),
      verbose=T
    )
    dge.cqn$offset<-cqn.obj$glm.offset
    dge.cqn<-estimateDisp(dge.cqn, design = design, robust = T)
    fit<-glmQLFit(dge.cqn, design, robust = T)  
    
    gr1<-groups[[contrasts[[c]][1]]]
    gr2<-groups[[contrasts[[c]][2]]]
    print(paste("Contrast:",c))
    print(paste("Group 1:",paste(gr1, collapse=", ")))
    print(paste("Group 2:",paste(gr2, collapse=", ")))
    print(contrasts[[c]])
    print(cn)
    
    # Estimate Differential Expression for standard model
    deg.cqn.et<-as.data.frame(topTags(exactTest(dge.cqn, pair = contrasts[[c]]), n=Inf))
    deg.cqn.qt<-as.data.frame(topTags(glmQLFTest(fit,contrast=cn), n=Inf))
  }
  deg.cqn.et$iter=i
  deg.cqn.qt$iter=i
  df.et<-bind_rows(df.et, deg.cqn.et)
  df.qt<-bind_rows(df.qt, deg.cqn.qt)
}
x<-data.frame(
  iter=rep(1:10, each=10),
  ID=rep(paste("G",1:10,sep=""),times=10),
  Value=rep(seq(0.1,1,0.1),times=10) + rep(seq(0.01,0.1,0.01),each=10)
)
x
a<-x$Value[order(x$Value)]
for(i in a){
  y<-x%>%filter(Value > a)
  print(y %>% 
    group_by(iter) %>%
    summarize(TimesCalledSig=n())
  )
}

for(i in 1:10){
  print(
    nrow(
      df.et %>% 
        filter(abs(logFC) > 2 & FDR < 0.05 & iter==i) %>%
        filter(WT_Avg > 2 | P6_Avg > 2) %>%
        filter(abs(WT_Avg - P6_Avg) > 2)
    )
  )
}


