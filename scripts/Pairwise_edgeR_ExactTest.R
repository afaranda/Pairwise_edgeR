library(edgeR)
library('AnnotationHub')
library(dplyr)
source('scripts/BuildDataMatrix.R')
source('scripts/Excel_Write_Functions.R')
source('scripts/Overlap_Comparison_Functions.R')
options(echo=T)

# Enter Working Directory and Load Raw Data
# setwd('/Users/afaranda/Desktop/LEC_Time_Series')
wd<-getwd()
results<-paste(wd,'results',sep='/')
data_dir<-paste(wd,'data',sep='/')

# Generate table of data files / sample information
# After generating the file table, update it with 
# Experimental Grouping information in excel/oocalc/emacs
ft<-hc_getFileTable(
  wd = data_dir,
  dirList=c(data_dir),
  filename = "HtSeq_GeneCountFiles.csv"
)

# Assemble Data Sets
ds<-hc_loadFiles(ft)
ft<-hc_identifierConsistency(ds, ft, idCol=1)

# Import "Gene Length" Annotation to add to dgelist
gt<-read.table(
  'data/gene_coding_lengths.txt', 
  header=T, quote="", sep="\t", 
  stringsAsFactors = F
)


# Import Annotations, append to "Gene Length" Table
ah<-AnnotationHub()
# Run Query to find proper Annotation Set:
# AnnotationHub::query(ah, pattern=c("EnsDb", "Mus musculus", "98"))
edb<-ah[['AH75036']]
gt<-merge(
  gt, AnnotationDbi::select(
    edb, keys=gt$gene_id, 
    columns = c("SYMBOL", "DESCRIPTION"), 
    keytype = "GENEID"
  ),
  by.x = 'gene_id', by.y='GENEID'
)
row.names(gt)<-gt$gene_id
rm(ah, edb)
detach(package:AnnotationHub, unload=T)
detach(package:ensembldb, unload=T)
detach(package:AnnotationFilter, unload=T)

# Assemble master "DGE List" Object. 
df=hc_buildDataFrame(ds, ft, return_matrix = F)
master<-DGEList( 
  counts = df,
  genes = gt[row.names(df),],
  samples = ft,
)

# Define Groups As a list of the form: "Group1=c('S1', 'S2', 'S3')
# Each list element is the name of the Group and points to a
# Vector of sample labels

groups=list(
  Control=c('1_S1', '2_S2', '4_S3'),
  Treatment=c('6_S4', '7_S5', '8_S6')
)

# A list of named contrasts; each element points to a vector with
# a pair of group labels. Positive fold changes will be associated
# with the second group listed. 
contrasts=list(
  Treatment_vs_Control=c('Control', 'Treatment')
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

# Calculate Differential Expression based on 
for( c in names(contrasts)){
  
  # dg1 = ag_master %>% filter(
  #   Experiment == pairwise[[p]]['Experiment'] &
  #     Group_1 == pairwise[[p]]['Group_1'] &
  #     Group_2 == pairwise[[p]]['Group_2']
  # )
  gr1<-groups[[contrasts[[c]][1]]]
  gr2<-groups[[contrasts[[c]][2]]]
  print(paste("Contrast:",c))
  print(paste("Group 1:",paste(gr1, collapse=", ")))
  print(paste("Group 2:",paste(gr2, collapse=", ")))
  
  # Define Grouping factor
  grp<-factor(
    c(
      rep(1, length(gr1)), 
      rep(2, length(gr2))
    ), 
    labels=c(contrasts[[c]][1], contrasts[[c]][2])
  )
  
  # Get subset of samples for this contrast
  dge<-DGEList(
    master[,c(gr1, gr2)],
    genes = gt[row.names(df),],
    group = grp
  )
  # Filter DGEList: Remove genes where fewer than two samples have a cpm > 2
  # Reccommended by EdgeR manual
  keep <- rowSums(cpm(dge) > 1) >= 2
  dge.filter <-dge[keep, ]
  dge.filter$samples$lib.size <- colSums(dge.filter$counts) # Fix library size after filtering
  
  # Calculate Normalization factors and dispersion estimates
  dge.filter <-calcNormFactors(dge.filter)
  dge.filter <-estimateCommonDisp(dge.filter, verbose=T)
  dge.filter <-estimateTagwiseDisp(dge.filter)
  
  # Calculate Differential Expression
  et<-exactTest(dge.filter )
  degSet<-topTags(et, n=Inf)@.Data[[1]]
  degSet$gene_id<-row.names(degSet) # Necessary for downstream dplyr processing
  
  degSet['Avg1']<-apply(
    rpkm(dge.filter)[,groups[[contrasts[[c]][1]]]], 1, mean
  )[degSet$gene_id]
  
  degSet['Avg2']<-apply(
    rpkm(dge.filter)[,groups[[contrasts[[c]][2]]]], 1, mean
  )[degSet$gene_id]
  
  degSet$Group_1 <-contrasts[[c]][1]
  degSet$Group_2 <-contrasts[[c]][2]
  
  createDEGSpreadSheet(
    C1 = c,                          # Name of the contrast
    dg1 = degSet,                       # Data Set for the contrast
    dg1.bioFun = bioSigRNASeq,       # Biological significance filter for dg1
    dg1.fdr = "FDR",                 # Statistic used to filter genes for dg1
    dg1.lfc = "logFC",               # Column in dg1 with log Fold Changes
    dg1.Avg1 = "Avg1",               # Column in dg1 with average value for Group_1
    dg1.Avg2 = "Avg2",               # Column in dg1 with average value for Group_2
    dg1.me = 2,                      # Min. expression for dg1.bioFun
    dg1.x = 23,                      # row number, co0rner of dg1 Summary table
    dg1.y = 2,                       # col number, corner of dg1 Summary table
    dg1.ds = "Pax6 Genes",           # short description for contrast C1 (dg1)
    template = "scripts/deg_template.xlsx",
    descPageName="Data Description", # Name of sheet to write summary tables
    wb = NULL,                       # Optionally pass a workbook object instead.
    pref = "" ,                      # Prefix for output file.
    fname=NULL,                      # Manually specify an output file name
    use_lfc = FALSE,                 # Whether to use logFC or Fold_Change
    cols=c(
      "gene_id",
      "SYMBOL", 
      "DESCRIPTION", 
      "logFC", 
      "PValue", 
      "FDR", 
      "Avg1", 
      "Avg2"
    ),
    sc_cols=c("PValue", "FDR")
  )
  fn<-paste(wd,"/results/",c,"_","DEG_Table.tsv", sep="")
  write.table(degSet, file=fn, sep="\t", quote=F, row.names = F )
}
  