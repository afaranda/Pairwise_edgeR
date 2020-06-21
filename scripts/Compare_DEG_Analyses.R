##############################################################################
#                                                                            #
#  File: Pairwise_Tests.R                                                    #
#  Author: Adam Faranda                                                      #
#  Created: June 18, 2020                                                    #
#  Purpose: Write Pairwise analysis results to an excel spreadsheet          #
#                                                                            #
##############################################################################
library(edgeR)
library(dplyr)
options(echo=T)

# Enter Working Directory and Load Raw Data
setwd('/Users/adam/Pairwise_edgeR')
source('scripts/BuildDataMatrix.R')
source('scripts/Excel_Write_Functions.R')
source('scripts/Overlap_Comparison_Functions.R')
wd<-getwd()
results<-paste(wd,'results',sep='/')
data_dir<-paste(wd,'data',sep='/')


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

# Define A list of named contrasts; each element points to a vector with
# a pair of group labels. Positive fold changes will be associated
# with the second group listed. 
contrasts=list(
  P6vsWT_Ex=c('Wildtype', 'Pax6cKO'),
  P6vsWT_Ql=c('Wildtype', 'Pax6cKO')
)

# Define A list of named contrasts; each element points to a string wigch
# briefly describes the contrast
contrast_descriptions<-list(
  P6vsWT_Ex="Pairwise Contrast between two conditions; Exact Test",
  P6vsWT_Ql="Pairwise Contrast between two conditions; Quasi Liklihood"
)

# Define A list of named contrasts; each element points to a vector with
# a pair of file names the first filename is the template to use for
# generating the spreadsheet. The second is the name of the file with
# differential expression results. 
contrast_files<-list(
  P6vsWT_Ex=paste(
      results,'Experiment_Pax6cKO_vs_Wildtype_Exact_Test_DEG.csv', sep='/'),
  P6vsWT_Ql=paste(
      results,'Experiment_Pax6cKO_vs_Wildtype_QLFTest_DEG.csv', sep='/')
)
comparisons<-list(
  Pax6vsWT_Stat=c('P6vsWT_Ex', 'P6vsWT_Ql')
)

# Compare Differential Expression Analyses

for( c in names(comparisons)){
  C1 = comparisons[[c]][1]
  C2 = comparisons[[c]][2]
  
  # Load in data the first data file
  dg1<-read.csv(
    contrast_files[[C1]],
    row.names = 1
  )
  dg1$Group_1<-contrasts[[C1]][1]
  dg1$Group_2<-contrasts[[C1]][2]
  

  Avg1 <- paste(contrasts[[C1]][1], 'Avg', sep='_')
  Avg2 <- paste(contrasts[[C1]][2], 'Avg', sep='_')
  cols.1 = c(
    'gene_id', 'SYMBOL', 'DESCRIPTION', 'logFC', 'PValue',
    'FDR', Avg1, Avg2, 'Group_1', 'Group_2'
  )
  print(head(dg1))
  
  # Load in data the second data file
  dg2<-read.csv(
    contrast_files[[C2]],
    row.names = 1
  )
  dg2$Group_1<-contrasts[[C2]][1]
  dg2$Group_2<-contrasts[[C2]][2]
  
  Avg1 <- paste(contrasts[[C2]][1], 'Avg', sep='_')
  Avg2 <- paste(contrasts[[C2]][2], 'Avg', sep='_')
  cols.2 = c(
    'gene_id', 'SYMBOL', 'DESCRIPTION', 'logFC', 'PValue',
    'FDR', Avg1, Avg2, 'Group_1', 'Group_2'
  )
  print(head(dg2))
  
  createMethodComparisonSpreadsheet(
    C1 =C1, C2 = C2,                 # Names of each contrast
    dg1 = dg1, dg2 = dg2,  # DEG sets to compare (Full Sets)
    dg1.bioFun = bioSigRNASeq,       # Biological significance filter for dg1
    dg1.fdr = "FDR",                 # Statistic used to filter genes for dg1
    dg1.me = 2,                      # Min. expression for dg1.bioFun
    dg1.x = 39,                      # row number, corner of dg1 Summary table
    dg1.y = 2,                       # col number, corner of dg1 Summary table
    dg2.bioFun = bioSigRNASeq,       # Biological significance filter for dg2
    dg2.fdr = "FDR",                 # Statistic used to filter genes for dg2
    dg2.me = 2,                      # Min. expression for dg2.bioFun
    dg2.x = 44,                      # row number, corner of dg2 Summary table
    dg2.y = 2,                       # col number, corner of dg2 Summary table
    vns.x = 50,                      # row number, corner of ven intersect
    vns.y = 2,                       # col number, corner of stat. sig intersect
    vnb.x = 56,                      # row number, corner of ven intersect
    vnb.y = 2,                       # col number, corner of stat. sig intersect  
    ssg.x = 62,                      # row number, corner of stat. sig intersect
    ssg.y = 2,                       # col number, corner of stat. sig intersect
    bsg.x = 66,                      # row number, corner of bio. sig intersect
    bsg.y = 2,                       # col number, corner of bio. sig intersect
    dg1.ds = "TMM",             # short description for contrast C1 (dg1)
    dg2.ds = "CQN",           # short description for contrast C1 (dg2)
    template="scripts/method_comp_template.xlsx",     # Name of spreadsheet template file
    descPageName="Data Description", # Name of sheet to write summary  tables
    wb = NULL,                       # Optionally pass a workbook object instead.
    pref = "TMM_vs_CQN" ,          # Prefix for output file.
    fname=NULL,                      # Manually specify an output file name
    idc = 'gene_id',              # Column in dg1 and dg2 with unique gene id
    annot = NULL
  )
}