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
  Wildtype=paste('WT',1:3,sep=''),
  Pax6cKO=paste('P6',1:3,sep='')
)

# Define A list of named contrasts; each element points to a vector with
# a pair of group labels. Positive fold changes will be associated
# with the second group listed. 
contrasts=list(
  P6vsWT_Exact=c('Wildtype', 'Pax6cKO'),
  P6vsWT_QLF=c('Wildtype', 'Pax6cKO')
)

# Define A list of named contrasts; each element points to a string wigch
# briefly describes the contrast
contrast_descriptions<-list(
  P6vsWT_Exact="Pairwise Contrast between two conditions; Exact Test",
  P6vsWT_QLF="Pairwise Contrast between two conditions; Quasi Liklihood"
)

# Define A list of named contrasts; each element points to a vector with
# a pair of file names the first filename is the template to use for
# generating the spreadsheet. The second is the name of the file with
# differential expression results. 
contrast_files<-list(
  P6vsWT_Exact=c(
    template='scripts/deg_template.xlsx',
    filename=paste(
      results,'Experiment_Pax6cKO_vs_Wildtype_Exact_Test_DEG.csv', sep='/')
  ),
  P6vsWT_QLF=c(
    template='scripts/deg_template.xlsx',
    filename=paste(
      results,'Experiment_Pax6cKO_vs_Wildtype_QLFTest_DEG.csv', sep='/')
  )
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
  
  #
  degSet<-read.csv(
    contrast_files[[c]]['filename'],
    row.names = 1
  )
  degSet$Group_1<-contrasts[[c]][1]
  degSet$Group_2<-contrasts[[c]][2]
  

  Avg1 <- paste(contrasts[[c]][1], 'Avg', sep='_')
  Avg2 <- paste(contrasts[[c]][2], 'Avg', sep='_')
  cols = c(
    'gene_id', 'SYMBOL', 'DESCRIPTION', 'logFC', 'PValue',
    'FDR', Avg1, Avg2, 'Group_1', 'Group_2'
  )

  createDEGSpreadSheet(
    C1 = c,                          # Name of the contrast
    dg1 = degSet,                    # Data Set for the contrast
    dg1.bioFun = bioSigRNASeq,       # Biological significance filter for dg1
    dg1.fdr = "FDR",                 # Statistic used to filter genes for dg1
    dg1.lfc = "logFC",               # Column in dg1 with log Fold Changes
    dg1.Avg1 = Avg1,                 # Column in dg1 with average value for Group_1
    dg1.Avg2 = Avg2,                 # Column in dg1 with average value for Group_2
    dg1.me = 2,                      # Min. expression for dg1.bioFun
    dg1.x = 23,                      # row number, corner of dg1 Summary table
    dg1.y = 2,                       # col number, corner of dg1 Summary table
    dg1.ds = contrast_descriptions[[c]], # short description for contrast C1 (dg1)
    template = "scripts/deg_template.xlsx",
    descPageName="Data Description", # Name of sheet to write summary tables
    wb = NULL,                       # Optionally pass a workbook object instead.
    pref = "" ,                      # Prefix for output file.
    fname=NULL,                      # Manually specify an output file name
    use_lfc = FALSE,                 # Whether to use logFC or Fold_Change
    cols=setdiff(                    # Names of columns to keep in final tables
      cols,
      c('Group_1', 'Group2')
    ),
    sc_cols=c("PValue", "FDR")
  )
}