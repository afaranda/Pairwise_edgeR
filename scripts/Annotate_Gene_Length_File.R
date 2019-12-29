library(edgeR)
library('AnnotationHub')
library(dplyr)
source('scripts/BuildDataMatrix.R')
options(echo=T)

# Enter Working Directory and Load Raw Data
# setwd('/Users/afaranda/Desktop/LEC_Time_Series')
wd<-getwd()
results<-paste(wd,'results',sep='/')
#data_dir<-paste(wd,'data',sep='/')
data_dir <-'~/Documents/LEC_Time_Series_HTSeq_Counts'

# Generate table of data files / sample information
# After generating the file table, update it with 
# Experimental Grouping information in excel/oocalc/emacs
ft<-hc_getFileTable(
  dirList=c(data_dir),
  filename = "HtSeq_GeneCountFiles.csv"
)

# Assemble Data Sets
ds<-hc_loadFiles(ft)
ft<-hc_identifierConsistency(ds, ft, idCol=1)

# Import "Gene Length" Annotation to add to dgelist
gt<-read.table(
  '~/Desktop/gene_coding_lengths.txt', 
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
detach(package:AnnotationHub, unload=T)
detach(package:ensembldb, unload=T)
detach(package:AnnotationFilter, unload=T)


# Assemble master "DGE List" Object. 
df=hc_buildDataFrame(ds, ft, return_matrix = F)

master<-DGEList( 
  counts = df,
  # genes = gt[row.names(df),],
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

# # Define Replicate Groups
# sampleGroups<-list(
#  	WT_0_Hour=c("S1", "S2", "S3"),
#  	WT_48_Hour=c("S7", "S8", "S9"),
#  	FN_0_Hour=c("S10", "S11", "S12"),
#  	FN_48_Hour=c("S13", "S14", "S15")
# )
#  
# pairwiseContrasts<-list(
# 	WT_0_Hour_vs_WT_48_Hour = c('WT_0_Hour', 'WT_48_Hour'),
# 	FN_0_Hour_vs_FN_48_Hour = c('FN_0_Hour', 'FN_48_Hour'),
# 	WT_0_Hour_vs_FN_0_Hour = c('WT_0_Hour', 'FN_0_Hour'),
# 	WT_48_Hour_vs_FN_48_Hour = c('WT_48_Hour', 'FN_48_Hour')
# )

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
  
  # Add Annotations
  degSet <- dplyr::inner_join(
    degSet, 
  )
}
  
# }
  
  
  
  dg1<-left_join(
    dg1, an[,c("MGI.symbol", "description")],
    by="MGI.symbol"
  )
  createDEGSpreadSheet(
    C1 = p,                          # Name of the contrast
    dg1 = dg1,                       # Data Set for the contrast
    dg1.bioFun = bioSigRNASeq,       # Biological significance filter for dg1
    dg1.fdr = "FDR",                 # Statistic used to filter genes for dg1
    dg1.lfc = "logFC",               # Column in dg1 with log Fold Changes
    dg1.Avg1 = "Avg1",               # Column in dg1 with average value for Group_1
    dg1.Avg2 = "Avg2",               # Column in dg1 with average value for Group_2
    dg1.me = 2,                      # Min. expression for dg1.bioFun
    dg1.x = 23,                      # row number, corner of dg1 Summary table
    dg1.y = 2,                       # col number, corner of dg1 Summary table
    dg1.ds = "Pax6 Genes",           # short description for contrast C1 (dg1)
    template = gsub(
      "-X-", p,
      "Pax6_-X-_deg_template.xlsx"
    ),
    descPageName="Data Description", # Name of sheet to write summary tables
    wb = NULL,                       # Optionally pass a workbook object instead.
    pref = "Pax6" ,                  # Prefix for output file.
    fname=NULL,                      # Manually specify an output file name
    use_lfc = FALSE,                 # Whether to use logFC or Fold_Change
    cols=c(
      "MGI.symbol", 
      "description", 
      "logFC", 
      "p_value", 
      "FDR", 
      "Avg1", 
      "Avg2"
    )
  )
}






  
	# Define a factor for group assignment in the DGEList
	group<-factor(
		c(
			rep(1, length(sampleGroups[[g[1]]])), 
			rep(2, length(sampleGroups[[g[2]]]))
		), 
		labels=c(g[1], g[2])
	)
	
	# Build DGEList
	y<-DGEList(
		counts=df[,samples], 
		group=group, 
		genes=df[,c("PrimaryKey", "GeneID", "Chr", "CodingLength")]
	)
	row.names(y)<-y$genes$PrimaryKey

	# Filter DGEList: Remove genes where fewer than two samples have a cpm > 1
	keep <- rowSums(cpm(y) > 1) >= 2
	y.filter <-y[keep, ]
	y.filter$samples$lib.size <- colSums(y.filter$counts) # Fix library size after filtering
	
	# Calculate Normalization factors and dispersion estimates
	y.filter <-calcNormFactors(y.filter)
	y.filter <-estimateCommonDisp(y.filter, verbose=T)
	y.filter <-estimateTagwiseDisp(y.filter)
	
	# Calculate Differential Expression
	et<-exactTest(y.filter, pair=levels(group))
	degSet<-topTags(et, n=25000)@.Data[[1]]

	# Add Absolute Fold Change column "FC"
	degSet$FC<-sapply(degSet$logFC, calcFC)
	print(paste("Rows in degSet:", nrow(degSet)))
	
	# Join biomart query results to table of differentially expressed genes
	# on the gene symbol attribute. 
	degSet<-merge(degSet, result, by.x='GeneID', by.y='external_gene_name', all.x=T)
	file<-paste(outDir,"/",
		paste(g, collapse="_vs_"),
		"_expressedTags-all.txt", sep=""
	)
	print(paste("Rows in degSet after biomaRt:", nrow(degSet)))
	
	# Reorder Columns
	degSet<-degSet %>% dplyr::select(PrimaryKey,
				GeneID, ensembl_gene_id, description,
				FC, logFC, logCPM, PValue, FDR, CodingLength, Chr
	)
	
	# Join counts per million (calculated prior to filtering genes)
	cpmy<-as.data.frame(cpm(y))
	names(cpmy)<-gsub("Count","Count.cpm", names(cpmy))
	cpmy$PrimaryKey<-row.names(y)
	degSet<-merge(x=degSet, y=cpmy, by='PrimaryKey')
	print(paste("Rows after join cpm:", nrow(degSet)))

	# Join read counts
	counts<-as.data.frame(y$counts)
	names(counts)<-gsub("Count", "Count.counts", names(counts))
	counts$PrimaryKey<-row.names(counts)
	degSet<-merge(x=degSet, y=counts, by='PrimaryKey')
	print(paste("Rows after join counts:", nrow(degSet)))

	# Join rpkm
	degSet<-merge(x=degSet, y=df[,c("PrimaryKey", gsub("_GeneCount", "_RPKM", samples))], by="PrimaryKey")
	print(paste("Rows after join RPKM:", nrow(degSet)))
	
	write.table(degSet,
		file=file,
		sep="\t", row.names=F, quote=F
	)
}


print(sessionInfo())