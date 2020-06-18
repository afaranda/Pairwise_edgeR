#################################################################################
# File: BuildHTSeq_DataMatrix.R                                                 #
# Purpose: Implement functions to traverse a set of data directories, find      #
#          count data files for individual samples and join them into a data    #
#          matrix                                                               #
# Created: April 30, 2019                                                       #
# Changes:                                                                      #
#         November 2019                                                         #
#          - Replaced Function that estimates exon-union lengths from a gtf     #
#                                                                               #
#         June 2019                                                             #
#          - Add functionality for processing StringTie results                 #
#          - update function names to distinguish htseq methods from stringtie  #
#          - functions prefixed with "hc_" process htseq-count data             #
#          - functions prefixed with "st_" process stringtie data               #
#                                                                               #
#         Jun 2020                                                              #
#          - Moved Stringtie Aggregation functions to a separate script         #
#          - Removed GTF Processing functions                                   #
#                                                                               #
#                                                                               #
# Author: Adam Faranda                                                          #
#################################################################################
library(dplyr)
library(data.table)

#################################################################################
#                 Functions to build HTSeq-Count Matrices                       #
#################################################################################

# Iterate over a list of data directories containing htseq-count output files, 
# return a two column table of datafiles in each directory that match a target 
# pattern

hc_getFileTable <-function(
  wd = ".",
	dirList, pattern='_GeneCount.txt', 
	filename="HtSeq_GeneCountFiles.csv"
){
	if(filename %in% list.files(wd)){
	  fn = paste(wd, filename, sep="/")
		ft<-read.csv(fn, stringsAsFactors = F)
	} 
	else{
		ft<-data.frame(
		  sample = character(),
			directory=character(),
			filename=character(),
			stringsAsFactors=F
		)
		
		for (d in dirList){
			 list.files(d, pattern)
			 x<-data.frame(
			  sample =  gsub(pattern, "", list.files(d, pattern)),
			 	directory = d,
			 	filename = list.files(d, pattern),
			 	stringsAsFactors=F
			 )
			 ft<-rbind(ft, x)
		}
		fn = paste(wd, filename, sep="/")
		write.csv(ft, fn, row.names=F)
	}
	return(ft)
}

# Given a file table Load datafiles into a list object
hc_loadFiles<-function(ft){
	dataSets<-list()
	for(i in 1:nrow(ft)){
		file<-paste(ft[i,2], ft[i,3], sep='/')
		print(ft[i,1])
		sample<-ft[i,1]
		print(file)
		dataSets[[sample]]<-read.table(
			file, header=F, stringsAsFactors=F, sep='\t'
		)
	}
	dataSets
}

# Check Identifier Consistency by iterating over each datafile
hc_identifierConsistency<-function(ds, ft, idCol=1){
	ft$rowCount<-0
	ft$uniqueID<-0
	check<-ds[[1]][,idCol]
	for (i in 1:nrow(ft)){
		ft$rowCount<-nrow(ds[[ft[i,1]]])
		ft$uniqueID<-length(
						unique(
							ds[[ft[i,1]]][,idCol]
						)
		)
		if (sum(check != ds[[ft[i,1]]][,idCol]) > 0){
			print(paste("ID mismatch between first sample and:", sample))
		}
	}
	ft
}

# Join count data columns into a matrix
hc_buildDataFrame<-function(ds, ft, idCol=1, measCol=2, return_matrix=T){
	df<-data.frame(Ensembl=ds[[1]][,idCol], stringsAsFactors=F)
	for( i in 1:nrow(ft)){
		if(!any(df$Ensembl != ds[[ft[i,1]]][,idCol])){
			df<-merge(df, ds[[ft[i,1]]][,c(idCol, measCol)], 
			by.x='Ensembl', by.y=1,
			sort = F)
			names(df)[grep('V2', names(df))]<-ft[i,1]
		}
		else{
			print(sum(df$Ensembl != ds[[ft[i,1]]][,idCol]))
			print(paste("Can't join sample:",sample, "ID Mismatch" ))
		}
	
	}
	df<- df %>% filter(!grepl("__",Ensembl))
	row.names(df)<-df$Ensembl
	df<-df[order(df[,"Ensembl"]),]
	if(return_matrix){
	 return(as.matrix(df[,!grepl("Ensembl", names(df))]))
	} else {
	  return(df[,!grepl("Ensembl", names(df))])
	}
}

# Compile expression data into a DGE list
buildDGE<-function(mat, ft=ft, ft_id="sample", gt=gt, gt_id="gene_id"){
  # Verify that all ft_id and gn_id are unique in their respective tables
  # and that both match the count matrix
  if(length(unique(ft[,ft_id])) == nrow(ft)){
    if(length(unique(gt[,gt_id])) == nrow(gt)){
      if(!any(colnames(mat) != ft[, ft_id])){
        if(!any(rownames(mat) != gt[,gt_id])){
          row.names(ft)<-ft[,ft_id]
          row.names(gt)<-gt[,gt_id]
          dgelist<-DGEList(
            counts = mat,
            samples = ft,
            genes = gt
          )
        } else if(length(setdiff(rownames(mat), genes[,gt_id])) == 0){
          print("Sorting rows by Gene ID ")
          row.names(ft)<-ft[,ft_id]
          row.names(genes)<-gt[,gt_id]
          gt<-genes[order(gt[,gt_id])]
          mat<-mat[order(rownames(mat)),]
          dgelist<-DGEList(
            counts = mat,
            samples = ft,
            genes = genes
          )
        } else {
          print("Genes given did not match data matrix")
          return(NULL)
        }
      } else if(length(setdiff(colnames(mat), ft[,ft_id])) == 0) {
        print("Sorting Columns by Sample ID and Genes by Gene ID")
        row.names(ft)<-ft[,id_col]
        row.names(genes)<-genes[,id_col]
        
        ft<-ft[order(ft[,ft_id]),]
        genes<-genes[order(genes[,id_col])]
        
        mat<-mat[,order(colnames(mat))]
        mat<-mat[order(rownames(mat)),]
        dgelist<-DGEList(
          counts = mat,
          samples = ft,
          genes = genes
        )
      } else {
        print("Columns in Sample data do not match columns in count matrix")
        return(NULL)
      }
    }
  }
}

# Drop HTSeq-Samples
hc_dropSamples<-function(mat, ft=ft, id_col=1,samples=c()){
  ft<-ft[!(ft[,id_col] %in% samples),]
  s<-ft[,id_col] 
  mat<-mat[,s]
  if(!any(colnames(mat)!= ft$sample)){
    return(list(ft, mat))
  } else {
    print("Column Matching Error")
  }
}
