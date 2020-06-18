#################################################################################
# File: BuildStringtieDataMatrix.R                                              #
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
#         June 2020                                                             #
#          - Removed GTF Processing functions                                   #
#          - Moved HTSeq aggregation functions to a separate script             #
#                                                                               #
# Author: Adam Faranda                                                          #
#################################################################################
library(dplyr)
library(data.table)
library(ballgown)
library(rtracklayer)

#################################################################################
#           Functions to build string-tie based matrices (TPM)                  #
#################################################################################

# Helper Function -- returns the last element in a delimited list
splitLast<-function(x, delim="/"){
  s<-strsplit(x, delim)
  s<-sapply(s, function(y) y[length(y)])
  s
}

st_getFileTable <-function(
  dirList,
  wd = ".",
  filename="StringTie_TPM_Files.csv"
){
  if(filename %in% list.files(wd)){
    ft<-read.csv(filename, stringsAsFactors=F)
  } 
  else{
    ft<-data.frame(
      sample=character(),
      directory=character(),
      filename=character(),
      gtffile=character(),
      stringsAsFactors=F
    )
    
    for (d in dirList){
      dr<-paste(d, list.files(d), sep = "/")
      
      x<-data.frame(
        sample = as.character(list.files(d)),
        directory = dr, 
        filename = character(length(dr)),
        gtffile =character(length(dr)),
        stringsAsFactors=F
      )
      for (r in 1:nrow(x)){
        d<-x[r, 'directory']
        fn<-paste(splitLast(d), ".txt", sep="")
        gt<-paste(splitLast(d), ".gtf", sep="")
        if(file.exists(paste(d, fn, sep="/"))){
          x[r, 'filename']<-fn
          x[r, 'gtffile'] <-gt
        }
      }
      ft<-bind_rows(ft, x)
    }
    ft<-ft%>% filter(filename != "")
    #write.csv(ft, filename, row.names=F)
  }

  return(ft)
}

# Returns a ballgown object based on data files listed in 
# the given file table (ft)
st_sewBallGown<-function(
  ft = ft,
  wd = ".",
  filename="bg.Rdata"
){
  if(filename %in% list.files(wd)){
    load(filename)
  }
  else{
    bg<-ballgown(ft$directory, pData=ft, meas="all")
    save(bg, file="bg.Rdata")
  }
  bg
}

# Imports GTF files for each sample into data frames, retruns 
# a list of data frames indexed by sample ID
st_loadFiles<-function(ft, fnCol='gtffile'){
  dataSets<-list()
  gtfNumCols<-c(
    "start", "end", "score",
    "cov", "FPKM", "TPM", "exon_number" 
  )
  for(i in 1:nrow(ft)){
    file<-paste(ft[i,2], ft[i,fnCol], sep='/')
    sample<-ft[i,1]
    print(file)
    if(grepl(".gtf",splitLast(file))){
      dataSets[[sample]]<-readGFF(file)
      dataSets[[sample]]$PK<-paste(
        dataSets[[sample]]$type,
        dataSets[[sample]]$start,
        dataSets[[sample]]$end,
        dataSets[[sample]]$transcript_id,
        sep=""
      )
      for(c in gtfNumCols){
        dataSets[[sample]][,c]<-as.numeric( dataSets[[sample]][,c])
      }
    }
    else if(grepl(".txt",splitLast(file))){
      dataSets[[sample]]<-read.table(file, header = T, stringsAsFactors = F, sep="\t")
    }
  }
  dataSets
}



# Calculate the sum of a measurement (FPKM, TPM) over genes in a data set
# returns a three column table: 
#           Gene ID, number of transcripts, summed measurement
sum_tr_over_gn<-function(df, g_idCol, measCol){
  df<-df[,c(g_idCol, measCol)]
  names(df)<-c("g_id", "meas")
  
  return(
    as.data.frame(
     df %>% 
        group_by(g_id) %>%
        summarise( Transcripts = n(), AvgMeasure = sum(as.numeric(meas)))
    )
  )
}

# Join Transcript level measurements (FPKM or TPM) into a data matrix
st_buildTranscriptMatrix<-function(ds, idCol=17, measCol=14){
  df<-data.frame(
    Ensembl=(ds[[1]] %>% filter(type == "transcript"))[,idCol],
    stringsAsFactors=F
  )
  for( i in 1:nrow(ft)){
    dg<-ds[[ft[i,1]]] %>% filter(type == 'transcript')
    if(
      length(intersect(df$Ensembl, dg[,idCol])) == length(unique(df$Ensembl))
    ){
      dg<-ds[[ft[i,1]]] %>% filter(type == 'transcript')
      df<-merge(df, dg[,c(idCol, measCol)], 
                by.x='Ensembl', by.y=1,
                sort = F)
      mcn<-names(ds[[ft[i,1]]][measCol])
      names(df)[grep(mcn, names(df))]<-ft[i,1]
      print(i)
    }
    else{
      print(paste("Can't join sample:",ft[i,1], "ID Mismatch" ))
    }
    
  }
  df<- df %>% filter(!grepl("__",Ensembl))
  row.names(df)<-df$Ensembl
  df<-df[order(df[,"Ensembl"]),]
  dm<-as.matrix(df[,!grepl("Ensembl", names(df))])
  list(ft, dm)
}

# Join Gene level measurements (FPKM or TPM) into a data matrix
# transcript level measurements are summed over their respective genes
st_buildGeneMatrix<-function(ds, g_idCol=9, measCol=14){
  df<-data.frame(
    Ensembl=unique(ds[[1]][,g_idCol]),
    stringsAsFactors=F
  )
  for( i in 1:nrow(ft)){
    dg<-ds[[ft[i,1]]] %>% filter(type == 'transcript')
    dg<-sum_tr_over_gn(dg, g_idCol=g_idCol, measCol=measCol)
    if(
      length(intersect(df$Ensembl, dg[,1])) == length(unique(df$Ensembl))
    ){
      
      df<-merge(df, dg[,c(1, 3)], 
                by.x='Ensembl', by.y=1,
                sort = F)
      mcn<-names(dg[3])
      names(df)[grep(mcn, names(df))]<-ft[i,1]
      print(i)
    }
    else{
      print(paste("Can't join sample:",ft[i,1], "ID Mismatch" ))
    }
    
  }
  row.names(df)<-df$Ensembl
  df<-df[order(df[,"Ensembl"]),]
  dm<-as.matrix(df[,!grepl("Ensembl", names(df))])
  list(ft, dm)
}


# Validate lengthTable Counts and lengths internally
#any(lt$Xscript_Exons > lt$Union_Exons)
#any(lt$Xscript_Length > lt$Union_Length)

# Iterate over genes and validate The number of unique exons in each 
# check<-data.frame(
#   gene_id=character(),
#   Ex = numeric(),
#   stringsAsFactors = F
# )
# for (gn in 5000:5500){
#   print(gn)
#   g<-gn_len$gene_id[gn]
#   check[gn,1]<-g
#   check[gn,2]<-length(
#     unique(
#       orig_gtf[orig_gtf$gene_id == g & orig_gtf$type == "exon", "exon_id"]
#       )
#     )
# }
# inner_join(
#   gn_len,
#   check,
#   by="gene_id"
# ) %>%   
#   mutate(ck=abs(Union_Exons - Ex)) %>%
#   summarize(sum(ck))


