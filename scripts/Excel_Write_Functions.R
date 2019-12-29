################################################################################
# File: Pax6Epi_Multi_Write_Excel.R                                            #
# Author: Adam Faranda                                                         #
# Created: June 28, 2019                                                       #
# Purpose: Functionds that write tables to spreadsheets                        #
#                                                                              #
################################################################################

# Setup Workspace   
options(echo=T)
library('openxlsx')
library('dplyr')
library('org.Mm.eg.db')
source("scripts/Overlap_Comparison_Functions.R")

# Create Styles for text and numeric columns
tt<-createStyle(
  numFmt="TEXT", halign = "left", textDecoration = "bold", #For table title
  wrapText = TRUE
)
tr<-createStyle
tb<-createStyle(numFmt = "0", halign = "center")     # Use for table contents
tc<-createStyle(numFmt="TEXT", halign = "center")    # Use with column headers
tl<-createStyle(numFmt="TEXT", halign = "left")      # Use with Gene Symbols
nm<-createStyle(numFmt="0.00", halign = "center")    # Use with data tables
sc<-createStyle(numFmt="0.00E+0", halign = "center") # Use with p-value, FDR



# general function to write a data table to a page
writeSheet<-function(
  wb, df, name,    # wb: workbook object, df: a data frame, name: sheet name
  tx_cols=c(      
    "MGI.symbol", 
    "description",                    # Columns that contain text
    "Agreement"
  ),
  nm_cols=c("logFC", "Avg1", "Avg2"), # Columns that contain numeric data
  sc_cols=c("p_value", "FDR")         # Columns requiring scientific notation 
){
  if(!(name %in% wb$sheet_names)){
    addWorksheet(wb, sheetName = name)
  }
  hed.cells<-expand.grid(
    row=1, 
    col=1:ncol(df)
  )
  txt.cells<-expand.grid(
    row=2:(nrow(df)+1), 
    col=colNum(df, tx_cols)
  )
  num.cells<-expand.grid(
    row=2:(nrow(df)+1), 
    col=colNum(df, nm_cols)
  )
  sci.cells<-expand.grid(
    row=2:(nrow(df)+1), 
    col=colNum(df, sc_cols)
  )
  addStyle(wb, name, rows=hed.cells$row, cols=hed.cells$col, style=tc)
  addStyle(wb, name, rows=txt.cells$row, cols=txt.cells$col, style=tl)
  addStyle(wb, name, rows=num.cells$row, cols=num.cells$col, style=nm)
  addStyle(wb, name, rows=sci.cells$row, cols=sci.cells$col, style=sc)
  
    wc<-data.frame(lapply(df, as.character), stringsAsFactors=F)
    wc<-rbind(wc, names(wc))
    wc<-apply(wc, 2, function(x) max(nchar(x), na.rm=T))
  
  setColWidths(
    wb, sheet = name, cols = 1:ncol(df), 
    widths = sapply(wc, function(x) min(x+4, 36))
  )
  
  writeData(wb, name, df)
}

# Function to Generate a Description Page: Summary statististics are passed
# as a nested list, C1 and C2 are contrast names and start_cell is coordinates of
# the top left corner of the first table (ignored if a template is used). 

############### Sample "descTables" nested list object ######################
#  Each element in the top level list is a table. For each table,           #
#  the following elements are included:                                     # 
#            Table: a data frame containing the data to be printed          #
#            corner: The top left coordinates of the table                  #
#            cn: boolean flag -- whether to print column headers            #
#            rn: boolean flag -- whether to print row names                 #
#            tc: boolean flag -- prin the title in corner, or above table   #
#############################################################################

# descTables<-list(
#   Table_One=list(Table=df, corner=c(1, 1), cn=F, rn=F, tc=F),
#   Table_Two=list(Table=df, corner=c(6, 1), cn=F, rn=F, tc=T),
#   Table_Three=list(Table=df, corner=c(11, 1), cn=F, rn=T, tc=F),
#   Table_Four=list(Table=df, corner=c(16, 1), cn=F, rn=T, tc=T),
#   Table_Five=list(Table=df, corner=c(21, 1), cn=T, rn=F, tc=F),
#   Table_Six=list(Table=df, corner=c(26, 1), cn=T, rn=F, tc=T),
#   Table_Seven=list(Table=df, corner=c(31, 1), cn=T, rn=T, tc=F),
#   Table_Eight=list(Table=df, corner=c(36, 1), cn=T, rn=T, tc=T)
# )

writeDescTables<-function(
  descTables=list(), name="Description",
  template=NULL, offset=0, wb=NULL
  
){
  # Setup Workbook
  if(is.null(wb) & is.null(template)){
    wb<-createWorkbook()
  } else if(is.null(wb) & !is.null(template)){
    
    # If a template is supplied -- populate tables in the list "descTables"
    # starting from the coordinates in "start_cell", 
    print(template)
    if(file.exists(template)){
    	  print("Opening template")
      wb<-loadWorkbook(template)
      if(!(name %in% wb$sheet_names)){
        addWorksheet(wb, sheetName = name)
      }
      for(t in names(descTables)){
        print(t)
        tbl<-descTables[[t]]
        if(!tbl$tc){
          rn<-ifelse(tbl$rn, 1, 0)
          cn<-ifelse(tbl$cn, 1, 0)
          # Generate Table Title
          removeCellMerge(
            wb, name, rows=tbl$corner[1], 
            cols=(tbl$corner[2] + offset):(ncol(tbl$Table)+rn)
          )
          mergeCells(
            wb, name, rows=tbl$corner[1], 
            cols=(tbl$corner[2] + offset):(ncol(tbl$Table)+rn)
          )
          cells<-expand.grid(
            rows=tbl$corner[1], 
            cols=(tbl$corner[2] + offset):(ncol(tbl$Table)+rn)
          )
          addStyle(
            wb, name, tt, rows=cells$rows, 
            cols= cells$cols,
            gridExpand = T
          )
          writeData(
            wb, name, as.data.frame(t), colNames = F, rowNames = F,
            startRow = (tbl$corner[1]), startCol = offset + tbl$corner[2]
          )
          
          # Generate Table Body
          cells<-expand.grid(
            rows=(tbl$corner[1] + 1):(tbl$corner[1] + nrow(tbl$Table)+cn),
            cols=(
              tbl$corner[2] + offset + rn):(ncol(tbl$Table)+rn+offset)
          )
          print(paste(min(cells$rows), min(cells$cols)))
          print(paste(max(cells$rows), max(cells$cols)))
          addStyle(
            wb, name, tb,
            rows=cells$rows,
            cols=cells$cols
          )
          writeData(
            wb, name, tbl$Table,
            startRow =(tbl$corner[1] + 1), 
            startCol =(tbl$corner[2] + offset),
            colNames = tbl$cn, rowNames = tbl$rn
          )
        } else {
          rn<-ifelse(tbl$rn, 1, 0)
          cn<-ifelse(tbl$cn, 1, 0)
          
          # Generate Table Body
          cells<-expand.grid(
            rows=(tbl$corner[1]):(tbl$corner[1] + nrow(tbl$Table)),
            cols=(rn+tbl$corner[2] + offset):(rn+ncol(tbl$Table)+offset)
          )
          print(paste(min(cells$rows), min(cells$cols)))
          print(paste(max(cells$rows), max(cells$cols)))
          addStyle(
            wb, name, tb,
            rows=cells$rows,
            cols=cells$cols
          )
          writeData(
            wb, name, tbl$Table,
            startRow =(tbl$corner[1]), 
            startCol =(tbl$corner[2] + offset),
            colNames = tbl$rn, rowNames = tbl$rn
          )
          # Generate Table Title 
          
          cells<-expand.grid(
            rows=tbl$corner[1], 
            cols=(tbl$corner[2] + offset)
          )
          addStyle(
            wb, name, tt, rows=cells$rows, 
            cols= cells$cols,
            gridExpand = T
          )
          writeData(
            wb, name, as.data.frame(t), colNames = F, rowNames = F,
            startRow = (tbl$corner[1]), startCol = offset + tbl$corner[2]
          )
          
          print("wtf")
        }
      }
      
      
    }
  }
  if(!(name %in% wb$sheet_names)){
    addWorksheet(wb, sheetName = name)
  }
  return(wb)
}



# Function that Generates a spreadsheet comparing results from two
# pairsiwe contrasts. 

createBioSigOverlapSpreadSheet<-function(
  C1 = "AvsB",  C2 = "CvsD",       # Names of ach contrast
  dg1 = dg1, dg2 = dg2,            # DEG sets to compare
  dg1.bioFun = bioSigRNASeq,       # Biological significance filter for dg1
  dg1.fdr = "FDR",                 # Statistic used to filter genes for dg1
  dg1.me = 2,                      # Min. expression for dg1.bioFun
  dg1.x = 30,                      # row number, corner of dg1 Summary table
  dg1.y = 2,                       # col number, corner of dg1 Summary table
  dg2.bioFun = bioSigRNASeq,       # Biological significance filter for dg2
  dg2.fdr = "FDR",                 # Statistic used to filter genes for dg2
  dg2.me = 2,                      # Min. expression for dg2.bioFun
  dg2.x = 35,                      # row number, corner of dg2 Summary table
  dg2.y = 2,                       # col number, corner of dg2 Summary table
  ssg.x = 41,                      # row number, corner of stat. sig intersect
  ssg.y = 2,                       # col number, corner of stat. sig intersect
  bsg.x = 45,                      # row number, corner of bio. sig intersect
  bsg.y = 2,                       # col number, corner of bio. sig intersect
  dg1.ds = "Pax6 Genes",           # short description for contrast C1 (dg1)
  dg2.ds = "Runx1 Genes",          # short description for contrast C1 (dg2)
  template="Comparisons.xlsx",     # Name of spreadsheet template file
  descPageName="Data Description", # Name of sheet to write summary  tables
  wb = NULL,                       # Optionally pass a workbook object instead.
  pref = "FuncTest" ,              # Prefix for output file.
  fname=NULL                       # Manually specify an output file name
){
  
  allResults<-query(dg1, dg2)
  bioResults<-query(
    dg1.bioFun(dg1, minExp = dg1.me),
    dg2.bioFun(dg2, minExp = dg2.me)
  )
  
  print(nrow(allResults))
  print(nrow(bioResults))
  
  stat.tables<-subsetTables(
    Contrast_1 = C1, Contrast_2 = C2,
    allResults, annot = an, unlog=T, descname = T
  )
  stat.inx<-tabulateOverlap(stat.tables, rename = T)
  
  bio.tables<-subsetTables(
    Contrast_1 = C1, Contrast_2 = C2,
    bioResults, annot = an, unlog=T, descname = T, stat = F
  )
  bio.inx<-tabulateOverlap(bio.tables, rename = T)
  
  print(stat.inx)
  print(bio.inx)
  
  # Set up list
  descTables = list(
    C1=list(
      Table=degSummary(dg1, fdr=dg1.fdr, minExp=dg1.me), corner=c(dg1.x, dg1.y), cn=T, rn=F, tc=F
    ),
    C2=list(
      Table=degSummary(dg2, fdr=dg2.fdr, minExp=dg2.me), corner=c(dg2.x, dg2.y), cn=T, rn=F, tc=F
    ),
    `Statistically Significant Intersection`=list(
      Table=stat.inx, corner=c(ssg.x, ssg.y), cn=T, rn=T, tc=T
    ),
    `Biologically Significant Intersection`=list(
      Table=bio.inx, corner=c(bsg.x, bsg.y), cn=T, rn=T, tc=T
    )
  )
  names(descTables)[grep("C1", names(descTables))]<-paste(dg1.ds, C1)
  names(descTables)[grep("C2", names(descTables))]<-paste(dg2.ds, C2)
  
  
  # Write Description Tables -- see script "Excel_Write_Functions.R"
  wb<-writeDescTables(
    template=template,          # Name of template file
    name = descPageName,        # Name of sheet with descriptive tables
    descTables = descTables
  )
  
  # Delete "Contrasts tab from directional subsets
  stat.tables["Contrasts"]<-NULL
  bio.tables["Contrasts"]<-NULL
  
  # Add statistically significant directional subsets to workbook object
  for(i in names(stat.tables)){
    writeSheet(
      wb, stat.tables[[i]], i,
      nm_cols=c("Change", "Avg")
    )
  }
  
  # Add biologically significant directional subsets to workbook object
  for(i in names(bio.tables)){
    writeSheet(
      wb, bio.tables[[i]], i,
      nm_cols=c("Change", "Avg")
    )
  }
  # Save workbooks to file  
  # fn<-paste(names(pairwise)[pairwise %in% comparisons[[c]]], collapse="_")
  if(is.null(fname)){
  	fn<-paste(C1, C2, sep="_")
  	fn<-paste(pref,fn, "DEG_Comparison.xlsx", sep="_")
  	print(fn)
  } else {
  	fn<-fname
  }	
  saveWorkbook(wb, file=fn, overwrite = T)
  wb
}

# Function that Generates a spreadsheet for results from a pariwise contrast
createDEGSpreadSheet<-function(
  C1 = "AvsB",                     # Name of the contrast
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
  template="Contrast.xlsx",        # Name of spreadsheet template file
  descPageName="Data Description", # Name of sheet to write summary tables
  wb = NULL,                       # Optionally pass a workbook object instead.
  pref = "" ,                      # Prefix for output file.
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
  ),
  nm_cols=c("Change", "Avg"),      # Columns being formated numerically
  sc_cols=c("p_value", "FDR")      # Columns with scientific notation formatting
){
  gr1<-paste(unique(dg1$Group_1),"Avg", sep="_")
  gr2<-paste(unique(dg1$Group_2),"Avg", sep="_")
  
  dg1<-dg1[,cols]
  names(dg1)[grep(dg1.Avg1, names(dg1))]<-gr1
  names(dg1)[grep(dg1.Avg2, names(dg1))]<-gr2
  
  allResults<-dg1.bioFun(dg1, minExp = 0, a1 = gr1, a2=gr2)
  bioResults<-dg1.bioFun(dg1, minExp = dg1.me, a1 = gr1, a2=gr2)
  
  print(nrow(allResults))
  print(nrow(bioResults))
  
  # Set up list
  descTables = list(
    C1=list(
      Table=degSummary(dg1, fdr=dg1.fdr, minExp=dg1.me, Avg1=gr1, Avg2=gr2), 
      corner=c(dg1.x, dg1.y), cn=T, rn=F, tc=F
    )
   )
  names(descTables)[grep("C1", names(descTables))]<-paste(dg1.ds, C1)
  
  
  
  # Write Description Tables -- see script "Excel_Write_Functions.R"
  wb<-writeDescTables(
    template=template,          # Name of template file
    name = descPageName,        # Name of sheet with descriptive tables
    descTables = descTables
  )
    
  tables = list(
  	`All Genes` = dg1,
    `Statistically Significant` = allResults,
    `Biologically Significant` = bioResults
  )
  # If needed, Convert logFC to Fold_Change
  if(!use_lfc){
  	for(i in names(tables)){
  		tables[[i]]<-unlog(tables[[i]], dg1.lfc)
  		names(tables[[i]])[grep(dg1.lfc, names(tables[[i]]))]<-"Fold_Change"
  	}
  }
  
  for(i in names(tables)){
    writeSheet(
      wb, tables[[i]], i,
      nm_cols=nm_cols,
      sc_cols=sc_cols
    )
  }
  
  # Save workbooks to file  
  # fn<-paste(names(pairwise)[pairwise %in% comparisons[[c]]], collapse="_")
  if(is.null(fname)){
  	fn<-C1
  	fn<-paste(pref,fn, "Differential_Expression.xlsx", sep="_")
  	print(fn)
  } else {
  	fn<-fname
  }	
  saveWorkbook(wb, file=fn, overwrite = T)
  wb
}












