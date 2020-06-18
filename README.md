## Pairwise_edgeR
This repository provides a set of scripts and templates for running the
exact test provided by edgeR on one or more pairs of RNA Sequencing data
to apply this workflow:

1. Clone this repository


2. Copy datafiles with readcounts for each sample into the 'data' folder.
    * HTSeq-count generates a two column table with feature ID's in the first
    column and read counts in the second column
    * Stringtie import for edgeR analysis is not currently implemented; but
    will be implemented using the tximport package


3. Make updates to the script "Pairwise_Tests.R"
    * update the variables "edb", "groups", and "contrasts" to
    reflect the desired sample assignment and contrasts
    * update any output file names for generated plots and data tables
    * If needed, run the first 36 lines interactively to generate a sample
    manifest file.  This file can be edited in a spreadsheet editor to add
    sample grouping information.

4. Run the script "Pairwise_Tests.R" to generate spreadsheets with 
Differential expression statistcs calculated by edgeR
    * All Samples provided are combined into a single DGEList
    * Low-count filtering, Normalization Factor and Dispersion estimates
    are calculated for the full DGEList.
    * The script generates DEG Tables (csv files) for the 
    Quasi-likelihood test and the Exact Test based on user specified 
    pairwise contrasts.
    * The script generates a BCV and MDS plot for the full DGEList
    
5. Use Write_DEG_Spreadsheet.R can be used to generate excel spreadsheets
    * Edit this script to indicate which DEG Tables (found in the results 
    folder) to generate spreadsheets for.
    * There are tons of dumb little parameters
    
    
    