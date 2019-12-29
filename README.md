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


3. Make updates to the script "Pairwise_edgeR_ExactTest.R"
    * update the variables "edb", "groups", and "contrasts" to
  reflect the desired sample assignment and contrasts
    * Update the call to function "createDEGSpreadSheet" to match the column
  configuration desired for outpout spreadsheets.
    * If needed, run the first 22 lines interactively to generate a sample
  manifest file.  This file can be edited in a spreadsheet editor to add
  sample grouping information.

4. Run the script "Pairwise_edgeR_ExactTest.R" to generate spreadsheets with Differential expression statistcs calculated by edgeR
    * Future verisions will also generate a pdf report that has dispersion and
    MDS plots and summary DE statistics
    
5. Manually edit and review experiment-specific fields in the generated
spreadsheet