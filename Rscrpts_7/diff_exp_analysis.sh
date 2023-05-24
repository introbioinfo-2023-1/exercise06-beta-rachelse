#!/usr/bin/env Rscript
  
load("DEseq2.RData")
library(DESeq2)

# Run pipeline for differential expression steps
# Fill the blank specified with ***
*** <- DESeq(***)

save.image(file = "DEseq2.RData")
