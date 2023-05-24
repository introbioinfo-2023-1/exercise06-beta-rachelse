#!/usr/bin/env Rscript
  
load("DEseq2.RData")
library(DESeq2)

# Extract differential expression results
# This is only necessary if the design formula is multi-factorial or, as in our case, the variable in the design formula has more than 2 levels.
# This is done with the "contrast" parameter.
# For "tissueType" perform "primary colorectal cancer" vs "normal-looking surrounding colonic epithelium"
# Fill the blanks specified with ***
deseq2Results <- results(deseq2Data, contrast=c("***", "***", "***"))

# View summary of results
summary(deseq2Results)

save.image(file = "DEseq2.RData")
