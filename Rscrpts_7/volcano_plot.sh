#!/usr/bin/env Rscript
  
load("DEseq2.RData")
library(DESeq2)
library(EnhancedVolcano)
pdf("volcano_plot.pdf")

# Draw a volcano plot using EnhancedVolcano
***(***, lab = rownames(***),x = 'log2FoldChange',y = 'pvalue')

dev.off()
save.image(file = "DEseq2.RData")
