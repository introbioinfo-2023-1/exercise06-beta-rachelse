#!/usr/bin/env Rscript
  
load("DEseq2.RData")
library(DESeq2)
libaray(ggplot2)

#Order the results talbe by the smallest p value
resOrdered <- deseq2Results[order(deseq2Results$pvalue),]

#Print the row for the gene with the second smallest p-value.
#Let's call the gene as geneA.
#Fill the blank (***)
***[***,]

#Use plotCounts(built in method of DEseq2) to plot the count of geneA
#Fill the blanks (***)
counts <- ***(deseq2Data, gene="***", intgroup=c("tissueType", "individualID"), returnData=TRUE)

#Customize the plot using ggplot2
ggplot(counts, aes(x=tissueType, y=count, colour=individualID, group=individualID)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + guides(colour=guide_legend(ncol=3)) + ggtitle("geneA")

#Save the plot using ggsave
ggsave(file="count-plot.jpg", width=20, height=15)

save.image(file = "DEseq2.RData")
