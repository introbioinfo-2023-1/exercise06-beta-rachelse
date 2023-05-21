# command06.sh
# 1. Join the `./result/GSM461177.tsv` and `./result/GSM461180.tsv` with the `Gene ID` column.
#    Only consider genes with gene names (Ignore transcripts of which gene names are "-").

# 2. Calculate the log2 fold change of all genes between wild-type sample
#    and the *pasilla* RNAi sample and find the top 25 genes and bottom 25 genes
#    based on log2 fold change. Save the top 25 genes and their log2 fold change to
#    ./result/pasilla.log2FC.top25.csv and those of bottom 25 genes to
#    ./result/pasilla.log2FC.bottom25.csv