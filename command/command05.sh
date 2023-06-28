# command05.sh
# 1. To prepare differential expression analysis, you should convert stringtie output into a format that can be used by DESeq2. 
#     Fortunately, stringtie provides code to generate input file of DESeq2. Use prepDE.py to convert the stringtie output into a format that can be used by DESeq2. 
#     Save the result files into ./result directory. (Result files: gene_count_matrix.csv, transcript_count_matrix.csv)
prepDE.py -i source/prep_deseq.txt -g result/gene_count_matrix.csv -t result/transcript_count_matrix.csv

# 2. With the output file of previous step, run pyDEseq2 using code ./source/run_pydeseq.py. 
#     To run the tool properly, you should specify options. Save the result file as deseq2_results.csv. (Result file: deseq2_results.csv)
python source/run_pydeseq.py -i result/gene_count_matrix.csv -o result/deseq2_results.csv -f condition -s 
python source/run_pydeseq.py -i result/gene_count_matrix.csv -o result/deseq2_results2.csv -f condition
# 3. From the result of pyDEseq2 (./result/deseq2_results.csv), save the number of overexpressed genes.
#    Save the result file as overexpressed_genes.txt. (Result file: overexpressed_genes.txt)
cat result/deseq2_results.csv | awk -F, 'FNR!= 1&& $3>0{print $1,$7}' | wc -l > result/overexpressed_genes.txt

# 4. Using the result file of pyDEseq2, make volcano plot. 
#     Save the plots colored by adjusted p-value 0.01 and 0.05 to each volcano_01.png and volcano_05.png. 
#     (Result file: volcano_01.png, volcano_05.png)

python source/make_volcano.py -i result/deseq2_results.csv -o result/volcano_05.png -c -x log2FoldChange -y pvalue
python source/make_volcano.py -i result/deseq2_results.csv -o result/volcano_01.png -c -x log2FoldChange -y pvalue -t 0.01

python source/make_volcano.py -i result/deseq2_results2.csv -o result/volcano_05_2.png -c -x log2FoldChange -y pvalue
python source/make_volcano.py -i result/deseq2_results2.csv -o result/volcano_01_2.png -c -x log2FoldChange -y pvalue -t 0.01

python source/make_volcano.py -i result/deseq2_results.csv -o result/volcano_05_3.png -c -x log2FoldChange -y padj
python source/make_volcano.py -i result/deseq2_results.csv -o result/volcano_01_3.png -c -x log2FoldChange -y padj -t 0.01

python source/make_volcano.py -i result/deseq2_results2.csv -o result/volcano_05_4.png -c -x log2FoldChange -y padj
python source/make_volcano.py -i result/deseq2_results2.csv -o result/volcano_01_4.png -c -x log2FoldChange -y padj -t 0.01