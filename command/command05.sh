# command05.sh
# 1. Run prepDE.py to convert the stringtie output into a format that can be used by DESeq2. 
#     Save the result files into ./result directory. (Result files: gene_count_matrix.csv, transcript_count_matrix.csv)
prepDE.py -i data/prep_deseq.txt -g result/gene_count_matrix.csv -t result/transcript_count_matrix.csv

# 2.
python run_pydeseq.py -i result/gene_count_matrix.csv -o result/deseq2_results.csv -f condition -s 

# 3. 
# python make_plot.py -i result/deseq2_results.csv -o result/volcano.png -c 1 