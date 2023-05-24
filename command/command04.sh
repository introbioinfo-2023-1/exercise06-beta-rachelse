# command04.sh
# 1. Merge the GTF files from previous step with StringTie. 
#     Save the merged GTF file as ./result/merged.gtf. (Result file: merged.gtf)
stringtie --merge -p 8 -G ./data/s_cerevisiae.genes.gtf -o result/merged.gtf result/SRR453567.gtf result/SRR453570.gtf

# 2. Compare final transcripts with reference annotation with gffcompare. 
#     Save the result as ./result/merged.gtf.tmap. (Result file: gcp.merged.gtf.tmap, gcp.merged.gtf.refmap)
gffcompare -r ./data/s_cerevisiae.genes.gtf -o db/gcp/gcp -G result/merged.gtf

# 3. From result/gcp.merged.gtf.tmap extract lines that are not exact match of intron chain. 
#     Save the extracted lines into result/gcp_unmatched.txt. (Result file: gcp_unmatched.txt)
awk '$3!="="{print $0}' result/gcp.merged.gtf.tmap > result/gcp_unmatched.txt

# 4. 
#    4.1. Requantify transcript expression with StringTie. 
#    Save gtf files in data/ballgown directory. (No result file)

#    4.2. Run prepDE.py to convert the stringtie output into a format that can be used by DESeq2. 
#    Save the result files into ./result directory. (Result files: gene_count_matrix.csv, transcript_count_matrix.csv)

stringtie -e -B -G ./result/merged.gtf -o data/ballgown/SRR453567_requant.gtf ./data/SRR453567.sorted.bam
stringtie -e -B -G ./result/merged.gtf -o data/ballgown/SRR453570_requant.gtf ./data/SRR453570.sorted.bam
prepDE.py -i data/prep_deseq.txt -g result/gene_count_matrix.csv -t result/transcript_count_matrix.csv