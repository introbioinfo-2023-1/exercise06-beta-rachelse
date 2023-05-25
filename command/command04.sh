# command04.sh
# 1. Merge the GTF files from previous step with StringTie. 
#     Save the merged GTF file as ./result/merged.gtf. (Result file: merged.gtf)

stringtie --merge -p 8 -G ./data/s_cerevisiae.genes.gtf -o result/merged.gtf result/batch2.gtf result/chem2.gtf result/batch1.gtf result/chem1.gtf

# 2. Requantify transcript expression with StringTie. 
#    Save gtf files in data/ballgown directory. (No result file)

stringtie -e -B -G ./result/merged.gtf -o data/batch2_requant.gtf ./data/batch2.sorted.bam
stringtie -e -B -G ./result/merged.gtf -o data/chem2_requant.gtf ./data/chem2.sorted.bam
stringtie -e -B -G ./result/merged.gtf -o data/batch1_requant.gtf ./data/batch1.sorted.bam
stringtie -e -B -G ./result/merged.gtf -o data/chem1_requant.gtf ./data/chem1.sorted.bam

# 3. Find the FPKM and TPM of CYS3 gene in two TSV files and fill in the ./result/CYS3.csv ##TODO
# awk 'BEGIN {OFS=",";print"GeneName,SampleID,FPKM,TPM"}{f=substr(FILENAME,8,9)} $0~/CYS3/{print "CYS3",f,$8,$9}' result/{batch,chem}*.tsv > result/CYS3.csv

# 4. Get top 5 highly expressed genes (based on TPM) from ./result/SRR453567.tsv and ./result/SRR453570.tsv. ##TODO
#    Save the lines of top 5 genes to ./result/SRR453567.top5.tsv and ./result/SRR453570.top5.tsv.
#    (Result files: SRR453567.top5.tsv, SRR453570.top5.tsv)

# cat result/batch2.tsv | awk 'FNR!=1 && $2!="-"{print $0}' | sort -t $'\t' -k9rn | head -5 > result/batch2.top5.tsv
# cat result/chem2.tsv | awk 'FNR!=1 && $2!="-"{print $0}' | sort -t $'\t' -k9rn | head -5 > result/chem2.top5.tsv

# cat result/batch1.tsv | awk 'FNR!=1 && $2!="-"{print $0}' | sort -t $'\t' -k9rn | head -5 > result/batch1.top5.tsv
# cat result/chem1.tsv | awk 'FNR!=1 && $2!="-"{print $0}' | sort -t $'\t' -k9rn | head -5 > result/chem1.top5.tsv