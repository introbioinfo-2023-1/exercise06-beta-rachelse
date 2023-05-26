# command04.sh
# 1. Merge the GTF files from previous step with StringTie. 
#     Save the merged GTF file as ./result/merged.gtf. (Result file: merged.gtf)

stringtie --merge -p 8 -G ./data/s_cerevisiae.genes.gtf -o result/merged.gtf result/batch2.gtf result/chem2.gtf result/batch1.gtf result/chem1.gtf

# 2. Requantify transcript expression with StringTie. 
#    Save gtf files in data/ballgown directory. (No result file)

stringtie -e -B -G ./result/merged.gtf -o data/ballgown/batch2_requant.gtf ./data/batch2.sorted.bam
stringtie -e -B -G ./result/merged.gtf -o data/ballgown/chem2_requant.gtf ./data/chem2.sorted.bam
stringtie -e -B -G ./result/merged.gtf -o data/ballgown/batch1_requant.gtf ./data/batch1.sorted.bam
stringtie -e -B -G ./result/merged.gtf -o data/ballgown/chem1_requant.gtf ./data/chem1.sorted.bam

# 3. Find the FPKM and TPM of CYS3 gene in four quantified gtf files and fill in the ./result/CYS3.csv file.

awk -F"\t" '{if (FILENAME~/batch/) {f=substr(FILENAME,15,6)} else {f=substr(FILENAME,15,5)}} $0~/CYS3/ && $3=="transcript"{print "CYS3",f,$9}' data/ballgown/{batch,chem}*.gtf |\
 tr -d "\";" | awk 'BEGIN {OFS=",";print "GeneName,SampleID,FPKM,TPM"}{print $1,$2,$12,$14}' > result/CYS3.csv

# 4. Get top 5 highly expressed transcripts (based on TPM) from requantified gtf files.
#    Save gene_id, gene_name and TPM of top 5 transcripts to ./result/batch1.top5.tsv and ./result/batch2.top5.tsv, ./result/chem1.top5.tsv, ./result/chem2.top5.tsv. 
#    (Result files: batch1.top5.tsv, batch2.top5.tsv, chem1.top5.tsv, chem2.top5.tsv)

awk -F"\t" '$3=="transcript"{print $9}' data/ballgown/batch1_requant.gtf | tr -d "\";" |\
    awk 'BEGIN{OFS="\t"} $0~/gene_name/ {print $2,$6,$12}' | sort -k3nr | head -5 > result/batch1.top5.tsv
awk -F"\t" '$3=="transcript"{print $9}' data/ballgown/batch2_requant.gtf | tr -d "\";" |\
    awk 'BEGIN{OFS="\t"} $0~/gene_name/ {print $2,$6,$12}' | sort -k3nr | head -5 > result/batch2.top5.tsv
awk -F"\t" '$3=="transcript"{print $9}' data/ballgown/chem1_requant.gtf | tr -d "\";" |\
    awk 'BEGIN{OFS="\t"} $0~/gene_name/ {print $2,$6,$12}' | sort -k3nr | head -5 > result/chem1.top5.tsv
awk -F"\t" '$3=="transcript"{print $9}' data/ballgown/chem2_requant.gtf | tr -d "\";" |\
    awk 'BEGIN{OFS="\t"} $0~/gene_name/ {print $2,$6,$12}' | sort -k3nr | head -5 > result/chem2.top5.tsv