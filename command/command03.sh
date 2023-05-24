# command03.sh
# 1. Download and extract the gzipped GTF file of S.cerevisiae from the link.
#    Convert the chromosome names of downloaded GTF file to have "chr" at the beginning.
#    e.g., 3R --> chr3R
#    Save the GTF file as ./data/s_cerevisiae.genes.gtf. (No result file)
wget -O ./data/saccharomyces_cerevisiae..genes.gtf.gz https://ftp.ensembl.org/pub/release-109/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf.gz
gunzip ./data/S_cerevisiae.genes.gtf.gz
awk '{if ($0~/^#/) {print $0} else {print "chr"$0}}' data/saccharomyces_cerevisiae.genes.gtf > data/s_cerevisiae.genes_.gtf 

# 2. Assemble and quantify transcripts for ./data/SRR453567.sorted.bam and ./data/SRR453570.sorted.bam with StringTie. 
#     Save the output files to the result directory. (Result files: SRR453567.tsv, SRR453567.gtf, SRR453570.tsv, SRR453570.gtf)

stringtie -p 8 -o ./result/SRR453567.gtf -A ./result/SRR453567.tsv ./data/SRR453567.sorted.bam -G ./data/s_cerevisiae.genes.gtf 
stringtie -p 8 -o ./result/SRR453570.gtf -A ./result/SRR453570.tsv ./data/SRR453570.sorted.bam -G ./data/s_cerevisiae.genes.gtf

# 3. Find the FPKM and TPM of PHO11 gene in two TSV files and fill in the ./result/pho11.csv
awk 'BEGIN {OFS=",";print"GeneName,SampleID,FPKM,TPM"}{f=substr(FILENAME,8,9)} $0~/PHO11/{print "PHO11",f,$8,$9}' result/SRR4535*.tsv > result/pho11.csv

# 5. Get top 5 highly expressed genes (based on TPM) from ./result/SRR453567.tsv and ./result/SRR453570.tsv.
#    Save the lines of top 5 genes to ./result/SRR453567.top5.tsv and ./result/SRR453570.top5.tsv.
#    (Result files: SRR453567.top5.tsv, SRR453570.top5.tsv)

cat result/SRR453567.tsv | awk 'FNR!=1 && $2!="-"{print $0}' | sort -t $'\t' -k9rn | head -5 > result/SRR453567.top5.tsv
cat result/SRR453570.tsv | awk 'FNR!=1 && $2!="-"{print $0}' | sort -t $'\t' -k9rn | head -5 > result/SRR453570.top5.tsv