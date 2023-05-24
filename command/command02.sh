# command02.sh
# 1. Convert the HISAT2-mapped SAM files to BAM format with `samtools view` command.
samtools view -Sb ./data/SRR453567.sam > ./data/SRR453567.bam
samtools view -Sb ./data/SRR453570.sam > ./data/SRR453570.bam
# 2. Sort the BAM files from previous step with `samtools sort` command.
samtools sort ./data/SRR453567.bam -o ./data/SRR453567.sorted.bam
samtools sort ./data/SRR453570.bam -o ./data/SRR453570.sorted.bam
# 3. Make index files for sorted BAM files from previous step with `samtools index` command.
samtools index ./data/SRR453567.sorted.bam
samtools index ./data/SRR453570.sorted.bam