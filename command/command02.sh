# command02.sh
# 1. Convert the HISAT2-mapped SAM files to BAM format with samtools view command. (No result file)
samtools view -Sb ./data/batch2.sam > ./data/batch2.bam
samtools view -Sb ./data/chem2.sam > ./data/chem2.bam

samtools view -Sb ./data/batch1.sam > ./data/batch1.bam
samtools view -Sb ./data/chem1.sam > ./data/chem1.bam

# 2. Sort the BAM files from previous step with samtools sort command.
#    Save the sorted BAM file as *.sorted.bam. (No result file)
samtools sort ./data/batch2.bam -o ./data/batch2.sorted.bam
samtools sort ./data/chem2.bam -o ./data/chem2.sorted.bam

samtools sort ./data/batch1.bam -o ./data/batch1.sorted.bam
samtools sort ./data/chem1.bam -o ./data/chem1.sorted.bam

# 3. Make index files for sorted BAM files from previous step with samtools index command. (No result file)
samtools index ./data/batch2.sorted.bam
samtools index ./data/chem2.sorted.bam

samtools index ./data/batch1.sorted.bam
samtools index ./data/chem1.sorted.bam