# command03.sh
# 1. Download and extract the gzipped GTF file of S.cerevisiae from the link.
#    Convert the chromosome names of downloaded GTF file to have "chr" at the beginning.
#    e.g., 3R --> chr3R
#    Save the GTF file as ./data/s_cerevisiae.genes.gtf. (No result file)
wget -O ./data/saccharomyces_cerevisiae.genes.gtf.gz https://ftp.ensembl.org/pub/release-109/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf.gz
gunzip ./data/saccharomyces_cerevisiae.genes.gtf.gz
awk '{if ($0~/^#/) {print $0} else {print "chr"$0}}' data/saccharomyces_cerevisiae.genes.gtf > data/s_cerevisiae.genes.gtf 
rm data/saccharomyces_cerevisiae.genes.gtf

# 2. Assemble transcripts for ./data/batch1.sorted.bam, ./data/batch2.sorted.bam, ./data/chem1.sorted.bam, and ./data/chem2.sorted.bam with StringTie. 
#     Save the output files(batch1.gtf, batch2.gtf, chem1.gtf, chem2.gtf) to the result directory. 
#     (Result files: batch1.gtf, batch2.gtf, chem1.gtf, chem2.gtf)

stringtie -p 8 -o ./result/batch2.gtf ./data/batch2.sorted.bam -G ./data/s_cerevisiae.genes.gtf 
stringtie -p 8 -o ./result/chem2.gtf ./data/chem2.sorted.bam -G ./data/s_cerevisiae.genes.gtf

stringtie -p 8 -o ./result/batch1.gtf ./data/batch1.sorted.bam -G ./data/s_cerevisiae.genes.gtf 
stringtie -p 8 -o ./result/chem1.gtf ./data/chem1.sorted.bam -G ./data/s_cerevisiae.genes.gtf
