# command01.sh
# 1. The link below is the HISAT2 index of Saccharomyces cerevisiae genome "sacCer3". 
#     Download the index to db/ directory and name as sc3.tar.gz.
#     Extract the tar.gz file and remove the downloaded sc3.tar.gz. (No result file)
#    - Link: https://cloud.biohpc.swmed.edu/index.php/s/Gsq4goLW4TDAz4E/download
wget -O ~/exercise06/db/sacCer3.tar.gz  https://cloud.biohpc.swmed.edu/index.php/s/Gsq4goLW4TDAz4E/download
tar -zxvf ./db/sacCer3.tar.gz -C ./db
rm ./db/sacCer3.tar.gz

# 2. We will use paired-end reads from four S.cerevisiae RNA-seq runs; batch(SRR453566, SRR453567) and chem(SRR453569, SRR453570).
#    Two samples are batch condition (glucose-excess) and the other are chem condition (glucose-limited). Download these FASTQ files in ./data directory. 
#    Name each file as batch1_1.fastq, batch1_2.fastq, batch2_1.fastq, batch2_2.fastq, chem1_1.fastq, chem1_2.fastq, chem2_1.fastq, and chem2_2.fastq. (No result file)
#    - batch1_1: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch1_chrI_1.fastq
#    - batch1_2: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch1_chrI_2.fastq 
#    - batch2_1: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch2_chrI_1.fastq 
#    - batch2_2: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem1_chrI_1.fastq 
#    - chem1_1: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem1_chrI_2.fastq 
#    - chem1_2: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch2_chrI_2.fastq 
#    - chem2_1: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem2_chrI_1.fastq 
#    - chem2_2: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem2_chrI_2.fastq 
wget -O data/batch2_1.fastq https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch2_chrI_1.fastq
wget -O data/batch2_2.fastq https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch2_chrI_2.fastq
wget -O data/chem2_1.fastq https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem2_chrI_1.fastq
wget -O data/chem2_2.fastq https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem2_chrI_2.fastq

wget -O data/batch1_1.fastq https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch1_chrI_1.fastq
wget -O data/batch1_2.fastq https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch1_chrI_2.fastq 
wget -O data/chem1_1.fastq https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem1_chrI_1.fastq 
wget -O data/chem1_2.fastq https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem1_chrI_2.fastq 

# 3. Map the paired-end reads of four samples to sc3 genome with HISAT2.
#    Save the alignment result as batch1.sam, batch2.sam, chem1.sam, and chem2.sam in data directory. (No result file)

hisat2 -x ./db/sc3/genome -1 ./data/batch2_1.fastq -2 ./data/batch2_2.fastq -S ./data/batch2.sam
hisat2 -x ./db/sc3/genome -1 ./data/chem2_1.fastq -2 ./data/chem2_2.fastq -S ./data/chem2.sam

hisat2 -x ./db/sc3/genome -1 ./data/batch1_1.fastq -2 ./data/batch1_2.fastq -S ./data/batch1.sam
hisat2 -x ./db/sc3/genome -1 ./data/chem1_1.fastq -2 ./data/chem1_2.fastq -S ./data/chem1.sam