# command01.sh
# 1. Download the HISAT2 index of Saccharomyces cerevisiae genome **"sacCer3"** to `./db` and extract the `tar.gz` file.
#    Remove the downloaded `sc3.tar.gz`. (No result file)
#    - Link: https://cloud.biohpc.swmed.edu/index.php/s/Gsq4goLW4TDAz4E/download
wget -O ~/exercise06/db/sacCer3.tar.gz  https://cloud.biohpc.swmed.edu/index.php/s/Gsq4goLW4TDAz4E/download
tar -zxvf ./db/sacCer3.tar.gz -C ./db
./db/sc3/make_sc3.sh
mkdir ./db/sc3/genome
mv *genome* ./db/sc3/genome

# 2. We will use reads from two *S.cerevisiae* RNA-seq runs; SRR453567(batch) and SRR453570(chem).
#    One sample is batch condition (glucose-excess) and the other is chem condition (glucose-limited).
#    Download these FASTQ files in `./data` directory. (No result file)
#    - https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch2_chrI_1.fastq 
#    - https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch2_chrI_2.fastq 
#    - https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem2_chrI_1.fastq 
#    - https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem2_chrI_2.fastq 
wget -O data/SRR453567_1.fastq https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch2_chrI_1.fastq
wget -O data/SRR453567_2.fastq https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch2_chrI_2.fastq
wget -O data/SRR453570_1.fastq https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem2_chrI_1.fastq
wget -O data/SRR453570_2.fastq https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem2_chrI_2.fastq

# 3. Map the paired-end reads of two samples to sc3 genome with HISAT2.
#    Save the alignment result at ./data/SRR453567.sam and ./data/SRR453570.sam

hisat2 -x ./db/sc3/genome -1 ./data/SRR453567_1.fastq -2 ./data/SRR453567_2.fastq -S ./data/SRR453567.sam
hisat2 -x ./db/sc3/genome -1 ./data/SRR453570_1.fastq -2 ./data/SRR453570_2.fastq -S ./data/SRR453570.sam
