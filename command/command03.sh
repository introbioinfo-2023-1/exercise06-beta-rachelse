# command03.sh
# 1. Download the HISAT2 index of Drossophila genome **"dm6"** to `./db` and extract the `tar.gz` file.
#    After extraction, `make_dm6.sh` in the `./db/dm6` folder and remove the downloaded `dm6.tar.gz`.

# 2. Download these FASTQ files in ./data directory and change the file extensions
#    from .fastqsanger to .fastq
#    https://zenodo.org/record/4541751/files/GSM461177_1_subsampled.fastqsanger
#    https://zenodo.org/record/4541751/files/GSM461177_2_subsampled.fastqsanger
#    https://zenodo.org/record/4541751/files/GSM461180_1_subsampled.fastqsanger
#    https://zenodo.org/record/4541751/files/GSM461180_2_subsampled.fastqsanger

# 3. Map the paired-end reads of two samples to dm6 genome with HISAT2.
#    Save the alignment result at ./data/GSM461177.sam and ./data/GSM461180.sam