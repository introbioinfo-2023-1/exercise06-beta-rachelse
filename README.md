| Organization              | Course                         | Exercise         | Semester    | Professor                                               | T.A.                                                                     |
| ------------------------- | ------------------------------ | ---------------- | ----------- | ------------------------------------------------------- | ------------------------------------------------------------------------ |
| Seoul National University | Introduction to Bioinformatics | RNA-seq analysis | Spring 2023 | Asst. Prof. M. Steinegger (martin.steinegger@snu.ac.kr) | Luna Jang (jse9512@snu.ac.kr) <br /> Rachel Kim (eunbelivable@snu.ac.kr) |

- For submission, please follow the instructions in [submission](#submission).

# Exercise06: RNA-seq analysis pipeline

## About the exercise

In this exercise, you will practice RNA-seq pipeline with Saccharomyces cerevisiae data. 
The pipeline will follow mapping RNA-seq reads to reference genome, assembling transcripts and then quantifying gene expression. In addition, you will learn how to identify differentially expressed gene and visualize the results.

## Tools for the exercise

### HISAT2 : [link for manual](http://daehwankimlab.github.io/hisat2/manual/)
HISAT2 is an alignment program for mapping NGS reads (whole-genome, transcriptome, and exome sequencing data) against the general human population (as well as against a single reference genome).
- Input
   - hisat2 index for the reference genome
   - Sequencing read files (usually paired-end; _1.fastq, _2.fastq)
- Output
   - Alignment in SAM format

### SAMtools: [link for manual](http://www.htslib.org/doc/samtools.html)
Samtools is a set of utilities that manipulate alignments in the SAM (Sequence Alignment/Map), BAM, and CRAM formats. It converts between the formats, does sorting, merging and indexing, and can retrieve reads in any regions swiftly.

In this exercise, we are going to use sub-commands `view`, `index`, and `sort`.

### StringTie: [link for manual](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)
StringTie is a fast and efficient assembler of RNA-Seq alignments into potential transcripts.
- Input: BAM file with RNA-Seq read mappings (must be sorted)
- Output
   - GTF files containing the assembled transcripts
   - GTF files containing requantified transcript expression with merged assemblies from multiple samples

### pyDEseq2: [link for manual](https://pydeseq2.readthedocs.io/en/latest/)
The pyDESeq2 package is python version of R package DEseq2. It is designed for normalization, visualization, and differential analysis of high dimensional count data. It makes use of empirical Bayes techniques to estimate priors for log fold change and dispersion, and to calculate posterior estimates for these quantities.

---
## command01.sh
1. The link below is the HISAT2 index of Saccharomyces cerevisiae genome **"sacCer3"**. Download the index to `db/` directory and name as `sc3.tar.gz`.
Extract the `tar.gz` file and remove the downloaded `sc3.tar.gz`. (No result file)
   - Link: https://cloud.biohpc.swmed.edu/index.php/s/Gsq4goLW4TDAz4E/download

2. We will use paired-end reads from four *S.cerevisiae* RNA-seq runs; batch(SRR453566, SRR453567) and chem(SRR453569, SRR453570).
   Two samples are batch condition (glucose-excess) and the other are chem condition (glucose-limited). Download these FASTQ files in `./data` directory. 
   Name each file as `batch1_1.fastq`, `batch1_2.fastq`, `batch2_1.fastq`, `batch2_2.fastq`, `chem1_1.fastq`, `chem1_2.fastq`, `chem2_1.fastq`, and `chem2_2.fastq`. (No result file)
   - batch1_1: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch1_chrI_1.fastq
   - batch1_2: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch1_chrI_2.fastq 
   - batch2_1: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch2_chrI_1.fastq 
   - batch2_2: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem1_chrI_1.fastq 
   - chem1_1: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem1_chrI_2.fastq 
   - chem1_2: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch2_chrI_2.fastq 
   - chem2_1: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem2_chrI_1.fastq 
   - chem2_2: https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem2_chrI_2.fastq 

   Output
     - `./data/batch1_1.fastq`, `./data/batch1_2.fastq`
     - `./data/batch2_1.fastq`, `./data/batch2_2.fastq`
     - `./data/chem1_1.fastq`, `./data/chem1_2.fastq`
     - `./data/chem2_1.fastq`, `./data/chem2_2.fastq`

3. Map the paired-end reads of four samples to **sc3** genome with **HISAT2**.
   Save the alignment result as `batch1.sam`, `batch2.sam`, `chem1.sam`, and `chem2.sam` in **data** directory.
   (No result file)
   - Input
      - `./data/batch1_1.fastq`, `./data/batch1_2.fastq`
      - `./data/batch2_1.fastq`, `./data/batch2_2.fastq`
      - `./data/chem1_1.fastq`, `./data/chem1_2.fastq`
      - `./data/chem2_1.fastq`, `./data/chem2_2.fastq`
   - Output
      - `./data/batch1.sam`
      - `./data/batch2.sam`
      - `./data/chem1.sam`
      - `./data/chem2.sam`

## command02.sh

1. Convert the HISAT2-mapped SAM files to BAM format with `samtools view` command. (No result file)
   - Input: `./data/batch1.sam`, `./data/batch2.sam`, `./data/chem1.sam`, `./data/chem2.sam`
   - Output: `./data/batch1.bam`, `./data/batch2.bam`, `./data/chem1.bam`, `./data/chem2.bam`

2. Sort the BAM files from previous step with `samtools sort` command.
   Save the sorted BAM file as `*.sorted.bam`. (No result file)
   - Input: `./data/batch1.bam`, `./data/batch2.bam`, `./data/chem1.bam`, `./data/chem2.bam`
   - Output: `./data/batch1.sorted.bam`, `./data/batch2.sorted.bam`, `./data/chem1.sorted.bam`, `./data/chem2.sorted.bam`

3. Make index files for sorted BAM files from previous step with `samtools index` command. (No result file)
   - Input: `./data/batch1.sorted.bam`, `./data/batch2.sorted.bam`, `./data/chem1.sorted.bam`, `./data/chem2.sorted.bam`
   - Output: `./data/batch1.sorted.bam.bai`, `./data/batch2.sorted.bam.bai`, `./data/chem1.sorted.bam.bai`, `./data/chem2.sorted.bam.bai`

## command03.sh

1. Download and extract the gzipped GTF file of *S.cerevisiae* from the link.

   - Link: https://ftp.ensembl.org/pub/release-109/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf.gz
   **Convert the chromosome names** of downloaded GTF file to have **"chr"** at the beginning.
   > e.g., 3R --> chr3R

   Save the GTF file as `./data/s_cerevisiae.genes.gtf`. (No result file)

   - Output: `./data/s_cerevisiae.genes.gtf`

2. Assemble transcripts for `./data/batch1.sorted.bam`, `./data/batch2.sorted.bam`, `./data/chem1.sorted.bam`, and `./data/chem2.sorted.bam` with **StringTie**. Save the output files(batch1.gtf, batch2.gtf, chem1.gtf, chem2.gtf) to the `result` directory. (Result files: **batch1.gtf**, **batch2.gtf**, **chem1.gtf**, **chem2.gtf**)
   - Input
      - BAM files
         - `./data/batch1.sorted.bam`
         - `./data/batch2.sorted.bam`
         - `./data/chem1.sorted.bam`
         - `./data/chem2.sorted.bam`
      - Guide GTF file: `./data/s_cerevisiae.genes.gtf`
   - Output
         - `./result/batch1.gtf`, `./result/batch2.gtf`, `./result/chem1.gtf`, `./result/chem2.gtf`
   > Use `-G` option to provide the guide GTF file.



## command04.sh
To requantify transcript expression, we will merge the reconstructed transcriptomes from four samples and then quantify the merged transcriptome.
1. Merge the GTF files from previous step with **StringTie**. Save the merged GTF file as `./result/merged.gtf`. (Result file: **merged.gtf**)
   - Input: 
     - `./result/batch1.gtf`, `./result/batch2.gtf`, `./result/chem1.gtf`, `./result/chem2.gtf`
     - Guide GTF file: `./data/s_cerevisiae.genes.gtf`
   - Output: `./result/merged.gtf`
  
   > - use `--merge` option of StringTie to merge the GTF files.

2. Requantify transcript expression with **StringTie**. Save gtf files in `data/ballgown` directory. (No result file)
   - Input: `./result/merged.gtf` `./data/batch1.sorted.bam`, `./data/batch2.sorted.bam`, `./data/chem1.sorted.bam`, `./data/chem2.sorted.bam`
   - Output: `./data/batch1_requant.gtf`, `./data/batch2_requant.gtf`, `./data/chem1_requant.gtf`, `./data/chem2_requant.gtf`
   > Use `-B` and `-e` options of StringTie to requantify transcript expression.
   > Use `merged.gtf` as a guide GTF file.   

3. This is the information of *CYS3* gene.

   | Systemic name | Standard name | Species        | Chromosome | Start  | End    |
   | ------------- | ------------- | -------------- | ---------- | ------ | ------ |
   | YAL012W       | CYS3          | *S.cerevisiae* | chrI       | 130799 | 131983 |

   Find the FPKM and TPM of *CYS3* gene in four gtf files and Save in the `./result/CYS3.csv` as following. (Result file: **CYS3.csv**)
   ``` CYS3.csv
   GeneName,SampleID,FPKM,TPM
   CYS3,batch1,0.0,0.0
   CYS3,batch2,0.0,0.0
   CYS3,chem1,0.0,0.0
   CYS3,chem2,0.0,0.0
   ```
   > You should only consider **transcript** features.

   - Input: `./result/batch1_requant.gtf`, `./result/batch2_requant.gtf`, `./result/chem1_requant.gtf`, `./result/chem2_requant.gtf`
   - Output: `./result/cys3.csv`

4. Get top 5 highly expressed **transcripts** (based on TPM) from requantified gtf files.
   Save gene_id, gene_name and TPM of top 5 transcripts to `./result/batch1.top5.tsv` and `./result/batch2.top5.tsv`, `./result/chem1.top5.tsv`, `./result/chem2.top5.tsv`. (Result files: **batch1.top5.tsv**, **batch2.top5.tsv**, **chem1.top5.tsv**, **chem2.top5.tsv**)

   > - Only consider genes with gene names (Ignore transcripts without gene names)
   > - Result should be sorted by TPM in descending order.
   ```e.g.
   MSTRG.1  ABCD  12345.5000
   MSTRG.2  BCDE  12344.4000
   MSTRG.3  CDEF  12343.3000
   MSTRG.4  DEFG  12342.2000
   MSTRG.5  EFGH  12341.1000
   ```
   
   - Input: `./data/ballgown/batch1_requant.gtf`, `./data/ballgown/batch2_requant.gtf`, `./data/ballgown/chem1_requant.gtf`, `./data/ballgown/chem2_requant.gtf`
   - Output: `./result/batch1.top5.tsv`, `./result/batch2.top5.tsv`, `./result/chem1.top5.tsv`, `./result/chem2.top5.tsv`


## command05.sh
In this stage, you will learn how to analyze Differential Expression (DE) with pyDESeq2. After the analysis, you will visualize the expression result with volcano plot. Since the analysis and visualization require knowledge of R or Python, we will provide the sourcecodes(run_pydeseq.py, make_volcano.py). Instead of writing the code, you should provide arguments to run the sourcecodes.
Using the data of previous step, we will compare the gene expressions between two conditions(batch and chem).

:bangbang: Note that due to small replications, the result of DE analysis may not be reliable and the plot looks poor. However, you don't have to care about this since the purpose of this step is to learn briefly how to analyze DE with pyDESeq2 and visualize the result.
:bangbang: This procedure is abbreviated excluding run DESeq2 and plot by yourself. If you are interested in the detailed procedure, please read over the sourcecodes in `./source` directory and result files.

1. To prepare differential expression analysis, you should convert stringtie output into a format that can be used by DESeq2. Fortunately, stringtie provides code to generate input file of DESeq2. Use `prepDE.py` to convert the stringtie output into a format that can be used by DESeq2. Save the result files into `./result` directory. (Result files: **gene_count_matrix.csv**, **transcript_count_matrix.csv**)
   - Input: `./source/prep_deseq.txt`
   - Output: `./result/gene_count_matrix.csv`, `./result/transcript_count_matrix.csv`
  
   > **copy and paste** following command to command file
      ```sh
      prepDE.py -i source/prep_deseq.txt
      ```
   > move output files to result directory or use `-g` and `-t` options of prepDE.py to specify output directory

2. With the output file of previous step, run pyDEseq2 using code `./source/run_pydeseq.py`. To run the tool properly, you should specify options. Save the result file as **deseq2_results.csv**. (Result file: **deseq2_results.csv**)
   - Input: `./result/gene_count_matrix.csv`
   - Output: `./result/deseq2_results.csv`
   - command usage:
      ```sh
      python source/run_pydeseq.py -i <inputfile> -o <output_file> -f <factor>
      ```
      - `-i`, `-o`, and `-f` are required.
      - `-f`: You should specify **factor to compare**. You can find appropriate factor among column names of `source/phenodata.csv` file. This is case sensitive.
      - `-s`: is optional. You can enhance False Positve Rate and visualization with this option.
      - `-h`: For details of options and usage, use `-h` option.
  
3. From the result of pyDEseq2 (`./result/deseq2_results.csv`), save the **number** of overexpressed genes.
   Save the result file as **overexpressed_genes.txt**. (Result file: **overexpressed_genes.txt**)
  > - No threshold for overexpression. Just extract all the gene ids of overexpressed genes.

4. Using the result file of pyDEseq2, make volcano plot. Save the plot as **volcano.png**. (Result file: **volcano.png**)
   - Input: `./result/deseq2_results.csv`
   - Output: `./result/volcano.png`
   - command usage:
      ```sh
      python source/make_volcano.py -i <inputfile> -o <output_file> -x <x_axis> -y <y_axis>
      ```
      - `-i`, `-o`, `-x`, and `-y` are mandatory.
      - `-x`: You should specify **x-axis** of plot. You can find appropriate factor among column names of the input file. This is case sensitive and please use exactly same name as input file.
      - `-y`: You should specify **y-axis** of plot. You can find appropriate factor among column names of the input file. This is case sensitive and please use exactly same name as input file.
      - `-c`: is optional. You can color significant genes(based on adjusted p-value) in red with this option.
      - `-t`: is optional. You can modify threshold to color significant genes based on adjusted p-value. Default is 0.05.
      - `-h`: For details of options and usage, use `-h` option.

   > Refer to the plot on lecture slide. (19 RNA-seq differential gene expression - page 27)


---

*Reference*
   - Nookaew I, Papini M, Pornputtapong N, Scalcinati G, Fagerberg L, Uhlén M, Nielsen J. A comprehensive comparison of RNA-Seq-based transcriptome analysis from reads to differential gene expression and cross-comparison with microarrays: a case study in Saccharomyces cerevisiae. Nucleic Acids Res. 2012 Nov 1;40(20):10084-97. doi: 10.1093/nar/gks804. Epub 2012 Sep 10. PMID: 22965124; PMCID: PMC3488244.

---
## Submission

To submit your result, follow these steps:

- Step 1. Clone this template repository to your working directory and execute "setup.sh"
- Step 2. Fill in the command used in the command0X.sh in the "command" directory. The commands should generate the result of step 3. The result can either be printed to the terminal or written to a file.
- Step 3. Save the result files for each command.
- Step 4. Add edited files to git and commit
   ```sh
   git add .
   git commit -m "COMMIT MESSAGE"
   ```
- Step 5. Submit your answers by pushing the changes.
   ```sh
   git push origin master
   ```