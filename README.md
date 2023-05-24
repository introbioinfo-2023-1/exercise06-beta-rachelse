| Organization | Course | Exercise | Semester | Professor | T.A. | 
| --- | --- | --- | --- | --- | --- |
| Seoul National University | Introduction to Bioinformatics | Gene finding & RNA & Gene expression | Spring 2023 | Asst. Prof. M. Steinegger (martin.steinegger@snu.ac.kr) | Luna Jang (jse9512@snu.ac.kr) <br /> Rachel Kim (eunbelivable@snu.ac.kr)|

For submission, please follow the instructions in [submission](#submission).

# Exercise06: Gene finding & RNAseq analysis & Gene expression

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

In this exercise, we are going to use sub-commands `view`, `index`, `sort` and `coverage`.

### StringTie: [link for manual](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)
StringTie is a fast and efficient assembler of RNA-Seq alignments into potential transcripts.
- Input: BAM file with RNA-Seq read mappings (must be sorted)
- Output
   - GTF file containing the assembled transcripts
   - Gene abundances in tab-delimted format (TSV)
   - GTF file with requantified transcript expression

### gffcompare: [link for manual](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)
The program is used to compare, merge, annotate and estimate accuracy of one or more GFF files to compare with a reference annotation. It can compare and evaluate accuracy of RNA-seq transcript assemblers. Also, it merge duplicate transcripts from multiple GTF/GFF files as well as classify transcripts as they relate to reference transcripts provided in a annotation file.

### pyDEseq2: [link for manual](https://pydeseq2.readthedocs.io/en/latest/)
The pyDESeq2 package is python version of R package DEseq2. It is designed for normalization, visualization, and differential analysis of high dimensional count data. It makes use of empirical Bayes techniques to estimate priors for log fold change and dispersion, and to calculate posterior estimates for these quantities.

### EnhancedVolcano: [link for manual](https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html)
Volcano plots represent a useful way to visualise the results of differential expression analyses. EnhancedVolcano is a highly-configurable function that produces publication-ready volcano plots. EnhancedVolcano (Blighe, Rana, and Lewis 2018) will attempt to fit as many labels in the plot window as possible, thus avoiding ‘clogging’ up the plot with labels that could not otherwise have been read. Other functionality allows the user to identify up to 5 different types of attributes in the same plot space via colour, shape, size, encircling, and shade parameter configurations.

---
## command01.sh
1. Download the HISAT2 index of Saccharomyces cerevisiae genome **"sacCer3"** to `./db` and extract the `tar.gz` file.
   Remove the downloaded `sc3.tar.gz`. (No result file)
   - Link: https://cloud.biohpc.swmed.edu/index.php/s/Gsq4goLW4TDAz4E/download

2. We will use reads from two *S.cerevisiae* RNA-seq runs; SRR453567(batch) and SRR453570(chem).
   One sample is batch condition (glucose-excess) and the other is chem condition (glucose-limited).
   Download these FASTQ files in `./data` directory. (No result file)
   - https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch2_chrI_1.fastq 
   - https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch2_chrI_2.fastq 
   - https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem2_chrI_1.fastq 
   - https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem2_chrI_2.fastq 

   Output
      - `./data/SRR453567_1.fastq`
      - `./data/SRR453567_2.fastq`
      - `./data/SRR453570_1.fastq`
      - `./data/SRR453570_2.fastq`

3. Map the paired-end reads of two samples to **sc3** genome with **HISAT2**.
   Save the alignment result at `./data/SRR453567.sam` and `./data/SRR453570.sam`.
   (No result file)
   - Input
      - `./data/SRR453567_1.fastq`, `./data/SRR453567_2.fastq`
      - `./data/SRR453570_1.fastq`, `./data/SRR453570_1.fastq`
   - Output
      - `./data/SRR453567.sam`
      - `./data/SRR453570.sam`

## command02.sh

1. Convert the HISAT2-mapped SAM files to BAM format with `samtools view` command. (No result file)
   - Input: `./data/SRR453567.sam`, `./data/SRR453570.sam`
   - Output: `./data/SRR453567.bam`, `./data/SRR453570.bam`

2. Sort the BAM files from previous step with `samtools sort` command.
   Save the sorted BAM file as `*.sorted.bam`. (No result file)
   - Input: `./data/SRR453567.bam`, `./data/SRR453570.bam`
   - Output: `./data/SRR453567.sorted.bam`, `./data/SRR453570.sorted.bam`

3. Make index files for sorted BAM files from previous step with `samtools index` command. (No result file)
   - Input: `./data/SRR453567.sorted.bam`, `./data/SRR453570.sorted.bam`
   - Output: `./data/SRR453567.sorted.bam.bai`, `./data/SRR453570.sorted.bam.bai`

## command03.sh

1. Download and extract the gzipped GTF file of *S.cerevisiae* from the link.

   - Link: https://ftp.ensembl.org/pub/release-109/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf.gz
   **Convert the chromosome names** of downloaded GTF file to have **"chr"** at the beginning.
   > e.g., 3R --> chr3R

   Save the GTF file as `./data/s_cerevisiae.genes.gtf`. (No result file)

   - Output: `./data/s_cerevisiae.genes.gtf`

2. Assemble and quantify transcripts for `./data/SRR453567.sorted.bam` and `./data/SRR453570.sorted.bam` with **StringTie**. Save the output files(SRR453567.gtf, SRR453570.gtv) and quantification results(SRR453567.tsv, SRR453570.tsv) to the `result` directory. (Result files: **StringTie output**)
   - Input
      - BAM files
         - `./data/SRR453567.sorted.bam`
         - `./data/SRR453570.sorted.bam`
      - Guide GTF file: `./data/s_cerevisiae.genes.gtf`
   - Output
      - GTF files for assembled transcripts
         - `./result/SRR453567.gtf`
         - `./result/SRR453570.gtf`
      - TSV files for gene abundances
         - `./result/SRR453567.tsv`
         - `./result/SRR453570.tsv`

3. This is the information of *PHO11* gene.

   SGD ID  | Systemic name | Standard name    | Species           | Chromosome | Start   | End
   ------------|--------|---------|-------------------|------------|---------|--------
   S000000094 | YAR071W     | PhO11 | *S.cerevisiae* | chrI      | 225460 | 226863 

   Find the FPKM and TPM of *PHO11* gene in two TSV files and Save in the `./result/pho11.csv` as following.
   ``` pho11.csv
   GeneName,SampleID,FPKM,TPM
   PHO11,SRR453567,0.0,0.0
   PHO11,SRR453570,0.0,0.0
   ```
   (Result file: **pho11.csv**)

   - Input: `./result/SRR453567.tsv`, `./result/SRR453570.tsv`
   - Output: `./result/PHO11.csv`

4. Get top 5 highly expressed genes (based on TPM) from `./result/SRR453567.tsv` and `./result/SRR453570.tsv`.
   Save the lines of top 5 genes to `./result/SRR453567.top5.tsv` and `./result/SRR453570.top5.tsv`.

   (Result files: **SRR453567.top5.tsv**, **SRR453570.top5.tsv**)

   > - Only consider genes with gene names (Ignore transcripts of which gene names are "-")
   > - [Sort](https://man7.org/linux/man-pages/man1/sort.1.html) command will be helpful.
   > - `*_top5.tsv` should not include the header lines.

## command04.sh
To requantify transcript expression, we will merge the reconstructed transcriptomes from two samples and then quantify the merged transcriptome.
1. Merge the GTF files from previous step with **StringTie**. Save the merged GTF file as `./result/merged.gtf`. (Result file: **merged.gtf**)
   - Input: `./result/SRR453567.gtf`, `./result/SRR453570.gtf`
   - Output: `./result/merged.gtf`
  
   > - use `--merge` option of StringTie to merge the GTF files.

2. Compare final transcripts with reference annotation with **gffcompare**. Save the result as `./result/merged.gtf.tmap`. (Result file: **gcp.merged.gtf.tmap**, **gcp.merged.gtf.refmap**)
   - Input: `./result/merged.gtf`, `./data/s_cerevisiae.genes.gtf`
   - Output: `./result/gcp.merged.gtf.tmap`, `./result/gcp.merged.gtf.refmap`
  
   > - use `-o db/gcp/gcp` option of gffcompare. 
   > This will save .refmap and .tmap files into result directory while saving .loci, .tracking and .stats files into db/gcp directory.

3. From `result/gcp.merged.gtf.tmap` extract lines that are not exact match of intron chain. Save the extracted lines into `result/gcp_unmatched.txt`. (Result file: **gcp_unmatched.txt**)
   - Input: `./result/gcp.merged.gtf.tmap`
   - Output: `./result/gcp_unmatched.txt`

   > [Transcript classification codes](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml#transfrag-class-codes)
   > Exact match will be marked as **=** in class_code column of gffcompare output.

4. To prepare differential expression analysis, you should convert stringtie output into a format that can be used by DESeq2. Use `prepDE.py` to convert the stringtie output into a format that can be used by DESeq2. For preparation, we need two steps. First, requantify transcript expression and then convert the output into a format that can be used by DESeq2.
   
   4.1. Requantify transcript expression with **StringTie**. Save gtf files in `data/ballgown` directory. (No result file)
   - Input: `./result/merged.gtf` `./data/SRR453567.sorted.bam`, `./data/SRR453570.sorted.bam`
   - Output: `./data/ballgown/SRR453567_requant.gtf`, `./data/ballgown/SRR453570_requant.gtf`
   > Use `-Be` option of StringTie to requantify transcript expression.

   4.2. Run prepDE.py to convert the stringtie output into a format that can be used by DESeq2. Save the result files into `./result` directory. (Result files: **gene_count_matrix.csv**, **transcript_count_matrix.csv**)
   - Input: `./data/prep_deseq.txt`, `./data/ballgown/SRR453567_requant.gtf`, `./data/ballgown/SRR453570_requant.gtf`
   - Output: `./result/gene_count_matrix.csv`, `./result/transcript_count_matrix.csv`
  
   > copy and paste following command to command file
      ```sh
      prepDE.py -i data/prep_deseq.txt
      ```
   > move output files to result directory or use `-g` and `-t` options of prepDE.py to specify output directory

## command05.sh
In this stage, you will learn how to analyze Differential Expression (DE) with pyDESeq2.


Differential Expression (DE) analysis is one of the most common types of analyses
when you work with RNA-seq data. Usually, we use public R libraries for DE analysis,
such as [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) or
[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html).
Although the normalization steps which is done with these libraries are important
in the DE analyses, we are going to do a quick and easier (also not normalized)
version of DE analysis in this exercise.

1. Join the `./result/SRR453567.tsv` and `./result/SRR453570.tsv` with the `Gene ID` column.
   Only consider genes with gene names (Ignore transcripts of which gene names are "-").
   (No result file)
   - Input: `./result/SRR453567.tsv`, `./result/SRR453570.tsv`
   - Output: Joined tsv file

   > You may use `join` command or `join.awk` from `exercise03`.

2. [Fold change](https://en.wikipedia.org/wiki/Fold_change) is the ratio of the
   expression of particular gene in condition 1 and 2 and we commonly use the log2 (base-2 logarithm)
   value of fold change in DE analysis. To prevent division by zero errors,
   pseudocounts are added to the denominator and numerator when calculating the fold change.

   ```
   log2_fold_change = log2((expression_from_condition1 + pseudo_count) / (expression_from_condition2 + pseudo_count))
   ```

   From the joined tsv file, calculate the log2 fold change of all genes between
   wild-type sample and the *PHO11* RNAi sample and find the top 25 genes and bottom 25 genes
   based on log2 fold change. Save the top 25 genes and their log2 fold change to
   `./result/PHO11.log2FC.top25.csv` and those of bottom 25 genes to `./result/PHO11.log2FC.bottom25.csv`.
   (Result files: **PHO11.log2FC.top25.csv** and **PHO11.log2FC.bottom25.csv**)
   Write the **gene names as column 1** and the **log2 fold change as column 2** like this:
   ```
   gene1,3.5
   gene2,3.4
   ...
   ```
   For the calculation, use **TPM** values and **1** as psuedocount.

   > AWK has a `log()` function for calculating natual logarithm.

   > You may generate intermediate files.

   > You can use [sort](https://man7.org/linux/man-pages/man1/sort.1.html) command for
   sorting log2 fold change values.


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

Please submit these files as results:

- StringTie output (command05.sh - Step 2.)
   - **SRR453567.gtf**, **SRR453570.gtf**, **SRR453567.tsv**, **SRR453567.tsv**
- FPKM and TPM values of *PHO11* gene: **PHO11.csv** (command05.sh - Step 3.)
- Histograms for *PHO11* gene region (command05.sh - Step 4.)
   - **SRR453567.PHO11.coverage**, **SRR453570.PHO11.coverage**
- Sample ID with RNAi: **sample_with_RNAi.txt** (command05.sh - Step 4.)
- Top 5 highly expressed genes: **SRR453567.top5.tsv**, **SRR453570.top5.tsv** (command05.sh - Step 5.)
- Top and bottom 25 DE genes: **PHO11.log2FC.top25.csv** and **PHO11.log2FC.bottom25.csv** (command06.sh)
