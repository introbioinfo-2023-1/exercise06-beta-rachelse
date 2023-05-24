

## About the exercise
- This exercise is based on McDonnell Genome Institute's material

Differential expression analysis is used to identify differences in the transcriptome (gene expression) across a cohort of samples. Often, it will be used to define the differences between multiple biological conditions (e.g. drug treated vs. untreated samples). There are many, many tools available to perform this type of analysis. In this course we will rely on a popular Bioconductor package DEseq2. We will then make various visualizations to help interpret our results.

## Dataset
- We start this exercise with given transcript quatification data because we already practiced how to prepare them in the last exercise.

For this analysis we will use the RNAseq data obtained from the EBI Expression Atlas (GXA). Specifically data set E-GEOD-50760 which corresponds to PMID: 25049118. This data consists of 54 samples from 18 individuals. Each individual has a primary colorectal cancer sample, a metastatic liver sample, and a normal sample of the surrounding colonic epithilium. The quantification data required to run differential expression analysis using DEseq2 are raw readcounts for either genes or transcripts. We will use the output from HTseq as a starting point.

## Rscripts
- Skeleton codes for each step are provided in Rscripts directory. 
- Some scripts have specified spots that you have to fill in to make them work.
- Understand the context and search what should be inserted in the spots.
- The working progress will be stored in DEseq2.RData, so the variables are shared between the scripts.
- You have to upload the scripts after you fill the blanks.

---

## Step 1.init_and_load_data.sh
- Result file: None
- No blanks in the script

This step is for downloading and reformatting data so that DEseq2 can handle it.

## Step 2. diff_exp_analysis.sh
Run the function DEseq() on the DEseq2 data set object that we prepared in step 1.

The function performs the following:
   1. Estimation of size factors.
   2. Estimation of dispersion.
   3. Negative Binomial GLM fitting and Wald statistic. 

- Result file: None
- Check *init_and_load_data.sh* and *extract_result.sh*, and fill the blanks with an appropriate variable name. 


## Step 3. extract_results.sh
Our dataset was prepared from three distinct tissue types. Let's extract the DEG result between two tissue types (primary colorectal cancer & normal-looking surrounding colonic epithelium) 

- Result file: None.
- Search the usage of *results*, a built-in method of DEseq2, and fill the blanks.

## 4. MAplot.sh
Within computational biology, an MA plot is an application of a Bland–Altman plot for visual representation of genomic data. The plot visualizes the differences between measurements taken in two samples, by transforming the data onto M (log ratio) and A (mean average) scales, then plotting these values.
- Result file
   1. *MA-plot.jpg* will be generated. Move it to result directory and upload.
   2. **Briefly** explain the meaning of MA-plot and interpret the result in *MA-plot.txt* in the results directory.
- Search the built-in method of DEseq2 for MA-plot and fill the blanks


## 5. sort_and_plotCounts.sh
Find a statistically significant gene and plot gene counts measured in various tissue types and patients.
- Result file: 
   1. *count-plot.jpg* will be generated. Move it to result directory and upload.
   2. **Briefly** write your interpetation about the plot in *plotCounts.txt* in result directory.
- Fill the blanks in the script
   - Hint
      1. Search how to index a table in R
      2. Search the usage of plotCounts

## 6. volano_plot.sh
A volcano plot is a type of scatter plot that shows statistical significance (P value) versus magnitude of change (fold change). It enables quick visual identification of genes with large fold changes that are also statistically significant. These may be the most biologically significant genes.
- Result file: 
   1. *volcano_plot.pdf* will be generated. Move it to result directory and upload. 
   2. Write 1)the ensembl identifier of gene with the largest statistical significance and 2) the ensembl identifier of gene with the largest fold change in *volcano_genes.tsv* in result directory. 
   ```
   # Format
   ENSGXXXXXXXXXX,ENSGYYYYYYYYYYYY
   ```
- Search how to use EnhancedVolcano and fill the blanks.


## 7. heatmap_and_clustering.sh
- No blanks in the script
- Result file: 
   1. *clustered_heatmap.pdf* will be generated. Move it to result directory and upload.
   2. Based on the generated plot, whose normal sample is likely the metastasis sample? Write the ID in *probably_prep_mistake.txt* in result directory.
   3. Give an idea to improve the visualization of *clustered_heatmap.pdf* in *visual_idea.txt* in result directory.
