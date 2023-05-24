#!/usr/bin/env python3

import os
import pickle as pkl
import pandas as pd
import argparse

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# 1. Parse arguments
parser = argparse.ArgumentParser(description='DE analysis using pyDESeq2')
# parser.add_argument('-i', '--input', help='Input count matrix', required=True)

# 2. Load data
# 2.1. Load count matrix and sample info into DataFrame
sample_info = pd.DataFrame({'SRR453567': {'condition': 'batch'}, 'SRR453570' : {'condition':'chem'}}).T
counts = pd.read_csv('result/gene_count_matrix.csv', index_col=0).T

genes_to_keep = counts.columns[counts.sum(axis=0) > 10]
counts_filtered = counts[genes_to_keep]

dds = DeseqDataSet(counts=counts_filtered, clinical=sample_info, design_factors='condition', refit_cooks=True, n_cpus=8)
dds.fit_size_factors()
print(dds.obsm['size_factors'])

# print(sample_info)
# print(counts_filtered)