#import packages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
from tqdm import tqdm
import statistics
from scipy import stats
from polyA_utils import tail_metrics_per_gene
import warnings
import argparse
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input polyA tail dataframe")
parser.add_argument("-m", help="minimum tail length to be included")
parser.add_argument("-o", help="out table name")
args = parser.parse_args()
file = args.i
min_len = args.m
outfile = args.o

#Load data
df = pd.read_csv(file, delimiter='\t')
df = df.loc[df['tail_len'] >= int(min_len)]
len_df = np.array(df['tail_len'])
genes_df = np.array(df['gene_id'])

#calculate tails per gene
per_gene_df = tail_metrics_per_gene(len_df, genes_df)

#write out
per_gene_df.to_csv(args.o, sep = "\t", header = True, index = True)
