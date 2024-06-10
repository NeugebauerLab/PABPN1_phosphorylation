#import required packages
import numpy as np
import pandas as pd
import csv
from tqdm import tqdm
import statistics
from scipy import stats
import scipy.stats as sci
import warnings
import math
import argparse

#set up parser
parser = argparse.ArgumentParser()

parser.add_argument('-i', help='isoform ids mapped to transcriptome from flair quantify output')
parser.add_argument('-p', help='polyA tails of reads file from get_tails.py')
parser.add_argument('-o', help='outfile name for isoforms df')

args = parser.parse_args()

#function to convert gene_ids to isoform IDs from flair output
def get_isoforms(isoform_file, polyA_file):
    file = open(isoform_file, 'r')
    polyA_df = pd.read_csv(polyA_file, sep = '\t')
    reads_dict = {}
    for line in tqdm(file):
        tnx = line.split('\t')[0]
        reads = line.split('\t')[1].split(',')
        for read in reads:
            reads_dict[read] = tnx
    mapped_df = pd.DataFrame.from_dict(reads_dict, orient = 'index').rename(columns = {0 : 'isoform_id'})
    mapped_df['r_name'] = mapped_df.index
    merged_df = polyA_df.merge(mapped_df, on = 'r_name')
    merged_df = merged_df[['isoform_id', 'tail_len', 'r_name']].rename(columns = {'isoform_id' : 'gene_id'})
    return(merged_df)

#run the function on inputs from parse
df = get_isoforms(args.i, args.p)

#write the output to a df
df.to_csv(args.o, sep = '\t', index = False)