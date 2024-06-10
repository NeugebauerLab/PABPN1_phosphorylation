#Merge mapped reads from remapping to the transcriptome with all isoforms
#Outputs ann file that can be used in get_tails.py

#Import required packages
import os
import sys
import pandas as pd
import argparse
import tqdm
import numpy as np

#Set up argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input mapped reads in bed format")
parser.add_argument("-a", help="isoforms .bed file from flair collapse")
parser.add_argument("-o", help="name for outfile")
args = parser.parse_args()

#Define arguments
reads = args.i
isoforms = args.a
out_name = args.o

#Define column names for downstream and load reads .bed into a dataframe
col_names = ['isoform', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
iso_col_names = ['chr', 'start', 'end', 'isoform', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']

df = pd.read_csv(reads, sep = '\t', header = None, names = col_names)

#Read in all isoforms as df
all_isoforms = pd.read_csv(isoforms, sep = '\t', header = None, names = iso_col_names)

#Merge the isoforms and reads dataframes
merged_df = df.merge(all_isoforms, left_on = 'isoform', right_on = 'isoform')
isos = merged_df['isoform']
merged_df.insert(loc=15, column='isoform2', value=isos)

#Write out the merged df into .bed file format
merged_df.to_csv(out_name, sep = '\t', header = None, index = False)