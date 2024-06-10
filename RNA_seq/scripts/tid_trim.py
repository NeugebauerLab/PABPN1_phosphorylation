#!/usr/bin/env python
#This is a script to trim transcription IDs so that they match the standard ENST.# format

import os
import pysam
import argparse
import pandas as pd
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("-q", help="input quant.sf file for transcript ID trimming")
parser.add_argument("-o", help="outfile name")
args = parser.parse_args()

df = pd.read_csv(args.q, sep = "\t")

df["Name"] = df["Name"].str.split("|").str[0]

df.to_csv(args.o, sep = "\t")