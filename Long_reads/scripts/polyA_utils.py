import pysam
from tqdm import tqdm
import numpy as np
from collections import Counter
import re
import pandas as pd
from multiprocessing import Pool
import time


def nonA(files, min_A_frac=0.85, tail_range=200, min_intron_len=100):
    if isinstance(files, str):
        files = [files]
        
    tail_len = []
    base_count_prox = [np.zeros(tail_range), np.zeros(tail_range), np.zeros(tail_range), np.zeros(tail_range)] # A, G, U, C
    base_count_dist = [np.zeros(tail_range), np.zeros(tail_range), np.zeros(tail_range), np.zeros(tail_range)]
    base_cov = np.zeros(tail_range)
    
    total_c = Counter()
    
    for file in files:
        bamfile = pysam.AlignmentFile(file, 'rb')
        for read in tqdm(bamfile.fetch(), total=bamfile.count()):

            _, three, _ = Convert_Read(read, with_clipped=True)

            c = Counter(three)
            if c and c['A'] > 0 and c['A']/sum(c.values()) >= min_A_frac: # if it's really a poly(A) tail
                if tail_range-len(three) > 0: # if the tail is shorter than the chosen window (most cases)
                    for i, nt in enumerate(['A', 'G', 'T', 'C']):
                        base_count_prox[i] += np.pad((np.array(list(three)) == nt).astype(int), 
                                                     (0, tail_range-len(three)), 'constant', constant_values=0)
                        base_count_dist[i] += np.pad((np.array(list(three)) == nt).astype(int), 
                                                     (tail_range-len(three), 0), 'constant', constant_values=0)
                else: # if the tail is longer than the window
                    for i, nt in enumerate(['A', 'G', 'T', 'C']):
                        base_count_prox[i] += (np.array(list(three)) == nt).astype(int)[:tail_range]
                        base_count_dist[i] += (np.array(list(three)) == nt).astype(int)[-tail_range:]

                base_cov[:len(three)] += 1
                tail_len.append(len(three))
                total_c.update(three)
        bamfile.close()    
    return(total_c, tail_len, base_count_prox, base_count_dist, base_cov)


def polyA_info_bed(file, min_A_frac=0.85, bed=None, min_intron_len=100):    
            
    tail_len = []
    genes = []
    r_name = []
    introns = []
    ends_5 = []
        
    bamfile = pysam.AlignmentFile(file, 'rb')
    total_reads = bamfile.count()

    if not bed:
        for read in tqdm(bamfile.fetch(), total=total_reads):
            if not read.is_secondary and not read.is_supplementary and read.query_name in features:
                three, _, n_introns, mapped_5 = Convert_Read(read, min_intron_len=min_intron_len)
                c = Counter(three)
                if c and c['A'] > 0 and c['A']/sum(c.values()) >= min_A_frac: # if it's really a poly(A) tail
                    tail_len.append(len(three))
                    introns.append(n_introns)
                    ends_5.append(mapped_5)
                    r_name.append(read.query_name)
        bamfile.close()
        return(np.array(tail_len), np.array(introns), np.array(ends_5), np.array(r_name))

    else:
        features = pd.read_csv(bed, delimiter='\t', header=None, usecols=[3, 15]).set_index(3).to_dict()[15]
        for read in tqdm(bamfile.fetch(), total=total_reads):
            if not read.is_secondary and not read.is_supplementary and read.query_name in features:
                three, _, n_introns, mapped_5 = Convert_Read(read, min_intron_len=min_intron_len)
                c = Counter(three)
                if c and c['A'] > 0 and c['A']/sum(c.values()) >= min_A_frac: # if it's really a poly(A) tail
                    tail_len.append(len(three))
                    genes.append(features[read.query_name])
                    introns.append(n_introns)
                    ends_5.append(mapped_5)
                    r_name.append(read.query_name)
        bamfile.close()
        return(np.array(tail_len), np.array(genes), np.array(introns), np.array(ends_5), np.array(r_name))

def get_max_spliced(hsap, int_list):
    n_should_be_spliced = []
    # vectorize for SPEED increase
    int_list_chrom = np.array(int_list['chrom'])
    int_list_start = np.array(int_list['start'])
    int_list_stop = np.array(int_list['stop'])
    int_list_strand = np.array(int_list['strand'])
    int_list_gene_id = np.array(int_list['gene_id'])

    hsap_n_spliced_introns = np.array(hsap['n_spliced_introns'])
    hsap_gene_id = np.array(hsap['gene_id'])
    hsap_chrom = np.array([i.split(' ')[0].strip('[').strip("'") for i in list(hsap['read_end_pos'])])
    hsap_pos = np.array([int(i.split(' ')[1].strip(']').strip("'")) for i in list(hsap['read_end_pos'])])

    # iterate over reads
    for i in range(len(hsap_gene_id)):
        # for the assigned geneID, get the intron positions and strand
        imask = int_list_gene_id==hsap_gene_id[i]
        strand = np.unique(int_list_strand[imask])
        if strand.size != 1:
            n_should_be_spliced.append(np.nan)
            continue
        # if on + strand: count introns that have start > read 5'-end
        elif strand[0] == '+':
            n_introns = np.sum(int_list_start[imask] > hsap_pos[i])
        # if on - strand: count introns that have end < read 5'-end
        elif strand[0] == '-':
            n_introns = np.sum(int_list_stop[imask] < hsap_pos[i])

        if n_introns < hsap_n_spliced_introns[i]:
            pass
        if hsap_chrom[i] != str(int_list_chrom[imask][0]):
            n_should_be_spliced.append(np.nan)
            continue
        # add number of potentially spliced introns to df
        n_should_be_spliced.append(n_introns)
        
    return(n_should_be_spliced)

def tail_metrics_per_gene(lengths, genes, min_reads=25):
    
    per_gene_res = dict()
    for i in tqdm(np.unique(genes)):
        temp = lengths[genes==i]
        
        if (len(temp) >= min_reads):
            per_gene_res[i] = [np.mean(temp), np.median(temp), np.quantile(temp, 0.25), np.quantile(temp, 0.75), pearson_skew_2(temp), len(temp)]
            
    per_gene = pd.DataFrame.from_dict(per_gene_res, orient='index', columns=['mean', 'median', "quantile_25", "quantile_75", 'skew', 'n'])
    return(per_gene)

                       
def reverse_complement(dna):
    complement = {'N': 'N', 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    return(''.join([complement[base] for base in dna[::-1]]))

def complement(dna):
    complement = {'N': 'N', 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    return(''.join([complement[base] for base in dna]))

def Convert_Read(read, min_intron_len=100):
    """
    Parses the cigar string to find where the clipped bases are. Returns 5'- and 3'- soft clipped bases.
    Input: pysam.fetch() object
    """

    # 3'-end clipped bases are marked in the bit vector, 5'-end bases are not

    if read.is_unmapped:
        return('')

    clipped_3 = ''
    clipped_5 = ''

    # if the read is on the - strand we have to get the complementary base
    if read.is_reverse:
        reverse = True
        read_seq = reverse_complement(read.get_forward_sequence()).upper()
        mapped_5 = [read.reference_name, read.reference_end-1]
    else:
        reverse = False
        read_seq = read.get_forward_sequence().upper()  # Sequence of the read
        mapped_5 = [read.reference_name, read.reference_start+1]

    j = 0  # Pos in the read sequence
    l = 0  # Pos in the ref position list
    intron_counter = 0
    CIGAR_Ops = Parse_CIGAR(read.cigarstring)
    op_index = 0
    while op_index < len(CIGAR_Ops):  # Each CIGAR operation
        op = CIGAR_Ops[op_index]
        desc, length = op[1], int(op[0])

        if desc in ['M', 'I']:  # Match or mismatch
            j += length
            l += length

        elif desc == 'N': # potential intron
            if length >= min_intron_len:
                intron_counter += 1

            
        elif desc in ['D', 'H']: 
            pass

        elif desc == 'S':  # Soft clipping
            if (not reverse and op_index == len(CIGAR_Ops) - 1) or (reverse and op_index == 0):  # Soft clipped at the 3'-end
                clipped_3 = read_seq[j:j+length] # clipped 3'-end bases for studying poly(A) tail properties
            else: # Soft clipped at 5'-end
                clipped_5 = read_seq[j:j+length] # clipped 5'-end bases for studying RT-stops
                
            j += length  # Update read index
            l += length

        else:
            print('Unknown CIGAR op encountered: %s'%(desc))
            break

        op_index += 1

    # return the vector and clipped bases if desired
    if reverse:
        return(reverse_complement(clipped_3), reverse_complement(clipped_5), intron_counter, mapped_5)
    else:
        return(clipped_3, clipped_5, intron_counter, mapped_5)



def Parse_CIGAR(cigar_string):
    """
    Parse a CIGAR string
    Args:
        cigar_string (string): CIGAR string
    Returns:
        ops (list): List of operations. Each op is of type
        (length, description) such as ('37', 'M'), ('10', 'I'), (24, 'D'), etc.
    """
    ops = re.findall(r'(\d+)([A-Z]{1})', cigar_string)
    return(ops)


def pearson_skew_1(X):
    return((np.mean(X)-mode(X)[0][0])/np.std(X))

def pearson_skew_2(X):
    return(3*(np.mean(X)-np.median(X))/np.std(X))