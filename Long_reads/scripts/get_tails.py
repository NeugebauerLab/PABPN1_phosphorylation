import argparse
from polyA_utils import polyA_info_bed
import sys

parser = argparse.ArgumentParser()

parser.add_argument('bam', help='The bam file containing poly(A)-containing reads.')
parser.add_argument('-bed', help='Bed file containing an annotation for each read.')
parser.add_argument('-i', help="Set this to get number of spliced introns and read 5'-end position in the output", action='store_true')
parser.add_argument('--minAFraction', help='Fraction of As the tail must include to be considered poly(A). Default is 0.9', type=float, metavar='float')
parser.add_argument('--minIntronLength', help='Minimum length a gap has to have to be considered an intron. Default is 100', type=int, metavar='int')
args = parser.parse_args()

# print(args)

if not args.minAFraction:
	args.minAFraction = 0.9

if not args.i:

	if args.bed:
		length, genes, _, _, r_name = polyA_info_bed(args.bam, bed=args.bed, min_A_frac=args.minAFraction)

		sys.stdout.write('gene_id\ttail_len\tr_name\n')
		for i in range(len(length)):
			sys.stdout.write(str(genes[i]) + '\t' + str(length[i]) + '\t' + str(r_name[i]) + '\n')
	else:
		length, _, _, r_name = polyA_info_bed(args.bam, min_A_frac=args.minAFraction)

		sys.stdout.write('tail_len\tr_name\n')
		for i in length:
			sys.stdout.write(str(i) + '\t' + str(r_name[i]) + '\n')

else:
	if not args.minIntronLength:
		args.minIntronLength = 100

	if args.bed:
		length, genes, introns, ends, r_name = polyA_info_bed(args.bam, bed=args.bed, min_A_frac=args.minAFraction, min_intron_len=args.minIntronLength)

		sys.stdout.write('gene_id\ttail_len\tn_spliced_introns\tread_end_pos\tr_name\n')
		for i in range(len(length)):
			sys.stdout.write(str(genes[i]) + '\t' + str(length[i]) + '\t' + str(introns[i]) + '\t' + str(ends[i]) + '\t' + str(r_name[i]) + '\n')
	else:
		length, introns, ends, r_name = polyA_info_bed(args.bam, min_A_frac=args.minAFraction, min_intron_len=args.minIntronLength)

		sys.stdout.write('tail_len\tn_spliced_introns\tread_end_pos\tr_name\n')
		for i in range(len(length)):
			sys.stdout.write(str(length[i]) + '\t' + str(introns[i]) + '\t' + str(ends[i]) + '\t' + str(r_name[i]) + '\n')