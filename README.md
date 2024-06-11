# PABPN1 Phosphorylation: RNA-seq analyses

# Preprocessing, mapping, splicing, and poly(A) tail lengths quantification (Long-read Sequencing)

Preprocessing, mapping, splicing, and poly(A) tail lengths quantification can be carried out with the attached Snakemake pipeline. Be sure to create/adjust the sample.yaml file as needed.

Make sure the following dependencies are installed:

[Prinseq](https://prinseq.sourceforge.net/)\
[Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)\
[minimap2](https://lh3.github.io/minimap2/minimap2.html)\
[Samtools](https://www.htslib.org/)\
[Bedtools](https://bedtools.readthedocs.io/en/latest/index.html)\
[Flair](https://flair.readthedocs.io/en/latest/)\
[QAPA](https://github.com/morrislab/qapa)\
[bam2bakR](https://github.com/simonlabcode/bam2bakR)\
[bakR](https://github.com/simonlabcode/bakR)


# Plotting and correlative analyses

Code for making plots and carrying out statistical and correlative comparisons (e.g. tail length vs. RNA stability) can be found the attached jupyter notebooks. 

Code is arranged according to Figure number.

 Figure 3 - Quantifying poly(A) tail lengths and non-A nucleotides in tails.\
 Figure 4 - Comparing poly(A) tail lengths by gene and splice isoform.\
 Figure 5 - RNA turnover from TimeLapse-seq.\
 Figure S5 - PacBio library characteristics, non-A nucleoties in tails, and APA.\
 Figure S6 - Comparing poly(A) tail lengths by gene and differential gene expression.\
 Figure S7 - RNA synthesis rates from TimeLapse-seq.

