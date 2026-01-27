import pandas as pd
import pyBigWig
import pybedtools
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import re
from datetime import datetime

import matplotlib as mpl
from matplotlib import font_manager

arial_path = "/media/scratch/fy2306/tools/fonts"
font_files = font_manager.findSystemFonts(fontpaths=arial_path)

for file in font_files:
    font_manager.fontManager.addfont(file)
    
mpl.rcParams['font.family'] = 'Arial'

cons_path = "/media/rna/fy2306/hg38.phastCons100way.bw"
mane_gtf_path = "/media/scratch/fy2306/genomes/hg38/MANE/release_1.4/MANE.GRCh38.v1.4.refseq_genomic.gtf"

def calculate_conserved_nucleotide(bw, mane_cds, threshold):
	chrom_sizes = bw.chroms()

	total_bases = 0
	total_conserved = 0

	for chrom, length in bw.chroms().items():
		scores = bw.values(chrom, 0, length)
		scores = np.nan_to_num(np.array(scores), nan=0.0)
		total_bases += len(scores)
		total_conserved += np.sum(scores >= threshold)
	# bigWig and pybigwig use 0-based half-open coordinates

	cds_conserved = 0

	for _, row in mane_cds.iterrows():
		chrom, start, end = row['chrom'], int(row['start']), int(row['end'])	# now 0-based start
		if chrom not in bw.chroms():
			continue
		vals = bw.values(chrom, start, end)
		vals = np.nan_to_num(np.array(vals), nan=0.0)                  
		cds_conserved += np.sum(vals >= threshold)

	coding_pct = 100 * cds_conserved / total_conserved
	noncoding_pct = 100 - coding_pct
	noncoding_conserved = total_conserved - cds_conserved
	
	return coding_pct, noncoding_pct, cds_conserved, noncoding_conserved

def plot_conserved_nucleotide(thresholds):

	# load in bw
	bw = pyBigWig.open(cons_path)

	# load in gtf
	mane = pd.read_csv(mane_gtf_path, sep='\t', comment='#', header=None, names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
	mane['transcript_id'] = mane['attribute'].str.extract(r'transcript_id "([^"]+)"')
	mane['tag'] = mane['attribute'].str.extract(r'tag "([^"]+)"')
	mane = mane[(mane['tag'] == 'MANE Select') & (mane['feature'].isin(['CDS', 'transcript']))]
	print(mane['transcript_id'].nunique())
	mane_cds = mane[mane['feature'] == 'CDS']	# MANE contains lncRNAs. feature = CDS automatically removes lncRNA.
	print(mane_cds['transcript_id'].nunique())
	print(len(mane_cds))
	# merge intervals
	mane_cds['start'] = mane_cds['start'] - 1
	print(len(mane_cds))
	mane_cds = pybedtools.BedTool.from_dataframe(mane_cds[['seqname', 'start', 'end']])
	mane_cds = mane_cds.sort()
	mane_cds = mane_cds.merge()
	mane_cds = mane_cds.to_dataframe(names=['chrom', 'start', 'end'])
	print(len(mane_cds))

	results = []

	for threshold in thresholds:
		print(threshold)
		print(f'start time: {datetime.now()}')
		coding_pct, noncoding_pct, coding_conserved, noncoding_conserved = calculate_conserved_nucleotide(bw, mane_cds, threshold)
		results.append({'threshold': threshold, 'coding': coding_pct, 'noncoding': noncoding_pct, 'coding_number': coding_conserved, 'noncoding_number': noncoding_conserved})

	df = pd.DataFrame(results, columns=['threshold', 'coding', 'noncoding', 'coding_number', 'noncoding_number'])
	df.to_csv(
		"/media/scratch/fy2306/projects/base_editing/script/coding_vs_noncoding/whole_genome_conservation.tsv",
		sep="\t",
		index=False
	)

	# stack plot
	plt.figure(figsize=(4,2.5))
	plt.stackplot(
		df['threshold'],
		df['coding'],
		df['noncoding'],
		labels=['Coding', 'Noncoding'],
		colors=['#FF0066', '#000000']
	)

	plt.legend().set_visible(False)
	plt.xlabel('Threshold for conserved base pair', fontsize=13)
	plt.ylabel(r'Conserved base pair %', fontsize=13)
	plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], fontsize=12)
	plt.xlim(df['threshold'].min(), df['threshold'].max())
	plt.yticks(fontsize=12)
	plt.ylim(0,101)
	plt.gca().spines['top'].set_visible(False)
	plt.gca().spines['right'].set_visible(False)
	# add text
	plt.text(0.8, df['coding'].iloc[8] / 2, 'Coding', color='white', fontsize=13,
			ha='center', va='center', fontweight='bold')
	plt.text(0.2, df['noncoding'].iloc[1] / 2 + df['coding'].iloc[1], 'Noncoding', color='white', fontsize=13,
			ha='center', va='center', fontweight='bold')
	plt.tight_layout()
	# plt.savefig("/media/scratch/fy2306/projects/base_editing/script/coding_vs_noncoding/whole_genome_conservation.stack.png", dpi=300, bbox_inches='tight')
	plt.savefig("/media/scratch/fy2306/projects/base_editing/plots/miscellaneous/conserved_basepair.whole_genome.pdf", 
				bbox_inches="tight",
				dpi=300,              
				transparent=True,
				format='pdf')
plot_conserved_nucleotide([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])