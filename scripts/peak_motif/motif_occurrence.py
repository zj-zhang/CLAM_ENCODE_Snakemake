"""
for RBPs with definitive motifs, use the motif
occurrences as a evaluation to benchmark clam and piranha
results

Author:
	ZZJ

Date:
	2019.7.3
"""

import sys
from collections import defaultdict
import os
#from pygr import seqdb
import pysam
import gzip
import re
from tqdm import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')	


def read_file_peak(peak_fn, fn_id, peak_dict, score_column=8):
	"""
	score_column (0-based):  
		6: signalValue - Measurement of overall (usually, average) enrichment for the region.
		7: pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
		8: qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
	"""
	with open(peak_fn, 'r') as f:
		for line in f:
			ele = line.strip().split()
			peak_id = ':'.join([ele[0], ele[1], ele[2], ele[5] ])
			peak_dict[peak_id][fn_id] = float(ele[score_column])
	return peak_dict


def read_clam(peak_dict, project):
	clam_dirs = [x for x in os.listdir(os.path.join("projects", project, "clam")) if x.startswith('peaks-')]
	for clam_dir in clam_dirs:
		peak_fn = os.path.join('projects', project, 'clam', clam_dir, 'narrow_peak.unique.bed')
		fn_id = ":".join([clam_dir.lstrip('peaks-'), 'CLAM'])
		peak_dict = read_file_peak(peak_fn, fn_id, peak_dict, 8)

		#peak_fn = os.path.join('projects', project, 'clam', clam_dir, 'narrow_peak.combined.bed')
		#fn_id = ".".join([clam_dir, 'clam_m'])
		#peak_dict = read_file_peak(peak_fn, fn_id, peak_dict, 8)
	return peak_dict


def read_piranha(peak_dict, project):
	p_dirs = [x for x in os.listdir(os.path.join("projects", project, "piranha"))]
	for p_dir in p_dirs:
		peak_fn = os.path.join('projects', project, 'piranha', p_dir, 'piranha_out.bed')
		fn_id = ":".join([p_dir, 'Piranha'])
		peak_dict = read_file_peak(peak_fn, fn_id, peak_dict, 6)
	return peak_dict


def fetch_seq(chr, start, end, strand, genome):
	def reverse_complement(dna):
		complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
		return ''.join([complement[base] for base in dna[::-1]])
	try:
		seq = genome.fetch(chr, int(start), int(end))
		if strand=='-':
			seq = reverse_complement(seq)
	except:
		raise Exception('pysam cannot fetch sequeces: {}'.format(":".join([chr, str(start), str(end), strand])))
	return seq


def match_motif(seq, motif_list):
	match_list = []
	for motif in motif_list:
		if not motif.startswith("("): motif = "(" + motif
		if not motif.endswith(")"): motif = motif + ")"
		match_list.append( str( len( re.findall(motif, seq) ) ) )
	return match_list


def parse_motif(peak_dict, motif_list, genome):
	motif_dict = {}
	for peak in tqdm(peak_dict):
		chrom, start, end, strand = peak.split(':')
		peak_seq = fetch_seq(chrom, start, end, strand, genome)
		motif_count = match_motif(peak_seq, motif_list)
		motif_dict[peak] = motif_count
	return motif_dict


def extract_experi_motif_info(fn_id, peak_dict, motif_dict, motif_len):
	target_peaks = []
	for peak in peak_dict:
		if fn_id in peak_dict[peak]:
			target_peaks.append(peak)
	expr_df = pd.DataFrame(index=np.arange(len(target_peaks)), columns=['peak', 'experi', 'method', 'pval', 'motif', 'peak_len'])
	for i in tqdm(range(len(target_peaks))):
		peak = target_peaks[i]
		experi, method = fn_id.split(":")
		pval = peak_dict[peak][fn_id]
		motif = int(motif_dict[peak][0])
		chrom, start, end, strand = peak.split(":")
		peak_len = int(end) - int(start)
		this_row = {'peak': peak, 'experi':experi, 'method':method,
			'pval':pval, 'peak_len':peak_len,
			'motif':motif}
		expr_df.loc[i] = this_row
	expr_df = expr_df.sort_values(by=['pval'], ascending=True)

	expr_df['motif_prop_peak'] = pd.Series(expr_df['motif']>0).cumsum() / np.arange(1, expr_df.shape[0]+1)
	expr_df['motif_covg_peak'] = pd.Series(expr_df['motif']*motif_len).cumsum() / pd.Series(expr_df['peak_len']).cumsum()

	return expr_df


def plot_lines(expr_df, ax1, ax2, label='label', linestyle='-', color='red', skip_first_num = 100):
	ax1.plot(np.arange(skip_first_num, expr_df.shape[0]), expr_df.motif_prop_peak[skip_first_num:], label=label, ls=linestyle, color=color)
	if ax2 is not None:
		ax2.plot(np.arange(skip_first_num, expr_df.shape[0]), expr_df.motif_covg_peak[skip_first_num:], label=label, ls=linestyle, color=color)
	return ax1, ax2


def reset_style():
	from matplotlib import rcParams
	rcParams['font.family'] = 'serif'
	rcParams['font.serif'] = ['Times New Roman']
	rcParams['axes.titlesize'] = 14
	rcParams['axes.labelsize'] = 14
	rcParams['lines.linewidth'] = 1.5
	rcParams['lines.markersize'] = 8
	rcParams['xtick.labelsize'] = 14
	rcParams['ytick.labelsize'] = 14
	rcParams['legend.fontsize'] = 14


def make_figure(peak_dict, motif_dict, save_fn, motif):

	fn_id_list = [
		('rep1_2:CLAM', 'red', '-'),
		('rep1:CLAM', 'purple', '-'),
		('rep2:CLAM', 'orange', '-'),
		('rep1:Piranha', 'blue', '-'),
		('rep2:Piranha', 'green', '-'),
	]
	expr_df_list = []
	for fn_id, color, linestyle in fn_id_list:
		expr_df = extract_experi_motif_info(fn_id, peak_dict, motif_dict, motif_len=5)
		expr_df_list.append(expr_df)


	plt.close()
	reset_style()
	fig = plt.figure(figsize=(6, 6))
	#ax1, ax2 = fig.add_subplot(121), fig.add_subplot(122)
	ax1 = fig.add_subplot(111)

	for i in range(len(expr_df_list)):
		ax1, _  = plot_lines(expr_df_list[i], ax1, None, 
			label=fn_id_list[i][0], 
			linestyle=fn_id_list[i][2],
			color=fn_id_list[i][1],
			skip_first_num=1000)

	ax1.legend(loc='upper right')
	ax1.set_xlabel('Number of peaks ranked by significance')
	ax1.set_ylabel('Proportion of peaks with "%s" motif'%motif)
	ax1.set_ylim(0.2, 0.6)

	plt.savefig(save_fn)

def run(project, genome_fn, save_fn, motif):

	#genome_fn='/mnt/isilon/xing_lab/zhangz4/genome-build/hg19/hg19.noRand.fa'
	genome = pysam.FastaFile(genome_fn)

	#project='GSE77629_RBFOX2_293T_eCLIP'

	peak_dict = defaultdict(lambda:defaultdict(float))
	peak_dict = read_clam(peak_dict, project)
	peak_dict = read_piranha(peak_dict, project)

	#motif_list = ["(GCATG)"]
	motif_list = [motif]
	motif_dict = parse_motif(peak_dict, motif_list, genome)

	# plot & save
	#save_fn = 'Benchmark_CLAM-RBFOX2_eCLIP.pdf'
	make_figure(peak_dict, motif_dict, save_fn)


if __name__ == '__main__':
	project, genome_fn, save_fn, motif = sys.argv[1:]
	run(project, genome_fn, save_fn, motif)	
