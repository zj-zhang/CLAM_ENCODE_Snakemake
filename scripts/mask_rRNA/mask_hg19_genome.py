## mask all rRNA, tRNA and snRNA from hg19
## Zijun Zhang
## 3.6.2018

from collections import defaultdict


def read_rna_annot(fn):
	annot = defaultdict(list)
	with open(fn, 'r') as f:
		for line in f:
			ele = line.strip().split()
			annot[ele[]]