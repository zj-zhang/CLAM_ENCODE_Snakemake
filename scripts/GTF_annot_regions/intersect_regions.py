'''
Assign peaks to genomic regions
Zijun Zhang
8.1.2018
'''

import sys
import os
import pybedtools
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict

### input arguments
peak_fp, genome, outfn = sys.argv[1], sys.argv[2], sys.argv[3]
par_dir = os.path.dirname(os.path.realpath(__file__))


### make pybedtools objects
peaks = pybedtools.BedTool(peak_fp)
ref_dict = {
	'exon': pybedtools.BedTool(os.path.join(par_dir, genome, 'exons.bed')),
	'3UTR': pybedtools.BedTool(os.path.join(par_dir, genome, '3UTRs.bed')),
	'5UTR': pybedtools.BedTool(os.path.join(par_dir, genome, '5UTRs.bed')),
	'cds': pybedtools.BedTool(os.path.join(par_dir, genome, 'cds.bed')),
	'intron': pybedtools.BedTool(os.path.join(par_dir, genome, 'introns.bed')),
	'proximal200': pybedtools.BedTool(os.path.join(par_dir, genome, 'proximal200_intron.bed')),
	'proximal500': pybedtools.BedTool(os.path.join(par_dir, genome, 'proximal500_intron.bed'))
}

get_total_length = lambda featureset: sum([x.end-x.start-1 for x in featureset])

### process reference for use
target = {
	"3UTR": ref_dict['3UTR'],
	"5UTR": ref_dict['5UTR'],
	"CDS": ref_dict['cds'],
	"other_exon": ref_dict['exon']-ref_dict['3UTR']-ref_dict['5UTR']-ref_dict['cds'],
	"px200_intron": ref_dict['proximal200'],
	"px500_intron": ref_dict['proximal500'].subtract(ref_dict['proximal200']),
	"distal_intron": ref_dict['intron'].subtract(ref_dict['exon']).subtract(ref_dict['proximal500'])
}

category_list = ['3UTR', '5UTR', 'CDS', 'other_exon', "px200_intron", "px500_intron", "distal_intron"]

### compute counts and length-normalized counts
len_dict = {x:get_total_length(target[x]) for x in category_list}

#count_dict = {x: (peaks+target[x]).count() for x in category_list}
intersect_handler_dict = {x: (peaks+target[x]) for x in category_list}
count_dict = {x: intersect_handler_dict[x].count() for x in category_list}


count_df = pd.DataFrame({
	'category':category_list, 
	'count':[count_dict[x] for x in category_list],
	'norm_count':[count_dict[x]/float(len_dict[x]) for x in category_list],
	})
count_df.norm_count /= sum(count_df.norm_count)

count_df.to_csv(outfn, sep='\t')
### plot
fig = plt.figure()
ax  =fig.add_subplot(111)
#ax.tick_params(axis='x',rotation=45)
ax2 = ax.twinx()
width = 0.25

count_df.plot(x='category', y='count', kind='bar', color='red', ax=ax, width=width, position=1)
count_df.plot(x='category', y='norm_count', kind='bar', color='blue', ax=ax2, width=width, position=0)

ax.set_ylabel('count')
ax2.set_ylabel('length-normalized count')
fig.tight_layout()
plt.xticks(rotation=45)
fig.savefig(outfn.rstrip('txt')+'png')
plt.close()

### output the annotations
annot_fn = outfn+'.annot'
annot_dict = defaultdict(list)
for x in intersect_handler_dict:
	for peak in intersect_handler_dict[x]:
		# chrom, start, end, gene, strand
		annot_dict[(peak[0], peak[1], peak[2], peak[3], peak[5])].append(x)

with open(annot_fn, 'w') as fo:
	for peak in annot_dict:
		fo.write("\t".join(peak)+'\t'+','.join(annot_dict[peak])+'\n')
