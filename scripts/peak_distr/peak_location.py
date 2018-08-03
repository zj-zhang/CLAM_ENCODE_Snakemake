## annotate peak locations, i.e. exon, intron or UTR
## Zijun Zhang
## 5.22.2018

import sys
#import multiprocessing
import pybedtools
import pandas as pd


gff = '/u/home/f/frankwoe/nobackup/hg19/gencode.v19.annotation.sorted.gtf'
bed = '/u/home/f/frankwoe/scratch/m6A_CLAM/src/annotate_peaks/all_peak.sorted.bed'

g = pybedtools.BedTool(gff).remove_invalid().saveas()


def featuretype_filter(feature, featuretype):
	if feature[2] == featuretype:
		return True
	return False


def subset_featuretypes(featuretype):
	result = g.filter(featuretype_filter, featuretype).saveas()
	return pybedtools.BedTool(result.fn)


def intersect_peaks_in_features(features_fn, inverse=False, write_both=True):
	"""
	Callback function to count reads in features
	"""
	return pybedtools.BedTool(bed).intersect(
							 b=features_fn,
							 stream=False, v=inverse, wo = write_both)


def record_peaks_in_features(intersection_fn):
	peak_set = set()
	with open(intersection_fn, 'r') as f:
		for line in f:
			ele = line.strip().split()
			peak_set.add(ele[3])
	return peak_set
	


featuretypes = ('exon', 'CDS', 'UTR')
exons = subset_featuretypes('exon')
cds = subset_featuretypes('CDS')
utrs = subset_featuretypes('UTR')

exons = exons.fn
cds_only = cds.subtract(utrs).sort().merge().remove_invalid().saveas().fn
utr_only = utrs.subtract(cds).sort().merge().remove_invalid().saveas().fn
cds_and_utr = utrs.intersect(cds).sort().merge().remove_invalid().saveas().fn


features = (
	['intron', exons, True], 
	['CDS_only', cds_only, False], 
	['UTR_only', utr_only, False], 
	['CDS_and_UTR', cds_and_utr, False])

results = {}
for feature in features:
	results[feature[0]] = intersect_peaks_in_features(feature[1], feature[2])


allpeak = record_peaks_in_features(bed)
peak_location_df = pd.DataFrame(0, index=allpeak, columns=['CDS_only', 'UTR_only', 'CDS_and_UTR', 'intron'])	

for feature in results:
	peakset = record_peaks_in_features(results[feature].fn)
	peak_location_df.loc[peakset, feature] += 1


peak_location_df.to_csv('peak_localization_df.csv.gz', sep='\t', compression='gzip')
pybedtools.helpers.cleanup(verbose=True, remove_all=True)

