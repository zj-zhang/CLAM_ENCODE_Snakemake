# make bigwig based on the sample info
# Zijun Zhang
# 5.20.2017
# revised 5.23.2017: changed to print ucsc track link
# revised 4.6.2018: changed to print CASS link

import os
import sys

par_dir = 'https://xinglab.cass.idre.ucla.edu/public/'
child_dir = 'zijun/custom_tracks/20180616_Gastro_RCNP_FGC1866/'

#track_name_dict = {}
#with open('track_name_dict.txt', 'r') as f:
#	for line in f:
#		ele = line.strip().split('\t')
#		track_name_dict[ele[0]] = ele[1]

tracks = os.listdir('/u/home/f/frankwoe/nobackup/RussCarstens_Lab/pipeline/RNAseq_Snakemake_pipeline/Gastro_RCNP_FGC1866/bigwig')
for track in tracks:
	url = '/'.join([par_dir, child_dir, track, track+'.bw'])
	#track_name = track_name_dict[track]
	track_name = track
	s1 = "track name='{0}' bigDataUrl={1} type=bigWig visibility=2".format(track_name, url)
	print s1


