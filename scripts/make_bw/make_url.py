# make bigwig based on the sample info
# Zijun Zhang
# 5.20.2017
# revised 5.23.2017: changed to print ucsc track link
# revised 4.6.2018: changed to print CASS link

import os
import sys

par_dir = 'https://xinglab.cass.idre.ucla.edu/public/'
child_dir = 'zijun/custom_tracks/20180625_GI_ESRP-CLIP_FGC1862/'

#track_name_dict = {}
#with open('track_name_dict.txt', 'r') as f:
#	for line in f:
#		ele = line.strip().split('\t')
#		track_name_dict[ele[0]] = ele[1]

tracks = os.listdir('/u/nobackup/yxing/NOBACKUP/frankwoe/RussCarstens_Lab/eclip/CLAM_ENCODE_Snakemake/projects/ESRP1_mm10_GI/bigwig')
for track in tracks:
	if 'peaks-' in track:
		url = '/'.join([par_dir, child_dir, track, 'unique_peak.bb'])
		s1 = "track name='{0}' bigDataUrl={1} type=bigBed visibility=2".format(track, url)
	else:
		url = '/'.join([par_dir, child_dir, track, 'unique_pos.bw'])
		#track_name = track_name_dict[track]
		track_name = track
		s1 = "track name='{0}' bigDataUrl={1} type=bigWig visibility=2".format(track_name, url)
	print(s1)


