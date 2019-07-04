## convert CLAM output into bigwig
## modified from CLAM evaluator back in 2016
## Zijun Zhang
## 11.6.2017
## revised 6.25.2018: add `genome` to arguments

import sys
import subprocess
import os

# get helper scripts path
cur_dir = os.path.dirname(os.path.realpath(__file__))
# bam_script: clam_bam2bed.py aligner_out is_stranded +/-
bam_script = os.path.join(cur_dir, 'clam_bam2bed.py')
# bam_script: clam_bam2bed.py score_column
bed_script = os.path.join(cur_dir, 'clam_bed2bedG.py')
# bedGraphToBigWig; ucsc binary
bedGtoBw = os.path.join(os.path.dirname(cur_dir), 'UCSC', 'bedGraphToBigWig')


# parse command-line args
unique_bam = sys.argv[1]
aligner_out = sys.argv[2]
out_bg_dir = sys.argv[3]
genome = sys.argv[4]
is_stranded = True if len(sys.argv)>=6 and sys.argv[5].lower().startswith('t') else False
verbose = False

chrom_sizes_fn = os.path.join(os.path.dirname(cur_dir), 'UCSC', genome+'.chrom.sizes')


# make unique bedGraph
if is_stranded:
	#cmd = 'bedtools genomecov -bg -split -strand + -ibam ' + unique_bam + ' -g /u/home/f/frankwoe/nobackup/programs/UCSC/hg19.chrom.sizes | LC_COLLATE=C sort -k1,1 -k2,2n > ' + out_bg_dir + '/unique_pos.bdg'
	cmd = 'python %s %s %s %s %s'%(bam_script, unique_bam, is_stranded, '+', 'False')+''' | awk '{print $1"\t"$2"\t"$3"\t"$5}' ''' + '| LC_COLLATE=C sort -k1,1 -k2,2n | python ' + bed_script + ' ' + chrom_sizes_fn + ' 4 | LC_COLLATE=C sort -k1,1 -k2,2n > ' + out_bg_dir + '/unique_pos.bdg'
else:
	#cmd = 'bedtools genomecov -bg -split -ibam ' + unique_bam + ' -g /u/home/f/frankwoe/nobackup/programs/UCSC/hg19.chrom.sizes | LC_COLLATE=C sort -k1,1 -k2,2n > ' + out_bg_dir + '/unique_pos.bdg'
	cmd = 'python %s %s %s %s %s'%(bam_script, unique_bam, is_stranded, '+', 'False')+''' | awk '{print $1"\t"$2"\t"$3"\t"$5}' ''' + '| LC_COLLATE=C sort -k1,1 -k2,2n | python ' + bed_script + ' ' + chrom_sizes_fn + ' 4 | LC_COLLATE=C sort -k1,1 -k2,2n > ' + out_bg_dir + '/unique_pos.bdg'
if verbose:
	print(cmd)
subprocess.call(cmd, shell=True)
#cmd = 'bedtools genomecov -bg -split -strand - -ibam ' + unique_bam + ' -g /u/home/f/frankwoe/nobackup/programs/UCSC/hg19.chrom.sizes | LC_COLLATE=C sort -k1,1 -k2,2n > ' + out_bg_dir + '/unique_neg.bdg'
cmd = 'python %s %s %s %s %s'%(bam_script, unique_bam, is_stranded, '-', 'False')+''' | awk '{print $1"\t"$2"\t"$3"\t"$5}' ''' + '| LC_COLLATE=C sort -k1,1 -k2,2n | python ' + bed_script + ' ' + chrom_sizes_fn + ' 4 | LC_COLLATE=C sort -k1,1 -k2,2n > ' + out_bg_dir + '/unique_neg.bdg'
if is_stranded:
	if verbose: print(cmd)
	subprocess.call(cmd, shell=True)

# make multiple bedGraph
if is_stranded:
	cmd = 'python %s %s %s %s'%(bam_script, aligner_out, is_stranded, '+')+''' | awk '{print $1"\t"$2"\t"$3"\t"$5}' ''' + '| LC_COLLATE=C sort -k1,1 -k2,2n | python ' + bed_script + ' ' + chrom_sizes_fn + ' 4 | LC_COLLATE=C sort -k1,1 -k2,2n > ' + out_bg_dir + '/multi_pos.bdg'
	if verbose: print(cmd)
	subprocess.call(cmd, shell=True)

	cmd = 'python %s %s %s %s'%(bam_script, aligner_out, is_stranded, '-')+''' | awk '{print $1"\t"$2"\t"$3"\t"$5}' ''' + '| LC_COLLATE=C sort -k1,1 -k2,2n | python ' + bed_script + ' ' + chrom_sizes_fn + ' 4 | LC_COLLATE=C sort -k1,1 -k2,2n > ' + out_bg_dir + '/multi_neg.bdg'
	if verbose: print(cmd)
	subprocess.call(cmd, shell=True)
else:
	cmd = 'python %s %s %s %s'%(bam_script, aligner_out, is_stranded, '+')+''' | awk '{print $1"\t"$2"\t"$3"\t"$5}' ''' + '| LC_COLLATE=C sort -k1,1 -k2,2n | python ' + bed_script + ' ' + chrom_sizes_fn + ' 4 | LC_COLLATE=C sort -k1,1 -k2,2n > ' + out_bg_dir + '/multi_pos.bdg'
	if verbose: print(cmd)
	subprocess.call(cmd, shell=True)

# take union to get combined bedGraph

cmd = 'bedtools unionbedg -i ' + out_bg_dir + '/multi_pos.bdg' + ' ' + out_bg_dir + '''/unique_pos.bdg | awk '{print $1"\t"$2"\t"$3"\t"$4+$5}' >''' + out_bg_dir + '/combined_pos.bdg'
if verbose: print(cmd)
subprocess.call(cmd, shell=True)

cmd = 'bedtools unionbedg -i ' + out_bg_dir + '/multi_neg.bdg' + ' ' + out_bg_dir + '''/unique_neg.bdg | awk '{print $1"\t"$2"\t"$3"\t"$4+$5}' >''' + out_bg_dir + '/combined_neg.bdg'
if is_stranded:
	if verbose: print(cmd)
	subprocess.call(cmd, shell=True)

# convert combined bedGraph to BigWig

cmd = bedGtoBw + ' ' + out_bg_dir + '/combined_pos.bdg ' + chrom_sizes_fn + ' ' + out_bg_dir + '/combined_pos.bw'
if verbose: print(cmd)
subprocess.call(cmd, shell=True)

cmd = bedGtoBw + ' ' + out_bg_dir + '/combined_neg.bdg ' + chrom_sizes_fn + ' ' + out_bg_dir + '/combined_neg.bw'
if is_stranded:
	if verbose: print(cmd)
	subprocess.call(cmd, shell=True)

# convert unique bedGraph to BigWig
cmd = bedGtoBw + ' ' + out_bg_dir + '/unique_pos.bdg ' + chrom_sizes_fn + ' ' + out_bg_dir + '/unique_pos.bw'
if verbose: print(cmd)
subprocess.call(cmd, shell=True)

cmd = bedGtoBw + ' ' + out_bg_dir + '/unique_neg.bdg ' + chrom_sizes_fn + ' ' + out_bg_dir + '/unique_neg.bw'
if is_stranded:
	if verbose: print(cmd)
	subprocess.call(cmd, shell=True)

# clean up
os.remove(os.path.join(out_bg_dir, "unique_pos.bdg"))
os.remove(os.path.join(out_bg_dir, "multi_pos.bdg"))
os.remove(os.path.join(out_bg_dir, "combined_pos.bdg"))
if is_stranded:
	os.remove(os.path.join(out_bg_dir, "unique_neg.bdg"))
	os.remove(os.path.join(out_bg_dir, "multi_neg.bdg"))
	os.remove(os.path.join(out_bg_dir, "combined_neg.bdg"))
