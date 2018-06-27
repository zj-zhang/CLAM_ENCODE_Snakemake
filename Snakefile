#!/usr/bin/env python

"""A pipeline for processing ENCODE eCLIP data using CLAM
Author: 
	Zijun Zhang <zj.z@ucla.edu>
Date: 
	3.6.2018
Histroy:
	6.24.2018: update to general eCLIP
"""

import os
import re
import itertools

###**--- UTILITY FUNCTION ---**###

def concat_star_fq(sample_name):
	fq1 = []
	fq2 = []
	par_dir = os.path.join('projects', PROJECT, 'reads')
	suffix = '.fastq.gz' if GZIP else '.fastq' 
	for sample_id in config['sample_dict'][sample_name]:
		fq1.append( os.path.join(par_dir, sample_name, sample_id+'_1'+suffix+'.adapterTrim.round2'+suffix) )
		fq2.append( os.path.join(par_dir, sample_name, sample_id+'_2'+suffix+'.adapterTrim.round2'+suffix) )
	if PAIRED_END:
		s = ','.join(fq1) + ' ' + ','.join(fq2)
	else:
		s = ','.join(fq1)
	#print(s)
	return s

###**--- IN-FILE CONFIG ---**###
	
PROJECT = os.environ.get("PROJECT")
CONFIG_FP = os.path.join("projects", PROJECT, 'config', 'config.yaml')

configfile: 
	#"config.template.yaml"
	CONFIG_FP

GZIP = True if 'gzip' not in config['reads'] else config['reads']['gzip']
PAIRED_END = False if 'paired_end' not in config['reads'] else config['reads']['paired_end']
FQ_SUFFIX = '.fastq.gz' if GZIP else '.fastq'
FQ_PATTERN = ['projects/{project}/reads/{sample_name}/{sample_id}_1'+FQ_SUFFIX,
			'projects/{project}/reads/{sample_name}/{sample_id}_2'+FQ_SUFFIX]

FQ_TRIM_PATTERN = [	'projects/{project}/reads/{sample_name}/{sample_id}_1'+FQ_SUFFIX+'.adapterTrim.round2'+FQ_SUFFIX,
		'projects/{project}/reads/{sample_name}/{sample_id}_2'+FQ_SUFFIX+'.adapterTrim.round2'+FQ_SUFFIX]
			
GENOME = config['genome']
FQ_CMD = '--readFilesCommand zcat' if GZIP else ''

SAMPLE_DICT = config['sample_dict']
COMPARISON_LIST = config['clam']['sample_comparison'].keys()
MAX_TAGS = config['clam']['max_tags']

INCLUDE_MREAD_ANALYSIS = True

###**--- SNAKE RULES ---**###

rule all:
	input:
		# pre-processing
		# ["projects/{project}/star/{sample_name}/Aligned.out.mask_rRNA.dup_removed.r2.bam".format(
			# project=PROJECT, sample_name=x) for x in SAMPLE_DICT
		# ],
		# clam unique peaks
		# ["projects/{project}/clam/peaks-{comparison}/narrow_peak.unique.bed".format(
			# project=PROJECT, comparison=x) for x in config['clam']['sample_comparison']
		# ],
		# clam combined peaks
		# ["projects/{project}/clam/peaks-{comparison}/narrow_peak.combined.bed".format(
			# project=PROJECT, comparison=x) for x in config['clam']['sample_comparison']
		# ],
		# compare peaks
		#["projects/{project}/clam/peaks-{comparison}/peak_num.png".format(
		#	project=PROJECT, comparison=x) for x in config['clam']['sample_comparison']
		#],
		# homer
		#expand(["projects/{project}/homer/{comparison}/clam_unique/homerResults.html","projects/{project}/homer/{comparison}/clam_rescue/homerResults.html"],
		#	project=PROJECT, comparison=COMPARISON_LIST),
		# archive
		"projects/{project}/archive/{project}.tar.gz".format(project=PROJECT)

rule download_fq:
	input:
		url_fn="projects/{project}/config/url.txt".format(project=PROJECT)
	output:
		["projects/{project}/reads/{sample_name}/{sample_id}_1{suffix}".format(project=PROJECT, sample_id=y, sample_name=x, suffix=FQ_SUFFIX)
			for x in SAMPLE_DICT.keys() for y in SAMPLE_DICT[x] ],
		["projects/{project}/reads/{sample_name}/{sample_id}_2{suffix}".format(project=PROJECT, sample_id=y, sample_name=x, suffix=FQ_SUFFIX)
			for x in SAMPLE_DICT.keys() for y in SAMPLE_DICT[x] ]

	shell:
		"cat {input} | python2 scripts/downloader/eclip_downloader.py"


rule cutadapt:
	input:
		sample=lambda wildcards: 
			expand(FQ_PATTERN, 
				project=PROJECT, 
				sample_name=wildcards.sample_name, 
				sample_id=SAMPLE_DICT[wildcards.sample_name])
	output:
		FQ_TRIM_PATTERN
	shell:
		"bash scripts/cutadapt/cutadapt_caller.sh {input.sample}"


rule star_map:
	input:
		lambda wildcards:
			expand(FQ_TRIM_PATTERN, project=PROJECT, sample_name=wildcards.sample_name,sample_id=SAMPLE_DICT[wildcards.sample_name])
	output:
		"projects/{project}/star/{sample_name}/Aligned.out.bam"
	params:
		fq_input = lambda wildcards: concat_star_fq(wildcards.sample_name),
		prefix = "projects/{project}/star/{sample_name}/",
		index=config['genome_build'][GENOME]['star_idx'],
		fq_cmd = FQ_CMD,
		max_hits = 100
	threads: 4
	shell:
		"""
STAR --genomeDir {params.index} \
--readFilesIn {params.fq_input}  --outSAMtype BAM Unsorted \
--outFileNamePrefix {params.prefix} \
--outFilterMultimapNmax {params.max_hits} \
--runThreadN 4 \
{params.fq_cmd} 
		"""


rule mask_rRNA:
	input:
		"projects/{project}/star/{sample_name}/Aligned.out.bam"
	output:
		"projects/{project}/star/{sample_name}/Aligned.out.mask_rRNA.bam"
	params:
		rna_annot = config['reads']['rna_annot']
	shell:
		"bash scripts/mask_rRNA/mask_rRNA.sh {input} {output} {params.rna_annot}"


rule collapse_dup:
	input:
		"projects/{project}/star/{sample_name}/Aligned.out.mask_rRNA.bam"
	output:
		"projects/{project}/star/{sample_name}/Aligned.out.mask_rRNA.dup_removed.r2.bam"
	params:
		star_bam="projects/{project}/star/{sample_name}/Aligned.out.bam",
		metric_file="projects/{project}/star/{sample_name}/dup_removal.metrics.txt",
		collapse_script="scripts/collapse_pcr/collapse_duplicates.py"
	shell:
		"""
python2 scripts/collapse_pcr/collapse_duplicates.py -b {input} -o {output} -m {params.metric_file}
rm {input} {params.star_bam}
		"""


### peak calling

rule clam_prep:
	input:
		align="projects/{project}/star/{sample_name}/Aligned.out.mask_rRNA.dup_removed.r2.bam"
	output:
		"projects/{project}/clam/{sample_name}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else \
			"projects/{project}/clam/{sample_name}/unique.sorted.bam"
	log:
		"projects/{project}/logs/clam/{sample_name}-prep.log"
	params:
		outdir="projects/{project}/clam/{sample_name}",
		tagger_method="start",
		max_tags=MAX_TAGS
	shell:
		"CLAM preprocessor -i {input.align} -o {params.outdir} --max-tags {params.max_tags} --read-tagger-method {params.tagger_method} >{log} 2>&1"


rule clam_em:
	input:
		"projects/{project}/clam/{sample_name}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else \
			"projects/{project}/clam/{sample_name}/unique.sorted.bam"
	output:
		"projects/{project}/clam/{sample_name}/realigned.sorted.bam"
	log:
		"projects/{project}/logs/clam/{sample_name}-em.log"
	params:
		outdir="projects/{project}/clam/{sample_name}",
		max_tags=MAX_TAGS,
		winsize=50
	shell:
		"CLAM realigner -i {input} -o {params.outdir} --winsize {params.winsize} --max-tags {params.max_tags} --read-tagger-method start >{log} 2>&1"
		

rule clam_callpeak:
	input:
		ip_bam=lambda wildcards: 
			expand("projects/{project}/clam/{ip_sample_name}/unique.collapsed.sorted.bam" if MAX_TAGS>0 else "projects/{project}/clam/{ip_sample_name}/unique.sorted.bam", 
			project=PROJECT,
			ip_sample_name=config['clam']['sample_comparison'][wildcards.comparison][0] ),
		control_bam=lambda wildcards: 
			expand("projects/{project}/clam/{con_sample_name}/unique.collapsed.sorted.bam" if MAX_TAGS>0 else "projects/{project}/clam/{con_sample_name}/unique.sorted.bam", 
			project=PROJECT,
			con_sample_name=config['clam']['sample_comparison'][wildcards.comparison][1] )
	output:
		"projects/{project}/clam/peaks-{comparison}/narrow_peak.unique.bed"
	log:
		"projects/{project}/logs/clam/{comparison}-callpeak.log"
	params:
		outdir = "projects/{project}/clam/peaks-{comparison}",
		gtf = config['genome_build'][GENOME]['gtf'],
		binsize = 50,
		qval_cutoff = 0.05,
		fold_change = '0.69',  # log(2)
		threads = 4,
		ip_bam = lambda wildcards: ','.join(expand("projects/{project}/clam/{ip_sample_name}/unique.collapsed.sorted.bam" if MAX_TAGS>0 else "projects/{project}/clam/{ip_sample_name}/unique.sorted.bam", 
			project=PROJECT,
			ip_sample_name=config['clam']['sample_comparison'][wildcards.comparison][0] )),
		control_bam = lambda wildcards: ','.join(expand("projects/{project}/clam/{con_sample_name}/unique.collapsed.sorted.bam" if MAX_TAGS>0 else "projects/{project}/clam/{con_sample_name}/unique.sorted.bam", 
			project=PROJECT,
			con_sample_name=config['clam']['sample_comparison'][wildcards.comparison][1] )),
		pool = lambda wildcards: config['clam']['sample_comparison'][wildcards.comparison][2] if len(config['clam']['sample_comparison'][wildcards.comparison])>2 else '',
	shell:
		"""
CLAM peakcaller -i {params.ip_bam}  -c {params.control_bam} \
-p {params.threads} \
-o {params.outdir} --gtf {params.gtf} --unique-only --binsize {params.binsize} \
--qval-cutoff 0.5 --fold-change 0.01 {params.pool} >{log} 2>&1
mv {output} {output}.all
awk '$9<{params.qval_cutoff} && $7>{params.fold_change}' {output}.all > {output}
		"""


rule clam_callpeak_mread:
	input:		
		ip_bam= lambda wildcards: expand([
			"projects/{project}/clam/{ip_sample_name}/unique.collapsed.sorted.bam" if MAX_TAGS>0 else "projects/{project}/clam/{ip_sample_name}/unique.sorted.bam", 
			"projects/{project}/clam/{ip_sample_name}/realigned.sorted.bam"], 
			project=PROJECT, ip_sample_name=config['clam']['sample_comparison'][wildcards.comparison][0]), 
		con_bam= lambda wildcards: expand([
			"projects/{project}/clam/{con_sample_name}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else "projects/{project}/clam/{con_sample_name}/unique.sorted.bam", 
			"projects/{project}/clam/{con_sample_name}/realigned.sorted.bam"],
			project=PROJECT, con_sample_name=config['clam']['sample_comparison'][wildcards.comparison][1])
	output:
		"projects/{project}/clam/peaks-{comparison}/narrow_peak.combined.bed"
	log:
		"projects/{project}/logs/clam/{comparison}-callpeak_mread.log"
	params:
		ip_ubam = lambda wildcards: ','.join(expand(
			"projects/{project}/clam/{ip_sample_name}/unique.collapsed.sorted.bam" if MAX_TAGS>0 else "projects/{project}/clam/{ip_sample_name}/unique.sorted.bam", 
			project=PROJECT, ip_sample_name=config['clam']['sample_comparison'][wildcards.comparison][0])),
		ip_mbam = lambda wildcards: ','.join(expand(
			"projects/{project}/clam/{ip_sample_name}/realigned.sorted.bam",
			project=PROJECT, ip_sample_name=config['clam']['sample_comparison'][wildcards.comparison][0])),
		con_ubam = lambda wildcards: ','.join(expand(
			"projects/{project}/clam/{con_sample_name}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else "projects/{project}/clam/{con_sample_name}/unique.sorted.bam", 
			project=PROJECT, con_sample_name=config['clam']['sample_comparison'][wildcards.comparison][1])),
		con_mbam = lambda wildcards: ','.join(expand(
			"projects/{project}/clam/{con_sample_name}/realigned.sorted.bam",
			project=PROJECT, con_sample_name=config['clam']['sample_comparison'][wildcards.comparison][1])),
		outdir="projects/{project}/clam/peaks-{comparison}",
		gtf=config['genome_build'][GENOME]['gtf'],
		binsize=50,
		qval_cutoff=0.05,
		fold_change='0.69',  # log(2)
		threads=4
	shell:
		"""
CLAM peakcaller -i {params.ip_ubam} {params.ip_mbam} -c {params.con_ubam} {params.con_mbam} \
-p {params.threads} \
-o {params.outdir} --gtf {params.gtf} --binsize {params.binsize} \
--qval-cutoff 0.5 --fold-change 0.01 >{log} 2>&1
mv {output} {output}.all
awk '$9<{params.qval_cutoff} && $7>{params.fold_change}' {output}.all > {output}
		"""

rule clipper:
	input:
		lambda wildcards: "projects/{project}/star/{sample_name}/Aligned.out.mask_rRNA.dup_removed.r2.bam".format(project=PROJECT, sample_name=config['clipper']['sample_comparison'])
	output:
		"projects/{project}/clipper/{comparison}/clipper_out.bed"
	params:
		threads=4,
		genome=GENOME
	shell:
		"""
samtools index {input}
clipper -b {input} -p {params.threads} -s {params.genome} -o {output}
		"""
		
		
rule compare_peaks:
	input:
		#clipper_peak="projects/{project}/clipper/{comparison}/clipper_out.bed",
		clam_upeak="projects/{project}/clam/peaks-{comparison}/narrow_peak.unique.bed",
		clam_mpeak="projects/{project}/clam/peaks-{comparison}/narrow_peak.combined.bed"
	output:
		clam_rescued_peak="projects/{project}/clam/peaks-{comparison}/narrow_peak.rescue.bed",
		plot_fn="projects/{project}/clam/peaks-{comparison}/peak_num.png"
	params:
		plot_script="scripts/report/plot_peak_num.R",
	shell:
		"""
bedtools intersect -v -a {input.clam_mpeak} -b {input.clam_upeak} > {output.clam_rescued_peak}
Rscript {params.plot_script} {input.clam_mpeak} {input.clam_upeak} {output.plot_fn}
		"""


### evaluations

rule homer_motif:
	input:
		peak_fn="projects/{project}/clam/peaks-{comparison}/narrow_peak.unique.bed"
	output:
		"projects/{project}/homer/{comparison}/clam_unique/homerResults.html"
	params:
		outdir="projects/{project}/homer/{comparison}/clam_unique",
		motif_len='5,6,7',
		genome=GENOME,
		nthread=4,
		size=100,
		motif_num=10
	log:
		"projects/{project}/logs/homer/log.homer.{comparison}.txt"
	shell:
		"findMotifsGenome.pl {input.peak_fn} {params.genome} {params.outdir} "\
		" -rna -len {params.motif_len} "\
		"-p {params.nthread} -size {params.size} -S {params.motif_num} >{log} 2>&1"


rule repeat_comp:
	input:
		peak_fn="projects/{project}/clam/peaks-{comparison}/narrow_peak.unique.bed"
	output:
		"projects/{project}/repeats/{comparison}/clam_unique/dist.png"
	params:
		genome=GENOME,
		outdir="projects/{project}/repeats/{comparison}/clam_unique",
		count_script="scripts/repeat_composition/RepeatsPie.py",
		plot_script="scripts/repeat_composition/plot.R"
	shell:
		"""
		python2 {params.count_script} {input.peak_fn} {params.genome} >{params.outdir}/dist.data
		Rscript {params.plot_script} {params.outdir}/dist.data {output}
		"""


rule homer_motif_rescue:
	input:
		peak_fn="projects/{project}/clam/peaks-{comparison}/narrow_peak.rescue.bed"
	output:
		"projects/{project}/homer/{comparison}/clam_rescue/homerResults.html"
	params:
		outdir="projects/{project}/homer/{comparison}/clam_rescue",
		motif_len='5,6,7',
		genome=GENOME,
		nthread=4,
		size=100,
		motif_num=10
	log:
		"projects/{project}/logs/homer/log.homer.{comparison}.rescue.txt"
	shell:
		"findMotifsGenome.pl {input.peak_fn} {params.genome} {params.outdir} "\
		" -rna -len {params.motif_len} "\
		"-p {params.nthread} -size {params.size} -S {params.motif_num} >{log} 2>&1"


rule repeat_comp_rescue:
	input:
		peak_fn="projects/{project}/clam/peaks-{comparison}/narrow_peak.rescue.bed"
	output:
		"projects/{project}/repeats/{comparison}/clam_rescue/dist.png"
	params:
		genome=GENOME,
		outdir="projects/{project}/repeats/{comparison}/clam_rescue",
		count_script="scripts/repeat_composition/RepeatsPie.py",
		plot_script="scripts/repeat_composition/plot.R"
	shell:
		"""
python2 {params.count_script} {input.peak_fn} {params.genome} >{params.outdir}/dist.data
Rscript {params.plot_script} {params.outdir}/dist.data {output}
		"""


rule make_bw:
	input:
		"projects/{project}/clam/{sample_name}/realigned.sorted.bam"
	output:
		"projects/{project}/bigwig/{sample_name}/foo.txt"
	params:
		genome=GENOME,
		is_stranded=False,
		mbam="projects/{project}/clam/{sample_name}/realigned.sorted.bam",
		ubam="projects/{project}/clam/{sample_name}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else \
			"projects/{project}/clam/{sample_name}/unique.sorted.bam",
		bw_dir="projects/{project}/bigwig/{sample_name}/",
		bw_script="scripts/make_bw/make_bigwig.py"
	shell:
		"""
mkdir -p {params.bw_dir}
python2 {params.bw_script} {params.ubam} {params.mbam} {params.bw_dir} {params.genome} {params.is_stranded}
echo "`date` done making bw" > {output}
		"""

rule make_peak_bb:
	input:
		unique_peak="projects/{project}/clam/peaks-{comparison}/narrow_peak.unique.bed",
		combined_peak="projects/{project}/clam/peaks-{comparison}/narrow_peak.combined.bed"
	output:
		unique_peak_bb="projects/{project}/bigwig/peaks-{comparison}/unique_peak.bb",
		combined_peak_bb="projects/{project}/bigwig/peaks-{comparison}/combined_peak.bb",
	params:
		outdir="projects/{project}/bigwig/peaks-{comparison}",
		sorted_unique="projects/{project}/bigwig/peaks-{comparison}/unique.sorted.bed",
		sorted_combined="projects/{project}/bigwig/peaks-{comparison}/combined.sorted.bed",
		bedToBigBed="scripts/UCSC/bedToBigBed",
		chrom_size="scripts/UCSC/%s.chrom.sizes"%GENOME,
	shell:
		"""
mkdir -p {params.outdir}
LC_ALL=c sort -k1,1 -k2,2n {input.unique_peak} > {params.sorted_unique}
{params.bedToBigBed} -type=bed6+4 {params.sorted_unique} {params.chrom_size} {output.unique_peak_bb}
rm {params.sorted_unique}
LC_ALL=c sort -k1,1 -k2,2n {input.combined_peak} > {params.sorted_combined}
{params.bedToBigBed} -type=bed6+4 {params.sorted_combined} {params.chrom_size} {output.combined_peak_bb}
rm {params.sorted_combined}
		"""

rule merged_bw:
	input:
		group_bw_tracks = lambda wildcards: ["{project}/bigwig/{sample_name}/foo.txt".format(project=PROJECT, sample_name=x) for x in config['ucsc_tracks']['merge'][wildcards.group]]
	output:
		"{project}/bigwig_merge/{group}/{group}.bw"
	params:
		script="scripts/merge_bw/merge_bw.py",
		genome_index=GENOME,
		group_bw_tracks = lambda wildcards: ["{project}/bigwig/{sample_name}/combined_pos.bw".format(project=PROJECT, sample_name=x) for x in config['ucsc_tracks']['merge'][wildcards.group]]
	shell:
		"""
python {params.script} {output} {params.genome_index} {input}
		"""



### generating reports and cleaning up

rule mapping_stats:
	input:
		align="projects/{project}/star/{sample_name}/Aligned.out.mask_rRNA.dup_removed.r2.bam"
	output:
		"projects/{project}/star/{sample_name}/mapping_stats.txt"
	shell:
		"python2 scripts/report/mapping_stat.py {input.align} > {output}"


rule report:
	input:
		# require mapping stats
		[ "projects/{project}/star/{sample_name}/mapping_stats.txt".format(
			project=PROJECT, 
			sample_name=x)
			for x in config['sample_dict']
		],
		# require peak comparison
		[ "projects/{project}/clam/peaks-{comparison}/peak_num.png".format(
			project=PROJECT, 
			comparison=x
			)
			for x in COMPARISON_LIST
		],
		# require unique homer
		[ "projects/{project}/homer/{comparison}/clam_unique/homerResults.html".format(
			project=PROJECT, 
			comparison=x )
			for x in COMPARISON_LIST
		],
		# require unique repeats
		[ "projects/{project}/repeats/{comparison}/clam_unique/dist.png".format(
			project=PROJECT, 
			comparison=x )
			for x in COMPARISON_LIST
		],
		# require rescued homer
		[ "projects/{project}/homer/{comparison}/clam_rescue/homerResults.html".format(
			project=PROJECT, 
			comparison=x )
			for x in COMPARISON_LIST
		],
		# require rescued repeats
		[ "projects/{project}/repeats/{comparison}/clam_rescue/dist.png".format(
			project=PROJECT, 
			comparison=x )
			for x in COMPARISON_LIST
		],
	output:
		"projects/{project}/reports/report_{project}.pdf".format(project=PROJECT)
	params:
		out_html="projects/{project}/reports/report_{project}.html".format(project=PROJECT)
	run:
		from scripts.report import generate_report
		import pdfkit
		pardir = os.getcwd()
		generate_report.generate_report(config['clam']['sample_comparison'], pardir, params.out_html, PROJECT, include_mread_analysis=INCLUDE_MREAD_ANALYSIS)
		pdfkit.from_file(params.out_html, output[0])

rule archive:
	input:
		# clam peak-calling
		clam_peak = [
			"projects/{project}/clam/peaks-{comparison}/narrow_peak.unique.bed".format( 
				project=PROJECT, comparison=x)
				for x in COMPARISON_LIST
				],
		# bigwig coverage and peak bigbed
		bigwig = [ "projects/{project}/bigwig/{sample_name}/foo.txt".format(
				project=PROJECT, sample_name=x )
				for x in config['sample_dict']
				],
		peak_bb = ["projects/{project}/bigwig/peaks-{comparison}/unique_peak.bb".format(
				project=PROJECT, comparison=x)
				for x in COMPARISON_LIST
				],
		# merged bigwig
		merged_bws = ["{project}/bigwig_merge/{group}/{group}.bw".format(project=PROJECT, group=x) \
			for x in config['ucsc_tracks']['merge']],
		# report
		report = "projects/{project}/reports/report_{project}.pdf".format(project=PROJECT),
	output:
		"projects/{project}/archive/{project}.tar.gz".format(project=PROJECT)
	params:
		project=PROJECT
	shell:
		"""
tar -czvf {output} projects/{params.project}/clam/peaks-* {input.report} projects/{params.project}/bigwig/*
## rm -rf projects/{params.project}/reads/* projects/{params.project}/reads/*
echo "`date` done archiving" > projects/{params.project}/foo.txt
		"""
