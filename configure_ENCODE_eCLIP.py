#!/usr/bin/evn python
"""
Read in the metadata file and write the configuration files
for Snakemake pipeline
Zijun Zhang
3.7.2018
"""

import os
import pandas as pd


def read_metadata(fn):
	data = pd.read_table(fn)
	data.index = data['File accession']
	return data

def table_to_content(t, prefix, project_dir):
	"""
	Args
		prefix (str): 'IP' or 'Inp'
	"""
	line_formatter = "{url}\t{path}\t{md5}\n"
	s = ''
	for j in range(t.shape[0]):
		rep = t['Biological replicate(s)'][j]
		pair = t['Paired end'][j]
		url = t['File download URL'][j]
		md5 = t['md5sum'][j]
		path = os.path.join(project_dir, 'reads', "%s_rep%s"%(prefix,rep), 'rep%s_%s.fastq.gz'%(rep, pair))
		s += line_formatter.format(url=url, path=path, md5=md5) 
	return s

def write_url_file(par_dir='/u/nobackup/yxing/NOBACKUP/frankwoe/CLAM_ENCODE-eCLIP_Snakemake/projects'):
	data = read_metadata('metadata_20180307.tsv')
	written_experiments = set()
	for i in range(data.shape[0]):
		row = data.iloc[i].to_dict()
		if 'mock input' in row['Experiment target']:
			continue
		if row['Experiment accession'] in written_experiments:
			continue
		exp = data.loc[ data['Experiment accession']==row['Experiment accession'] ]
		target = exp['Experiment target'][0].split('-')[0]
		accession = exp['Experiment accession'][0]
		cell = exp['Biosample term name'][0]
		#project_dir = os.path.join(par_dir, '%s_%s_%s'%(target, cell, accession))
		project_dir = os.path.join(par_dir, '%s_%s_%s'%(target, cell, accession), 'projects','%s_%s_%s'%(target, cell, accession))
		if not os.path.isdir(project_dir):
			os.makedirs(project_dir)
		if not os.path.isdir(os.path.join(project_dir,'config')):
			os.mkdir(os.path.join(project_dir,'config'))
		if not os.path.isdir(os.path.join(project_dir,'reads')):
			os.mkdir(os.path.join(project_dir,'reads'))
		url_fh  = open(os.path.join(project_dir,'config', 'url.txt'), 'w')
		#config_fh  = open(os.path.join(project_dir,'config', 'config.yaml'), 'w')
		s1 = ''
		s2 = ''
		s1 = table_to_content(exp, 'IP', project_dir)
		url_fh.write( s1 )
		try:
			control_id = list(set([x.split('/')[2] for x in exp['Controlled by'].values]))
		except:
			#print exp['Controlled by'].values
			continue
		control = data.loc[control_id]
		s2 = table_to_content(control, 'Inp', project_dir)
		url_fh.write( s2 )
		written_experiments.add(accession)
		url_fh.close()
		#print(project_dir)
		#break
		print('%s_%s_%s'%(target, cell, accession))
	return 0

if __name__ == '__main__':
	write_url_file()