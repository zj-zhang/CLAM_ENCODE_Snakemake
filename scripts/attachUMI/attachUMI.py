'''pre-processing for eCLIP paired-end reads
Take the first N bp off from second mate and attach that 
to the read name
Zijun Zhang
June 21, 2018
'''

import os
import sys
import gzip

in_dir = sys.argv[1]
out_dir = sys.argv[2]
lane = sys.argv[3]
barcode_len = int(sys.argv[4]) if len(sys.argv)>4 else 10

first_file_list = sorted([x for x in os.listdir(in_dir) if '%s_1'%lane in x])
second_file_list = sorted([x for x in os.listdir(in_dir) if '%s_2'%lane in x])

assert len(first_file_list) == len(second_file_list)

for fr_fn, se_fn in zip(first_file_list, second_file_list):
	print(fr_fn, se_fn)
	i = 0
	with gzip.GzipFile(os.path.join(out_dir,fr_fn), 'wb') as fr_o, \
			gzip.GzipFile(os.path.join(out_dir,se_fn), 'wb') as se_o, \
			gzip.open(os.path.join(in_dir,fr_fn), 'rb') as fr_i, \
			gzip.open(os.path.join(in_dir,se_fn), 'rb') as se_i:
		while True:
			i += 1
			if not i % 10**5:
				print(i)
			fr_header = str(fr_i.readline())
			fr_seq = str(fr_i.readline())
			fr_foo = str(fr_i.readline())
			fr_qual = str(fr_i.readline())
			se_header = str(se_i.readline())
			se_seq = str(se_i.readline())
			se_foo = str(se_i.readline())
			se_qual = str(se_i.readline())
			if len(fr_header)<5:
				break
			# modify original lines
			barcode = se_seq[0:barcode_len]
			se_seq = se_seq[barcode_len::]
			se_qual = se_qual[barcode_len::]
			fr_header = '@' + barcode + ':' + fr_header.lstrip('@')
			se_header = '@' + barcode + ':' + se_header.lstrip('@')			
			# write out
			fr_o.write(''.join([fr_header, fr_seq, fr_foo, fr_qual]))
			se_o.write(''.join([se_header, se_seq, se_foo, se_qual]))
			#break
