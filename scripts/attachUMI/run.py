
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
	#print(fr_fn, se_fn)
	if os.path.isfile(os.path.join(out_dir, fr_fn)):
		continue
	qname = fr_fn.split('.')[0]
	cmd = "qsub -N {} -l h_data=5G,h_rt=12:00:0 submit.sh {} {} {} {}".format(qname, os.path.join(in_dir, fr_fn), os.path.join(in_dir, se_fn), out_dir, barcode_len)
	print(cmd)
	os.system(cmd)