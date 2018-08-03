"""Collapse PCR duplicates based on barcode and align-locus
for ENCODE eCLIP data

Zijun Zhang
3.5.2018
"""

from collections import Counter, defaultdict
import itertools
from optparse import OptionParser
import sys
import os
import pysam
import operator

def stranded_read_start(read):
	if read.is_reverse:
		return read.positions[-1]
	else:
		return read.pos

		
def output_metrics(metrics_file, total_count, removed_count):
	with open(metrics_file, 'w') as metrics:
		metrics.write("\t".join(["randomer", "total_count", "removed_count"]) + "\n")
		for barcode in total_count.keys():
			metrics.write("\t".join(map(str, [barcode, total_count[barcode], removed_count[barcode]])) + "\n")

def collapse_reads(result_dict):
	collapsed_read_dict = {}
	read_to_len = defaultdict(int)
	for loc in result_dict:
		# use read2, as this will be used for peak-calling
		qname_to_rlen = {qname:read2.rlen for qname,(read1,read2) in result_dict[loc].items()}
		max_rlen_qname = max(qname_to_rlen.iteritems(), key=operator.itemgetter(1))[0]
		max_rlen = max(qname_to_rlen.values())
		for read_qname in qname_to_rlen:
			if max_rlen_qname > read_to_len[read_qname]:
				collapsed_read_dict[read_qname] = max_rlen_qname
				read_to_len[read_qname] = max_rlen
	return collapsed_read_dict


def barcode_collapse(in_bam, out_bam):
	number_of_unmapped_mate_pairs = 0
	different_chroms = 0
	number_of_unpaired = 0
	removed_count = Counter()
	total_count = Counter()
	
	# (chr,start,end,strand,randomer) => { qname1: [original_read_qname, rlen],
	# qname2: [original_read_qname2, rlen2], ...}
	result_dict = defaultdict(dict) 
	
	number_of_processed = 0
	samfile1= pysam.Samfile(in_bam, 'rb') 
	samfile2 = pysam.Samfile(in_bam, 'rb') 
	samfile_read1 = itertools.islice(samfile1, 0, None, 1)
	samfile_read2 = itertools.islice(samfile2, 0, None, 1)
	samfile_read2.next()  # shift read2 by 1 alignment
	#for read1, read2 in itertools.izip(samfile_read1, samfile_read2):
	while True:
		try:
			read1 = samfile_read1.next()
			read2 = samfile_read2.next()
		except StopIteration:
			break
		number_of_processed += 1
		if not number_of_processed%10**6:
			print(number_of_processed)
		if not read1.qname == read2.qname:
			number_of_unpaired += 1
			continue
		if read1.is_unmapped and read2.is_unmapped:
			#Both reads don't map, don't even both saving them.
			continue
		if (not read1.is_unmapped and read2.is_unmapped) or (read1.is_unmapped and read2.is_unmapped):
			number_of_unmapped_mate_pairs += 1
			continue
		if read1.rname != read2.rname:
			different_chroms += 1
			continue
		#if the read order is swapped swap everything before running.
		if not read1.is_read1:
			read1, read2 = read2, read1
		randomer = read1.qname.split(":")[0]
		start = stranded_read_start(read1)
		stop = stranded_read_start(read2)
		strand = "-" if read1.is_reverse else "+"
		unique_location = (read1.rname, start, stop, strand, randomer)
		# shift one location before operation
		try:
			samfile_read1.next()
			samfile_read2.next()
		except StopIteration:
			pass
		
		total_count[randomer] += 1
		if unique_location in result_dict and read2.opt('NH')==1:
			removed_count[randomer] += 1
			continue
		
		# chr, start, end, strand, randomer => AlignedSegment1, AlignedSegment2 
		#result_dict[(read1.rname, start, stop, strand, randomer)][read2.qname] = (read1, read2)
		result_dict[unique_location][read2.qname] = (None, read2)
		
	collapsed_read_dict = collapse_reads(result_dict)
	with pysam.Samfile(out_bam, 'wb', template=samfile1) as out:
		for loc in result_dict:
			target_qname_set = set()
			for qname, (read1, read2) in result_dict[loc].items():
				#print qname, read2
				target_qname = collapsed_read_dict[qname]
				target_qname_set.add(target_qname)
			for target_qname in list(target_qname_set):
				if target_qname in result_dict[loc]:
					out.write(result_dict[loc][target_qname][1]) ## only write read2
		samfile1.close()
		samfile2.close()
	return total_count, removed_count

if __name__ == "__main__":
	description=""""Paired End randomer aware duplciate removal algorithm."""
	usage="""Assumes paired end reads are adjecent to each other in output file 
			(ie *MUST* only provide *unsorted* bams)
			 """
	parser = OptionParser(usage=usage, description=description)
	parser.add_option("-b", "--bam", dest="bam", help="bam file to barcode collapse")
	parser.add_option("-o", "--out_file", dest="out_file")
	parser.add_option("-m", "--metrics_file", dest="metrics_file")
	(options, args) = parser.parse_args()

	if not (options.bam.endswith(".bam")):
		raise TypeError("%s, not bam file" % options.bam)

	total_count, removed_count = barcode_collapse(options.bam, options.out_file)
	output_metrics(options.metrics_file, total_count, removed_count)

	sys.exit(0)
