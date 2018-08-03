__author__ = "Zijun Zhang"
__version__ = "0.0.1"
__status__ = "Debug"

import sys
import pysam, itertools
from collections import defaultdict
from time import strftime


def stranded_read_start(read):
	if read.is_reverse:
		return read.positions[-1]
	else:
		return read.pos


def barcode_collapse(in_bam):
	number_of_unmapped_mate_pairs = 0
	number_of_unpaired = 0
	different_chroms = 0
	result_dict = {}
	location_to_reads = defaultdict(list)
	read_to_locations = defaultdict(list)
	processed = 0
	with pysam.Samfile(in_bam, 'r') as samfile1:
		with pysam.Samfile(in_bam, 'r') as samfile2:
			samfile_read1 = itertools.islice(samfile1, 0, None)
			samfile_read2 = itertools.islice(samfile2, 0, None)
			samfile_read2.next()
			#for read1, read2 in itertools.izip(samfile_read1, samfile_read2):
			while True:
				try:
					processed += 1
					#if processed == 1000000:
					#	break
					if not processed % 1000000:
						print_time_stamp("%i unique reads processed." % processed)
					read1 = samfile_read1.next()
					read2 = samfile_read2.next()
					if not read1.qname == read2.qname:
						number_of_unpaired += 1
						continue
					if read1.is_unmapped and read1.is_unmapped:
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
					#read1.is_read1
					strand = "-" if read1.is_reverse else "+"
					if read1.opt('NH') == 1:
						unique_location = (read1.rname, start, stop, strand, randomer)
						if unique_location in result_dict:
							result_dict[(read1.rname, start, stop, strand, randomer)][2] += 1
							continue
						result_dict[(read1.rname, start, stop, strand, randomer)] = [read1, read2, 1]
					else:
						multiple_location = (read1.rname, start, stop, strand, randomer)
						location_to_reads[multiple_location].append([read1, read2, randomer])
						read_to_locations[(read1.qname, randomer)].append(multiple_location)
					
					samfile_read1.next()
					samfile_read2.next()
				except StopIteration:
					break
	print_time_stamp("Processed unique reads, found %i unpaired reads." % number_of_unpaired)
	return result_dict, location_to_reads, read_to_locations


def search_read_subg(qname, randomer, read_to_locations, location_to_reads, seen):
	queue = [(qname, randomer)]
	tmp_subg = set()
	subg = []
	subg_seq = {}
	while len(queue) > 0:
		new_queue = []
		for read in queue:
			this_locs = read_to_locations[read]
			for loc in this_locs:
				for read1, read2, this_randomer in location_to_reads[loc]:
					if this_randomer==randomer and (not (read1.qname, this_randomer) in tmp_subg):
						this_read = (read1.qname, this_randomer)
						new_queue.append(this_read)
						subg_seq[this_read] = [read1, read2]
		map(seen.add, queue)
		map(tmp_subg.add, queue)
		queue = list(set(new_queue))	
	subg = list(tmp_subg)
	return subg, subg_seq, seen


def collapse_subg(subg, subg_seq):
	read_pair = [None, None]
	max_len = 0
	for read in subg:
		read1, read2 = subg_seq[read]
		this_len = read1.query_length + read2.query_length
		if  this_len > max_len:
			read_pair = [read1, read2]
			max_len = this_len
	return read_pair + [len(subg)]

	
def collapse_multi_reads(location_to_reads, read_to_locations):
	seen = set()
	multi_reads = []
	processed = 0
	for qname, randomer in read_to_locations:
		processed += 1
		if not processed % 100000:
			print_time_stamp("%i multiple reads process" % processed)
		if (qname, randomer) in seen:
			continue
		subg, subg_seq, seen = search_read_subg(qname, randomer, read_to_locations, location_to_reads, seen)
		read_pair = collapse_subg(subg, subg_seq)
		multi_reads.append(read_pair)
	return multi_reads

def rev_complement(seq):
	complement_letter = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
	complement_seq = [complement_letter[x] for x in seq]
	return ''.join(complement_seq[::-1])
	
def output_fa(unique_reads, multi_reads, read_to_locations, location_to_reads, fa_file):
	output_f = open("%s_1.fasta" % (fa_file), "w")
	output_r = open("%s_2.fasta" % (fa_file), "w")
	for unique_location in unique_reads:
		read1, read2, number_of_duplicates = unique_reads[unique_location]
		seq_f = rev_complement(read1.query_sequence) if read1.is_reverse else read1.query_sequence
		seq_r = rev_complement(read2.query_sequence) if read2.is_reverse else read2.query_sequence
		output_f.write(">%s_%s/1\n%s\n" % (read1.qname, number_of_duplicates, seq_f))
		output_r.write(">%s_%s/2\n%s\n" % (read2.qname, number_of_duplicates, seq_r))
	for read in multi_reads:
		read1, read2, number_of_duplicates = read
		seq_f = rev_complement(read1.query_sequence) if read1.is_reverse else read1.query_sequence
		seq_r = rev_complement(read2.query_sequence) if read2.is_reverse else read2.query_sequence
		output_f.write(">%s_%s/1\n%s\n" % (read1.qname, number_of_duplicates, seq_f))
		output_r.write(">%s_%s/2\n%s\n" % (read2.qname, number_of_duplicates, seq_r))
	output_f.close()
	output_r.close()

def print_time_stamp(msg):
	current_time='[' + strftime("%Y-%m-%d %H:%M:%S") + '] '
	print >> sys.stderr, current_time + msg

if __name__ == "__main__":
	in_bam, out_fa = sys.argv[1], sys.argv[2]
	unique_reads, location_to_reads, read_to_locations = barcode_collapse(in_bam)
	multi_reads=collapse_multi_reads(location_to_reads, read_to_locations)

	output_fa(unique_reads, multi_reads, read_to_locations, location_to_reads, out_fa)
	