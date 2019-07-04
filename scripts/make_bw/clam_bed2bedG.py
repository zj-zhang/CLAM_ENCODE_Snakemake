import sys

# requires sorted bed as piped input; does not support strand information

def init(chr, start, end, score):
	return chr, [start,end], [score] * (end - start)

def pile_up(score_a, score_b):
	new_len = max(len(score_a), len(score_b))
	score = [0] * new_len
	score = [score[i] + score_a[i] if i<len(score_a) else score[i] for i in range(new_len)]
	score = [score[i] + score_b[i] if i<len(score_b) else score[i] for i in range(new_len)]
	return score

def write_out(chr, interval, score):
	if len(score)==0:
		return
	out_score = score[0]
	left = interval[0]
	right = interval[0]
	for i in range(len(score)):
		if score[i] == out_score:
			right += 1
		else:
			print('\t'.join([chr, str(left), str(right), str(out_score)]))
			left = right
			right = left + 1
			out_score = score[i]
		if right >= interval[1]:
				break
	if left!=right:
		print('\t'.join([chr, str(left), str(right), str(out_score)]))


def read_chrom_sizes(fn):
	chrom_sizes = {}
	with open(fn, 'r') as f:
		for line in f:
			chrom, size = line.strip().split()
			chrom_sizes[chrom] = int(size)
	return chrom_sizes

if __name__ == '__main__':
	chrom_sizes_fn = sys.argv[1]
	score_column = int(sys.argv[2]) or None

	init_flag = True
	chrom_sizes =read_chrom_sizes(chrom_sizes_fn)

	for line in sys.stdin:
		line = line.rstrip()
		ele = line.split('\t')
		chr, start, end = ele[0], int(ele[1]), int(ele[2])
		score = 1 if score_column is None else float(ele[score_column-1])
		if init_flag:
			last_chr, last_interval, last_score = init(chr, start, end, score)
			init_flag = False
			continue
		if chr == last_chr and start < last_interval[1]:
			next_interval = [start, max(end, last_interval[1]) ]
			next_score = pile_up(last_score[ (start-last_interval[0]):: ], [score] * (end - start))
			
			this_interval = [last_interval[0], start]
			this_score = last_score[ 0 : (start-last_interval[0]) ]
			
			if this_interval[0] > chrom_sizes[chr]-2:
				print('1', last_chr, this_interval, file=sys.stderr)
				this_interval[0] = chrom_sizes[chr]-2
			if this_interval[1] > chrom_sizes[chr]-1:
				print('1', last_chr, this_interval, file=sys.stderr)
				this_interval[1] = chrom_sizes[chr]-1

			write_out(chr, this_interval, this_score)
			
			last_interval = next_interval
			last_score = next_score
		else:
			if last_interval[0] > chrom_sizes[last_chr]-2:
				print('2', last_chr, last_interval, file=sys.stderr)
				last_interval[0] = chrom_sizes[last_chr]-2
			if last_interval[1] > chrom_sizes[last_chr]-1:
				print('2', last_chr, last_interval, file=sys.stderr)
				last_interval[1] = chrom_sizes[last_chr]-1
			write_out(last_chr, last_interval, last_score)
			last_chr, last_interval, last_score = init(chr, start, end, score)
			continue

	if last_interval[0] > chrom_sizes[last_chr]-2:
		print('3', last_chr, last_interval, file=sys.stderr)
		last_interval[0] = chrom_sizes[last_chr]-2
	if last_interval[1] > chrom_sizes[last_chr]-1:
		print('3', last_chr, last_interval, file=sys.stderr)
		last_interval[1] = chrom_sizes[last_chr]-1
	write_out(last_chr, last_interval, last_score)
			