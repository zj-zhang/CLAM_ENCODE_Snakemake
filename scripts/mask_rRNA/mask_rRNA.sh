#!/bin/bash
## mask alignments with >90% in rRNA,tRNA,snRNA in a bam file
## Zijun Zhang
## 3.6.2018

bam=$1
out=$2
annot=$3
bedtools intersect -f 0.90 -abam $bam -b $annot -v > $out