#!/bin/bash
#$ -l h_data=2G,h_rt=24:00:0
#$ -S /bin/bash
#$ -R y
#$ -V
#$ -cwd
#$ -j y
#$ -m bea
#$ -M zzj.zju@gmail.com
#


source activate py36
echo "PROJECT=$1"
./run_snakemake.sh $1 
