# make snakemake softlinks for parallel submission of pipeline jobs
# Zijun Zhang
# 1.3.2018
project=$1
ln -sf `readlink -f Snakefile` $project/Snakefile
ln -sf `readlink -f scripts` $project/scripts
ln -sf `readlink -f run_snakemake.sh` $project/run_snakemake.sh
ln -sf `readlink -f submit.sh` $project/submit.sh
ln -sf `readlink -f clusterconfig.yaml` $project/clusterconfig.yaml
ln -sf `readlink -f re_touch.sh` $project/re_touch.sh
