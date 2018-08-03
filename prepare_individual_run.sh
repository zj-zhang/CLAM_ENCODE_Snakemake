#!/bin/bash
FOLDER=$1
if [ ! -d "projects/$FOLDER" ]; then
  echo "invalid project $FOLDER"
  exit 0
fi
cd projects/$FOLDER
echo `pwd`
if [ -d "projects/$FOLDER" ]; then
  #mkdir -p projects/$FOLDER
  mv reads projects/$FOLDER/
  #rm -r reads
  mv config projects/$FOLDER/
  ln -sf `readlink -f ../../Snakefile` Snakefile
  ln -sf `readlink -f ../../scripts` scripts
  ln -sf  `readlink -f ../../clusterconfig.yaml` clusterconfig.yaml
  ln -sf  `readlink -f ../../config.template.yaml` config.template.yaml
  ln -sf  `readlink -f ../../submit.sh` submit.sh
  ln -sf `readlink -f ../../run_snakemake.sh` run_snakemake.sh
#else
  #rm -r projects/$FOLDER/sra
  #rm -r projects/$FOLDER/reads
  #rm -r projects/$FOLDER/star
  rm *.o*
  rm *.po*
fi
qsub -N $FOLDER -l h_rt=24:00:0,h_data=1.5G submit.sh $FOLDER
