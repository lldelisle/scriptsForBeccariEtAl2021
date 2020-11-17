#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 40G 
#SBATCH --cpus-per-task 1
#SBATCH --time 1:00:00
#SBATCH --array=1-4
#SBATCH --job-name Cooler
#SBATCH --chdir /scratch/ldelisle/Cool/

path="$PWD/"
pathForTableWithSamples="/home/ldelisle/softwares/scriptsForBeccariEtAl2021/Hi-C_fromValidPairsToCoolAndTADs/table.txt"
pathForPairs="${path}/"
pathForScripts="/home/ldelisle/softwares/scriptsForBeccariEtAl2021/scripts/"

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
exists=$(conda info --envs | awk '$1=="hicexplorer3.5.2"{print}' | wc -l)
if [ $exists -ne 1 ]; then
  conda create -y -n hicexplorer3.5.2 hicexplorer=3.5.2 pygenometracks=3.6
fi
conda activate hicexplorer3.5.2

sample=$(cat $pathForTableWithSamples | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
validPairs=$(cat $pathForTableWithSamples | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
script=$(cat $pathForTableWithSamples | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $3}')
binsize_kb=$(cat $pathForTableWithSamples | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $4}')
genome=$(cat $pathForTableWithSamples | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $5}')

pathForSizes="/home/ldelisle/genomes/fasta/${genome}.fa.fai"

if [ ! -e $pathForSizes ]; then
  echo "$pathForSizes does not exists"
  exit 1
fi
if [ ! -e ${pathForScripts}/${script} ]; then
  echo "${pathForScripts}/${script} does not exists"
  exit 1
fi

bash ${pathForScripts}/${script} ${pathForPairs}${validPairs} ${pathForSizes} ${binsize_kb}000 ${sample}_${binsize_kb}kb.cool

mkdir toGEO
cp ${sample}_${binsize_kb}kb.cool toGEO/
