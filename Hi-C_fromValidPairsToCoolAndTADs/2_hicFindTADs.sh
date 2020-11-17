#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 1
#SBATCH --time 6:00:00
#SBATCH --job-name hicFindTADs
#SBATCH --chdir /scratch/ldelisle/Cool/

path="$PWD/"
bins="5 20 40"

nbOfThreads=1

module purge

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
exists=$(conda info --envs | awk '$1=="hicexplorer3.5.2"{print}' | wc -l)
if [ $exists -ne 1 ]; then
  conda create -y -n hicexplorer3.5.2 hicexplorer=3.5.2 pygenometracks=3.6
fi
conda activate hicexplorer3.5.2

mkdir -p tads
for bin in $bins; do
  for matrix in *_${bin}kb.cool; do
    name=$(basename $matrix .cool)
    size=120
    if [ ! -e tads/${name}.${size}kb_domains.bed ]; then
      echo "finding tads for $name with one window"
      hicFindTADs -m ${matrix} --minBoundaryDistance 100000 --correctForMultipleTesting fdr --outPrefix tads/${name}.${size}kb \
        --minDepth ${size}000 --maxDepth $((2 * ${size}000)) --step $((2 * ${size}000)) &> tads/${name}.${size}kb.log &
    fi
    if [ ! -e tads/${name}.default_domains.bed ]; then
      echo "finding tads for ${name} with default param"
      hicFindTADs -m ${matrix} --minBoundaryDistance 100000 --correctForMultipleTesting fdr --outPrefix tads/${name}.default &> tads/${name}.default.log &
    fi
  done
done
wait
