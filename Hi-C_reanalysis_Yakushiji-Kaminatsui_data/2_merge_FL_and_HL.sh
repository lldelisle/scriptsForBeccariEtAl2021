#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 4G 
#SBATCH --cpus-per-task 1
#SBATCH --time 30:00
#SBATCH --job-name merge
#SBATCH --chdir /scratch/ldelisle/HiCUP

for sample in HH20_HL_WT_CHi-C HH20_FL_WT_CHi-C; do
  cat ${sample}/${sample}.validPairs_nonSorted.txt.gz >> HH20_HL_FL_WT_CHi-C.validPairs_nonSorted.txt.gz
done
