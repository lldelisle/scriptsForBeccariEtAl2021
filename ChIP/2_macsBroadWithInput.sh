#!/bin/bash

#SBATCH -o slurm-%A_%2a.out
#SBATCH -e slurm-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 32
#SBATCH --time 8:00:00
#SBATCH --array=1-5
#SBATCH --job-name macsInput
#SBATCH --chdir /scratch/ldelisle/ChIP/

# This script run macs2 with input using extsize 200 broad peak

path="$PWD/"
pathForTable="/home/ldelisle/softwares/scriptsForBeccariEtAl2021/ChIP/tableInputs.txt"
pathForScripts="/home/ldelisle/softwares/scriptsForBeccariEtAl2021/scripts/"

pathForFasta="/home/ldelisle/genomes/fasta/"
genome=mm10

module purge

sample=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
input=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')

pathResults=${path}/${sample}/

if [ ! -e $pathResults ]; then
  echo "$pathResults does not exists"
  exit 1
fi

cd $pathResults

if [ ! -e ${sample}_macs_200_q30_with_input_broad_peaks.broadPeak ]; then
  if [ ! -e ${pathResults}${sample}_mapped_sorted_q30.bam ]; then
    echo "${pathResults}${sample}_mapped_sorted_q30.bam does not exists"
    exit 1
  fi
  if [ ! -e ${path}/${input}/${input}_mapped_sorted_q30.bam ]; then
    echo "${path}/${input}/${input}_mapped_sorted_q30.bam does not exists"
    exit 1
  fi
  macs2 callpeak -t ${pathResults}${sample}_mapped_sorted_q30.bam -c ${path}/${input}/${input}_mapped_sorted_q30.bam -n ${sample}_macs_200_q30_with_input_broad -f BAM --nomodel --extsize 200 --broad 2> ${pathResults}${sample}_macs_200_q30_with_input_broad.log
  cp ${pathResults}${sample}_macs_200_q30_with_input_broad.log ${path}reports/
fi
if [ ! -e ${path}bedGraphs/${sample}_macs_200_q30_norm.bedGraph.gz ]; then
  nreads=$(cat ${pathResults}${sample}_macs_200_q30_with_input_broad.log | grep "tags after filtering in treatment" | awk '{print $NF}')
  tempFile=$(mktemp)
  cat ${sample}_macs_200_q30_with_input_treat_pileup.bdg | awk -v n=$nreads -v OFS="\t" '$4!=0{$4=$4/n*1e6; print}' > $tempFile
  bash ${pathForScripts}/fromMacs2BdgToSimplifiedBdgAndBw.sh ${tempFile} ${path}bedGraphs/${sample}_macs_200_q30_norm "macs2 extended 200bp of ${sample} by million reads" ${pathForFasta}${genome}.fa.fai &
fi
wait

mkdir -p toGEO
cp ${path}bedGraphs/${sample}_macs_200_q30_norm.bw toGEO/
cp ${sample}_macs_200_q30_with_input_broad_peaks.broadPeak toGEO/