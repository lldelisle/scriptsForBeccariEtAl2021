#!/bin/bash

#SBATCH -o slurm-%A_%2a.out
#SBATCH -e slurm-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 30G
#SBATCH --cpus-per-task 36
#SBATCH --time 8:00:00
#SBATCH --array=1-7
#SBATCH --job-name ChIP
#SBATCH --chdir /scratch/ldelisle/ChIP/

# This script remove adapters 
# make alignment with Bowtie2
# select MAPQ30 alignments

path="$PWD/"
pathForTableWithSRA="/home/ldelisle/softwares/scriptsForBeccariEtAl2021/ChIP/sraTable_ChIP.txt"
pathForFastq="${path}/fromGEO/"
pathForIndex="/home/ldelisle/genomes/bowtie2/"
pathForScripts="/home/ldelisle/softwares/scriptsForBeccariEtAl2021/scripts/"

nbOfThreads=36

module purge
module load gcc/7.4.0 # required for bowtie2, samtools, star and bedtools
module load bowtie2/2.3.5
module load samtools/1.9 
module load picard/2.19.0
module load bedtools2/2.27.1
# cutadapt is version 1.16 is working with python 3.6.1 built with intel 17.0.2
# macs2 is version 2.1.1.20160309 working with python 2


sample=$(cat $pathForTableWithSRA | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
genome=$(cat $pathForTableWithSRA | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $5}')
adapterSeq="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
fastqFile=${sample}.fastq.gz

indexPath=${pathForIndex}${genome}

mkdir -p ${path}/$sample

pathResults=${path}/${sample}/
echo $sample
cd $pathResults
if [ ! -e ${path}reports/${sample}_report-cutadapt.txt ]; then
  cd $path
  mkdir -p fromGEO
  module load sra-toolkit/2.9.6
  bash ${pathForScripts}/getFastqFromSRAtableAndIndex_fast_parallel_secured.sh $pathForTableWithSRA $SLURM_ARRAY_TASK_ID
  cd $pathResults

  if [ -z $adapterSeq ]; then
    cutadapt -m 15 -j $nbOfThreads -q 30 -o ${pathResults}${sample}-cutadapt.fastq.gz "${pathForFastq}${fastqFile}" > ${pathResults}${sample}_report-cutadapt.txt
  else
    cutadapt -m 15 -j $nbOfThreads -a $adapterSeq -q 30 -o ${pathResults}${sample}-cutadapt.fastq.gz "${pathForFastq}${fastqFile}" > ${pathResults}${sample}_report-cutadapt.txt
  fi
  mkdir -p ${path}reports/
  cp ${pathResults}${sample}_report-cutadapt.txt ${path}reports/
fi

if [ ! -e ${path}reports/${sample}_mapping_stats.txt ];then
  bowtie2 -p $nbOfThreads -x $indexPath -U ${sample}-cutadapt.fastq.gz 2> ${pathResults}${sample}_mapping_stats.txt  | samtools view --threads $nbOfThreads -Su - | samtools sort --threads $nbOfThreads -o ${pathResults}${sample}_mapped_sorted.bam
  mkdir -p ${path}reports/
  cp ${pathResults}${sample}_mapping_stats.txt ${path}reports/
fi

if [ ! -e ${pathResults}${sample}_mapped_sorted_q30.bam ]; then
  samtools view --threads $nbOfThreads -b ${pathResults}${sample}_mapped_sorted.bam -q 30 > ${pathResults}${sample}_mapped_sorted_q30.bam
fi

wait


