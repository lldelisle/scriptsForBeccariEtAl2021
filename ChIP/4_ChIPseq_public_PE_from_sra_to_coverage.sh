#!/bin/bash

#SBATCH -o slurm-%A_%2a.out
#SBATCH -e slurm-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 32
#SBATCH --time 8:00:00
#SBATCH --array=1-2
#SBATCH --job-name ChIP
#SBATCH --chdir /scratch/ldelisle/ChIP/

# This script remove adapters 
# make alignment with Bowtie2
# select MAPQ30 alignments
# Coverage with macs2

path="$PWD/"
pathForTableWithSRA="/home/ldelisle/softwares/scriptsForBeccariEtAl2021/ChIP/sraTable_ChIP_PE.txt"
pathForIndex="/home/ldelisle/genomes/bowtie2/"
pathForFasta="/home/ldelisle/genomes/fasta/"
pathForScripts="/home/ldelisle/softwares/scriptsForBeccariEtAl2021/scripts/"

nbOfThreads=32

module purge
module load gcc/7.4.0 # required for bowtie2, samtools
module load bowtie2/2.3.5
module load samtools/1.9
# cutadapt is version 1.16 is working with python 3.6.1 built with intel 17.0.2
# macs2 is version 2.1.1.20160309 working with python 2


sample=$(cat $pathForTableWithSRA | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
genome=$(cat $pathForTableWithSRA | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $5}')
adapterSeq="TruSeq"

indexPath=${pathForIndex}${genome}

mkdir -p ${path}/$sample

pathResults=${path}/${sample}/
echo $sample
cd $pathResults
if [ ! -e ${path}reports/${sample}_report-cutadapt_PE.txt ]; then
  cd $path
  mkdir -p fromGEO
  module load sra-toolkit/2.9.6
  bash ${pathForScripts}/getPEFastqFromSRAtableAndIndex_fast_parallel_fasterq_secured.sh $pathForTableWithSRA $SLURM_ARRAY_TASK_ID

  fastqR1="${path}/fromGEO/${sample}_R1.fastq.gz"
  fastqR2="${path}/fromGEO/${sample}_R2.fastq.gz"

  if [ ! -s $fastqR1 ]; then
    echo "FASTQ IS EMPTY"
    exit 1
  fi
  cd $pathResults
  if [ $adapterSeq = "TruSeq" ]; then
    cutadapt -j $nbOfThreads -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 30 -m 15 -o ${pathResults}${sample}-cutadapt_R1.fastq.gz -p ${pathResults}${sample}-cutadapt_R2.fastq.gz $fastqR1 $fastqR2 > ${pathResults}${sample}_report-cutadapt_PE.txt
  else
    echo "YOU NEED TO WRITE THE CODE"
    exit 1
  fi
  mkdir -p ${path}reports/
  cp ${pathResults}${sample}_report-cutadapt_PE.txt ${path}reports/
fi

if [ ! -e ${path}reports/${sample}_mapping_stats.txt ];then
  bowtie2 -p $nbOfThreads -x $indexPath -1 ${pathResults}${sample}-cutadapt_R1.fastq.gz -2 ${pathResults}${sample}-cutadapt_R2.fastq.gz 2> ${pathResults}${sample}_mapping_stats.txt  | samtools view --threads $nbOfThreads -Su - | samtools sort --threads $nbOfThreads -o ${pathResults}${sample}_mapped_sorted.bam
  mkdir -p ${path}reports/
  cp ${pathResults}${sample}_mapping_stats.txt ${path}reports/
fi

if [ ! -e ${pathResults}${sample}_mapped_sorted_q30.bam ]; then
  samtools view --threads $nbOfThreads -b ${pathResults}${sample}_mapped_sorted.bam -q 30 > ${pathResults}${sample}_mapped_sorted_q30.bam
fi

mkdir -p ${path}bedGraphs
if [ ! -e ${pathForFasta}${genome}.fa.fai ]; then
  samtools faidx ${pathForFasta}${genome}.fa
fi

if [ ! -e ${path}bedGraphs/${sample}_macs_default_q30_peaks.narrowPeak ]; then
  macs2 callpeak -t ${pathResults}${sample}_mapped_sorted_q30.bam -n ${pathResults}${sample}_macs_default_q30 --call-summits -f BAMPE -B 2> ${pathResults}${sample}_macs_default_q30.log
  cp ${pathResults}${sample}_macs_default_q30.log ${path}reports/
fi
if [ ! -e ${path}bedGraphs/${sample}_macs_default_q30_norm.bedGraph.gz ]; then
  nreads=$(cat ${pathResults}${sample}_macs_default_q30.log | grep "fragments after filtering in treatment" | awk '{print $NF}')
  tempFile=$(mktemp)
  cat ${pathResults}${sample}_macs_default_q30_treat_pileup.bdg | awk -v n=$nreads -v OFS="\t" '$4!=0{$4=$4/n*1e6; print}' > $tempFile
  bash ${pathForScripts}/fromMacs2BdgToSimplifiedBdgAndBw.sh ${tempFile} ${path}bedGraphs/${sample}_macs_default_q30_norm "macs2 default of ${sample} by million reads" ${pathForFasta}${genome}.fa.fai &
fi
wait

mkdir -p toGEO
cp ${path}bedGraphs/${sample}_macs_default_q30_norm.bw toGEO/
cp ${path}bedGraphs/${sample}_macs_default_q30_peaks.narrowPeak toGEO/
