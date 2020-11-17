#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 20G 
#SBATCH --cpus-per-task 28
#SBATCH --time 72:00:00
#SBATCH --array=1-2
#SBATCH --job-name HiCUP
#SBATCH --chdir /scratch/ldelisle/HiCUP

path="$PWD/"
pathForScripts="/home/ldelisle/softwares/scriptsForBeccariEtAl2021/scripts/"
pathForTableWithSRA="/home/ldelisle/softwares/scriptsForBeccariEtAl2021/Hi-C_reanalysis_Yakushiji-Kaminatsui_data/table.txt"
genome="galGal6"
restSeq="^GATC"
restName="DpnII"
pathForB2Index="/home/ldelisle/genomes/bowtie2/$genome"
pathForFasta="/home/ldelisle/genomes/fasta/${genome}.fa"
pathForHiCUP="/home/ldelisle/softwares/HiCUP-0.8.0/"
pathForSizes="${pathForFasta}.fai"


module purge
module load gcc/7.4.0
module load openblas/0.3.6-openmp
module load r/3.6.0
module load bowtie2/2.3.5
module load samtools/1.9
module load htslib/1.9
pathForR=$(which R)
pathForBowtie2=$(which bowtie2)

pathForDigest="${path}/${genome}_digester_$restName.txt.gz"

#Digest 
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
  # install dependencies for HiCUP >= v0.8
  echo "if (!\"devtools\" %in% installed.packages()){
  install.packages(\"devtools\", repos = \"https://stat.ethz.ch/CRAN/\")
}
devtools::install_github(\"lldelisle/usefulLDfunctions\")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor(\"tidyverse\")
safelyLoadAPackageInCRANorBioconductor(\"plotly\")
" > dependencies.R
  Rscript dependencies.R
  
  # produce the digester file
  if [ ! -e $pathForDigest ]; then
    ${pathForHiCUP}hicup_digester --re1 $restSeq,$restName --genome $genome --zip --outdir $path $pathForFasta
    mv ${path}/Digest_${genome}_${restName}* $pathForDigest
  fi
fi

# Get the sample name from the table
sample=$(cat $pathForTableWithSRA | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')

# Download the fastq
mkdir -p fromGEO
module load sra-toolkit/2.9.6
bash ${pathForScripts}/getPEFastqFromSRAtableAndIndex_fast_parallel_fasterq_secured.sh $pathForTableWithSRA $SLURM_ARRAY_TASK_ID

fastqFileR1="${path}/fromGEO/${sample}_R1.fastq.gz"
fastqFileR2="${path}/fromGEO/${sample}_R2.fastq.gz"

mkdir -p $sample
cd $sample

# Check if an output bam exists
inputBAM=$(find . -name "*.hicup.bam")

if [ -z $inputBAM ]; then
  # Run hicup
  ${pathForHiCUP}hicup --bowtie2 $pathForBowtie2 --digest $pathForDigest --format Sanger --index $pathForB2Index --keep --threads 28 --zip --r $pathForR $fastqFileR1 $fastqFileR2
  # Update the inputBAM variable
  inputBAM=$(find . -name "*.hicup.bam")
fi

pathForDigestNGZ="${path}/${genome}_digester_$restName.txt"

if [ ! -e $pathForDigestNGZ ]; then
  gunzip -c $pathForDigest > $pathForDigestNGZ
fi

# Convert the bam to validPair in juicebox format
# The python script contrary to the provided converter
# Keep the fragment id and use the middle of the fragment
if [ ! -e ${sample}.validPairs_nonSorted.txt.gz ]; then
  python ${pathForScripts}/fromHicupToJuicebox_withTempFolder.py  --fragmentFile $pathForDigestNGZ --colForChr 1 --colForStart 2 --colForEnd 3 --colForID 4 --lineToSkipInFragmentFile 2 --useMid $inputBAM | gzip > ${sample}.validPairs_nonSorted.txt.gz
fi

inputValidPairs=${sample}.validPairs_nonSorted.txt.gz
pairs=${sample}.validPairs.csort.gz

if [ ! -e $pathForSizes ]; then
  samtools faidx $pathForFasta
fi

if [ ! -e $pairs ]; then
  # sort and index the pairs with cooler and tabix
  cooler csort -i tabix -c1 3 -c2 7 -p1 4 -p2 8 -o $pairs $inputValidPairs $pathForSizes
fi

mkdir -p toGEO
cp $pairs toGEO/