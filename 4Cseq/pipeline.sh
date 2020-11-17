#! /bin/bash

# Define all paths
gitHubDirectory="/home/ldelisle/softwares/scriptsForBeccariEtAl2021/"
pathForBBCFutils="/home/ldelisle/softwares/" # BBCFutils is found at git clone https://github.com/bbcf/bbcfutils.git -b standalone
genomeDirectory="/scratch/ldelisle/BeccariEtAl2021/genomes/"
FourCDirectory="/scratch/ldelisle/BeccariEtAl2021/4C_likeHTSstation/"
FourCAnalysis="analysisBeccari2021"

# Check that the sra table is available:
if [ ! -e ${gitHubDirectory}/4Cseq/sraTable.txt ]; then
  echo "sraTable.txt does not exists. If the paper is published, please send an email to lucille.delisle@epfl.ch"
  exit 1
fi


# First we need to create all the directories
mkdir -p ${genomeDirectory}/Wt
mkdir -p ${FourCDirectory}

# We need seqtk:
if [ $(which seqtk | wc -l) -ne 1 ]; then
  echo "When looking for seqtk, got"
  which seqtk
  exit 1
fi

# To prepare the genome creation
# Get all numbered chrs + XYM from mm10 ucsc
for i in {1..19} X Y M; do
  wget "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/chr${i}.fa.gz" -O ${genomeDirectory}/chr${i}.fa.gz
done
# Check the download
# Download the md5sums
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/md5sum.txt -O ${genomeDirectory}/md5sum.txt
# Get only the one downloaded
cat ${genomeDirectory}/md5sum.txt | grep -v -P "random|chrUn" > ${genomeDirectory}/md5sum_mychr.txt
# Change dir to do md5sum
currentDir="$PWD"
cd ${genomeDirectory}
n=$(md5sum -c md5sum_mychr.txt | grep -v OK | wc -l)
cd ${currentDir}

# Check it did not failed
if [ $n -gt 0 ]; then
  echo "Download of chr from mm10 failed"
  exit 1
fi

# Concatenate and make them capital letter length 60
for i in {1..19} X Y M; do
  cat ${genomeDirectory}/chr${i}.fa.gz >> ${genomeDirectory}/mm10_UCSC.fa.gz
  echo "chr${i}" >> ${genomeDirectory}/listOfChrs.txt
done
# This can be long:
seqtk seq -U ${genomeDirectory}/mm10_UCSC.fa.gz | seqtk subseq -l 60 - ${genomeDirectory}/listOfChrs.txt > ${genomeDirectory}/Wt/mm10.fa

# You also need to put mm10_rmsk.bed.gz which can be obtained through UCSC website in tools table browser variations and repeat.
if [ ! -e ${genomeDirectory}/mm10_rmsk.bed.gz ]; then
  echo "You don't have mm10_rmsk.bed.gz you need to download it through the table browser in UCSC"
  exit 1
fi

# To prepare the mapseq step:
# Update the init_bein_minilims.py
cp $gitHubDirectory/scripts/init_bein_minilims.py ${FourCDirectory}
sed -i "s#/mnt/BBCF-raw#${FourCDirectory}#g" ${FourCDirectory}/init_bein_minilims.py
sed -i "s#/home/leleu/htsstation#${pathForBBCFutils}#g" ${FourCDirectory}/init_bein_minilims.py

# Begin to launch jobs:
# prepare wt genome for HTSstation
jidMakeGenome=$(sbatch --chdir $genomeDirectory $gitHubDirectory/4Cseq/slurm_scripts/00_makeGenomeFor4CHTS_DpnNla.sh $gitHubDirectory | awk '{print $NF}')

# Get fastq from sra
jidGetFastq=$(sbatch --chdir $FourCDirectory $gitHubDirectory/4Cseq/slurm_scripts/01_getDemultiplexedFastq.sh $FourCAnalysis $gitHubDirectory | awk '{print $NF}')

# Launch the mapping
jidMapping4C=$(sbatch --chdir $FourCDirectory --dependency afterok:${jidMakeGenome},${jidGetFastq} $gitHubDirectory/4Cseq/slurm_scripts/02_mapping.sh $FourCAnalysis $gitHubDirectory $genomeDirectory | awk '{print $NF}')

# Prepare all files for the 4C
jidpre4C=$(sbatch --chdir $FourCDirectory --dependency afterok:${jidMapping4C} $gitHubDirectory/4Cseq/slurm_scripts/03_prepare4C_DpnNla.sh $FourCAnalysis $gitHubDirectory $genomeDirectory | awk '{print $NF}')

# Launch the 4C
sbatch --chdir $FourCDirectory --dependency afterok:${jidpre4C} $gitHubDirectory/4Cseq/slurm_scripts/04_4Cseq.sh $FourCAnalysis $gitHubDirectory
