#This script will create a folder named by the sample to create temporary files and will put the fastq.gz in fromGEO folder with the name of the sample.fastq.gz
pathForSRATable=$1
index=$2
#module purge Will be in the main script
#module load sra-tools/2.8.2-5
SRAnumbers=(`cat $pathForSRATable | awk -v i=$index 'NR==i{print $3}'|tr "," " "`)
SlotNB=(`cat $pathForSRATable | awk -v i=$index 'NR==i{print $4}'|tr "," " "`)
sample=`cat $pathForSRATable | awk -v i=$index 'NR==i{print $2}'`
totNB=`echo ${SlotNB[@]} | awk '{for(i=1;i<=NF;i++){S+=$i}}END{print S}'`
if [ -e fromGEO/${sample}.fastq.gz ]; then
  nbReads=`zcat fromGEO/${sample}.fastq.gz | wc -l | awk '{print $1/4}'`
  if [ $nbReads = $totNB ]; then
    echo "fastq already correct"
    exit 0
  else
    echo "fastq incorrect"
    rm fromGEO/${sample}.fastq.gz
  fi
fi
mkdir -p tmp_$sample
if [ -e tmp_${sample}/failed ]; then
 rm tmp_${sample}/failed
fi
nbSRA=${#SRAnumbers[@]}
for ((i=0; i<$nbSRA;i++))
do
  echo "sraName=${SRAnumbers[$i]}">tmp_${sample}/bash_${i}.sh
  echo "nbSpots=${SlotNB[$i]}">>tmp_${sample}/bash_${i}.sh
  echo "if [ -e tmp_${sample}/\${sraName}.fastq.gz ]; then">>tmp_${sample}/bash_${i}.sh
  echo "  # We first check if the existing fastq is full">>tmp_${sample}/bash_${i}.sh
  echo "  nbSpotsWritten=\`zcat tmp_${sample}/\${sraName}.fastq.gz | wc -l | awk '{print \$1/4}'\`">>tmp_${sample}/bash_${i}.sh
  echo "  if [ ! \$nbSpots = \$nbSpotsWritten ]; then">>tmp_${sample}/bash_${i}.sh
  echo "    echo \"The file present for \$sraName is not correct\"">>tmp_${sample}/bash_${i}.sh
  echo "    rm tmp_${sample}/\${sraName}*">>tmp_${sample}/bash_${i}.sh
  echo "  else">>tmp_${sample}/bash_${i}.sh
  echo "    echo \"The file present for \$sraName is correct\"">>tmp_${sample}/bash_${i}.sh
  echo "    doNotDownload=\"TRUE\"">>tmp_${sample}/bash_${i}.sh
  echo "  fi">>tmp_${sample}/bash_${i}.sh
  echo "fi">>tmp_${sample}/bash_${i}.sh
  echo "if [ -z \$doNotDownload ]; then">>tmp_${sample}/bash_${i}.sh
  echo "  # The quickest way to get fastq from sra is 1) Downlaod sra 2) extract fastq.">>tmp_${sample}/bash_${i}.sh
  echo "  echo \"Will download \$sraName\"">>tmp_${sample}/bash_${i}.sh
  echo "  # Guess URL from the sra name:">>tmp_${sample}/bash_${i}.sh
  echo "  URL=\`echo \$sraName | awk '{print \"ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/\"substr(\$1,1,3)\"/\"substr(\$1,1,6)\"/\"\$1\"/\"\$1\".sra\"}'\`">>tmp_${sample}/bash_${i}.sh
  echo "  if [ -e tmp_${sample}/\${sraName}.sra ]; then">>tmp_${sample}/bash_${i}.sh
  echo "    rm tmp_${sample}/\${sraName}.sra">>tmp_${sample}/bash_${i}.sh
  echo "  fi">>tmp_${sample}/bash_${i}.sh
  echo "  wget \$URL -q -O tmp_${sample}/\${sraName}.sra">>tmp_${sample}/bash_${i}.sh
  echo "  if [ -s tmp_${sample}/\${sraName}.sra ]; then">>tmp_${sample}/bash_${i}.sh
  echo "    # If it worked (the sra file is not empty), we extract the reads">>tmp_${sample}/bash_${i}.sh
  echo "    fastq-dump --log-level fatal --accession \$sraName --defline-seq '@\$sn[_\$rn]/\$ri' --defline-qual '+' --gzip -O tmp_${sample} tmp_${sample}/\${sraName}.sra">>tmp_${sample}/bash_${i}.sh
  echo "  else">>tmp_${sample}/bash_${i}.sh
  echo "    # Else, we need to extract on the flow">>tmp_${sample}/bash_${i}.sh
  echo "    #Need to change home because fastq-dump will download in home/ncbi/public... and then there is no more place in home.">>tmp_${sample}/bash_${i}.sh
  echo "    export HOME=`pwd`">>tmp_${sample}/bash_${i}.sh
  echo "    fastq-dump --log-level fatal --accession \${sraName} --defline-seq '@\$sn[_\$rn]/\$ri' --defline-qual '+' --gzip -O tmp_${sample}">>tmp_${sample}/bash_${i}.sh
  echo "  fi">>tmp_${sample}/bash_${i}.sh
  echo "  nbSpotsWritten=\`zcat tmp_${sample}/\${sraName}.fastq.gz | wc -l | awk '{print \$1/4}'\`">>tmp_${sample}/bash_${i}.sh
  echo "  if [ ! \$nbSpots = \$nbSpotsWritten ]; then">>tmp_${sample}/bash_${i}.sh
  echo "    echo \"The number of spots expected:${nbSpots} is different from what is in the fastq:${nbSpotsWritten}\"">>tmp_${sample}/bash_${i}.sh
  echo "    touch tmp_${sample}/failed">>tmp_${sample}/bash_${i}.sh
  echo "  fi">>tmp_${sample}/bash_${i}.sh
  echo "fi">>tmp_${sample}/bash_${i}.sh
  bash tmp_${sample}/bash_${i}.sh &
done
wait
if [ -e tmp_${sample}/failed ]; then
  exit 1
fi
cat tmp_${sample}/*.fastq.gz > fromGEO/${sample}.fastq.gz
rm -r tmp_${sample}
