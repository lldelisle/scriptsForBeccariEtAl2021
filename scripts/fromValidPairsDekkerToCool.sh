# Generate a cool corrected matrix from a validPair from Dekker pipeline
# Require cooler
validPairs=$1
sizeFile=$2
binSize=$3
outputFile=$4
# We create a temporary directory
TEMPORARYDIR=`mktemp -d -t`
if [[ ${validPairs} =~ \.gz$ ]]; then
  gunzip -c ${validPairs} > ${TEMPORARYDIR}/input.txt
else
  if [[ ! "${validPairs}" = "/"* ]]; then
    validPairs="$PWD/${validPairs}"
  fi
  ln -s ${validPairs} ${TEMPORARYDIR}/input.txt
fi
# We remove duplicates
awk '{
  if(chrB1 != $2 || chrB2 != $8 || posB1 != $3 || posB2 != $9){
    print
  }
  chrB1=$2
  chrB2=$8
  posB1=$3
  posB2=$9
}' ${TEMPORARYDIR}/input.txt > ${TEMPORARYDIR}/unique.txt

cooler csort -i tabix -c1 2 -c2 8 -p1 3 -p2 9 -o ${TEMPORARYDIR}/unique.csort.gz ${TEMPORARYDIR}/unique.txt ${sizeFile}

cooler makebins -o ${TEMPORARYDIR}/bins.txt ${sizeFile} ${binSize}

assembly=$(basename $sizeFile .fa.fai)

cooler cload tabix --assembly ${assembly} -c2 8 -p2 9 ${TEMPORARYDIR}/bins.txt ${TEMPORARYDIR}/unique.csort.gz ${outputFile}

cooler balance --cis-only ${outputFile}
