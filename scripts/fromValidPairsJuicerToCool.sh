# Generate a cool corrected matrix from a validPair from Dekker pipeline
# Require cooler
validPairs=$1
sizeFile=$2
binSize=$3
outputFile=$4
# We create a temporary directory
TEMPORARYDIR=`mktemp -d -t`

cooler csort -i tabix -c1 3 -c2 7 -p1 4 -p2 8 -o ${TEMPORARYDIR}/unique.csort.gz ${validPairs} ${sizeFile}

cooler makebins -o ${TEMPORARYDIR}/bins.txt ${sizeFile} ${binSize}

assembly=$(basename $sizeFile .fa.fai)

cooler cload tabix --assembly ${assembly} -c2 7 -p2 8 ${TEMPORARYDIR}/bins.txt ${TEMPORARYDIR}/unique.csort.gz ${outputFile}

cooler balance --cis-only ${outputFile}
