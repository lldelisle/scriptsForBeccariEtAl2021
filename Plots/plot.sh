#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 30G
#SBATCH --cpus-per-task 1
#SBATCH --time 8:00:00
#SBATCH --job-name ChIP
#SBATCH --chdir /scratch/ldelisle/Plots/

pathForHiC="/scratch/ldelisle/Cool/"
pathForAnnot="$PWD/"
pathFor4Cbedgraph="/scratch/ldelisle/BeccariEtAl2021/4C_likeHTSstation/analysisBeccari2021/toGEO/"
pathForChIP="/scratch/ldelisle/ChIP/"
pathForScripts="/home/ldelisle/softwares/scriptsForBeccariEtAl2021/scripts/"

# Get the last ensembl annotation for mouse
wget ftp://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.gtf.gz

# Select genes in the region protein coding:
zcat Mus_musculus.GRCm38.101.gtf.gz | awk -F "\t" '
$1=="15" && $4 > 94200000 && $5 < 97000000 && # Only entries in the region
$9~/protein_coding/{ # Only protein_coding genes
  print
}' > ${pathForAnnot}/genes.gtf

# Select genes protein coding and convert to bed9:
cat genes.gtf | awk -F "\t" '
$3=="gene"{ # Only entries describing the gene
  split($9,a,"\"") # Cut the metadata at quote sign
  for(i in a){
    if(a[i]=="; gene_name " || a[i]=="gene_name "){
      name=a[i+1] # The gene name is after gene_name "
    }
  }
  rgb="204,204,204"
  if (name == "Dbx2"){
    rgb="0,0,102"
  }
  print $1"\t"$4 - 1"\t"$5"\t"name"\t"0"\t"$7"\t"$4 - 1"\t"$4 - 1"\t"rgb # Print in bed9 format 0-based
}' > ${pathForAnnot}/genes.bed9

# Get the last ensembl annotation for chicken galGal6
wget ftp://ftp.ensembl.org/pub/release-101/gtf/gallus_gallus/Gallus_gallus.GRCg6a.101.gtf.gz

# Select genes in the region protein coding:
zcat Gallus_gallus.GRCg6a.101.gtf.gz | awk -F "\t" '
$1=="1" && $4 > 30000000 && $5 < 32000000 && # Only entries in the region
$9~/protein_coding/{ # Only protein_coding genes
  print
}' > ${pathForAnnot}/genes_galGal6.gtf

# Select genes protein coding and convert to bed9:
cat ${pathForAnnot}/genes_galGal6.gtf | awk -F "\t" '
$3=="gene"{ # Only entries describing the gene
  split($9,a,"\"") # Cut the metadata at quote sign
  for(i in a){
    if(a[i]=="; gene_name " || a[i]=="gene_name "){
      name=a[i+1] # The gene name is after gene_name "
    }
  }
  rgb="204,204,204"
  if (name == "Dbx2"){
    rgb="0,0,102"
  }
  print $1"\t"$4 - 1"\t"$5"\t"name"\t"0"\t"$7"\t"$4 - 1"\t"$4 - 1"\t"rgb # Print in bed9 format 0-based
}' > ${pathForAnnot}/genes_galGal6.bed9


# Write bed files for annotated regions:
echo -e "chr15\t95600982\t95601996\t1
chr15\t95637600\t95639316\t2
chr15\t95647795\t95648894\t3" > ${pathForAnnot}/DLE.bed

echo -e "chr15\t95638744\t95641449\tVISTAmm1571" > ${pathForAnnot}/VISTA.bed

# Filter the narrowPeak to keep only the best summit:
for target in Hoxa13 Hoxd13; do
  bash ${pathForScripts}/filterNarrowPeak.sh ${pathForChIP}toGEO/DFL_E11_${target}_macs_default_q30_peaks.narrowPeak ${pathForAnnot}/DFL_E11_${target}_macs_default_q30_peaks_filtered.narrowPeak
done

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
exists=$(conda info --envs | awk '$1=="hicexplorer3.5.2"{print}' | wc -l)
if [ $exists -ne 1 ]; then
  conda env create -n hicexplorer3.5.2 hicexplorer=3.5.2 pygenometracks=3.6
fi
conda activate hicexplorer3.5.2

# Generation of ini files used in multiple figures
ini_file=Dbx2_4C.ini
echo "[4C Dbx2 DFL]
file = ${pathFor4Cbedgraph}Dbx2_E12_DFL_neq3.bedGraph.gz
title = 4C Dbx2
color = #3298cb
use_middle = true
min_value = 0
max_value = 20
alpha = 0.8
height = 3

[4C Dbx2 PFL]
file = ${pathFor4Cbedgraph}Dbx2_E12_PFL_neq3.bedGraph.gz
color = #346733
use_middle = true
alpha = 0.8
overlay_previous = share-y
show_data_range = false

[4C Dbx2 PFL]
file = ${pathFor4Cbedgraph}Dbx2_E12_PFL_neq3.bedGraph.gz
second_file = ${pathFor4Cbedgraph}Dbx2_E12_DFL_neq3.bedGraph.gz
operation = min(file, second_file)
color = #000099
use_middle = true
alpha = 1
overlay_previous = share-y
show_data_range = false
" > ${ini_file}

ini_file=acethylation_chip.ini
echo "[H3K27ac DFL]
file = ${pathForChIP}toGEO/DFL_E12_WT_H3K27ac_neq2.bw
title = E12 H3K27ac
color = #3298cb
min_value = 0
max_value = 1.2
number_of_bins = 1500
alpha = 0.8
height = 3

[H3K27ac PFL]
file = ${pathForChIP}toGEO/PFL_E12_WT_H3K27ac_rep1_macs_200_q30_norm.bw
color = #346733
number_of_bins = 1500
alpha = 0.8
overlay_previous = share-y
show_data_range = false

[H3K27ac DFL both]
file = ${pathForChIP}toGEO/DFL_E12_WT_H3K27ac_neq2.bw
second_file = ${pathForChIP}toGEO/PFL_E12_WT_H3K27ac_rep1_macs_200_q30_norm.bw
title = E12 H3K27ac
color = #0c0c78
overlay_previous = share-y
operation = min(file, second_file)
number_of_bins = 1500
alpha = 1
" > ${ini_file}

for sample in DFL_E12_WT_H3K27ac_rep1 DFL_E12_WT_H3K27ac_rep2 PFL_E12_WT_H3K27ac_rep1; do
  echo "[broad peaks $sample]
file = ${pathForChIP}toGEO/${sample}_macs_200_q30_with_input_broad_peaks.broadPeak
title = $sample
file_type = bed
labels = false
display = collapsed
" >> ${ini_file}
  if [[ $sample = *"DFL"* ]]; then
    echo "color = #3298cb
border_color = #3298cb
" >> ${ini_file}
  else
    echo "color = #346733
border_color = #346733
" >> ${ini_file}
  fi
done

# For each figure the ini file is generated followed by the plot command
ini_file="fig2top.ini"
echo "[scalebar]
where = top
x_center = 96800000
size = 250000
file_type = scalebar

[tads]
file = ${pathForHiC}toGEO/ES_chr15_5kb.120kb_domains.bed
color = #66cc99
border_color = none
display = interleaved
labels = false

[spacer]
height = 0.2

[Hi-C ES]
file = ${pathForHiC}toGEO/ES_chr15_5kb.cool
min_value = 0
max_value = 0.01
depth = 1000000

[insulation score]
file = ${pathForHiC}tads/ES_chr15_5kb.120kb_tad_score.bm
type = line
file_type = bedgraph
height = 3

[spacer]

[genes bed9]
file = ${pathForAnnot}genes.bed9
title = genes ensembl
color = bed_rgb
display = interleaved
file_type = bed
arrowhead_included = true
height = 1.5

[spacer]
" > ${ini_file}

cat Dbx2_4C.ini >> ${ini_file}

pgt --tracks ${ini_file} --region chr15:94200000-97000000 -o ${ini_file/ini/pdf} --dpi 500

ini_file="fig2bottom.ini"
echo "[scalebar]
where = top
x_center = 95787000
size = 10000
file_type = scalebar
" > $ini_file
cat Dbx2_4C.ini >> ${ini_file}
echo "[spacer]
height = 0.2

[genes]
file = ${pathForAnnot}genes.gtf
height = 1
color = #0c0c78
prefered_name = gene_name
merge_transcripts = true
style = UCSC
all_labels_inside = true

[spacer]
" >> ${ini_file}
for sample in ES_chr15_5kb E14_Cortex_chr15_5kb E12_DL_40kb; do
  echo "[tads]
file = ${pathForHiC}toGEO/${sample}.120kb_domains.bed
title = $sample
color = #99cccc
border_color = #999999
display = interleaved
labels = false
height = 0.8

[spacer]
height = 0.2
" >> ${ini_file}
done
cat acethylation_chip.ini >> ${ini_file}
for target in Hoxa13 Hoxd13; do
  echo "[${target} E11 DFL]
file = ${pathForChIP}toGEO/DFL_E11_${target}_macs_default_q30_norm.bw
title = DFL_E11_${target}
color = #003366
summary_method = max
number_of_bins = 1500
min_value = 0
max_value = 3.1
height = 3

[peak ${target} E11 DFL]
file = ${pathForAnnot}/DFL_E11_${target}_macs_default_q30_peaks_filtered.narrowPeak
title = peaks (no input)
file_type = bed
labels = false
display = collapsed
color = binary
line_width = 0.1
min_value = 0
" >> ${ini_file}
done
echo "[DLE]
file = ${pathForAnnot}/DLE.bed
title = DLE
color = #990000
display = collapsed
border_color = none

[VISTA]
file = ${pathForAnnot}/VISTA.bed
title = VISTA
color = #0c0055
border_color = none
labels = false

">> ${ini_file}


pgt --tracks ${ini_file} --region chr15:95,520,000-95,797,000 -o ${ini_file/ini/pdf}


ini_file="fig6D.ini"
echo "[scalebar]
where = top
x_center = 96338000
size = 125000
file_type = scalebar
" > $ini_file
for sample in ES_chr15_5kb E14_Cortex_chr15_5kb E12_DL_40kb; do
  echo "[tads]
file = ${pathForHiC}toGEO/${sample}.120kb_domains.bed
title = $sample
color = #99cccc
border_color = #999999
display = interleaved
labels = false
height = 0.8

[spacer]
height = 0.2
" >> ${ini_file}
done
cat Dbx2_4C.ini >> ${ini_file}
for viewpoint in Nell2 Ano6; do
  echo "[spacer]
height = 0.2

[4C ${viewpoint} DFL]
file = ${pathFor4Cbedgraph}${viewpoint}_E12_DFL.bedGraph.gz
title = 4C ${viewpoint}
color = #3298cb
use_middle = true
min_value = 0
max_value = 30
alpha = 1
height = 3
" >> ${ini_file}
done
echo "[spacer]
height = 0.2

[genes bed9]
file = ${pathForAnnot}genes.bed9
title = genes ensembl
color = bed_rgb
display = interleaved
file_type = bed
arrowhead_included = true
height = 1.5
" >> ${ini_file}
for viewpoint in DLE1 DLE2; do
  echo "[spacer]
height = 0.2

[4C ${viewpoint} DFL]
file = ${pathFor4Cbedgraph}${viewpoint}_E12_DFL.bedGraph.gz
title = 4C ${viewpoint}
color = #3298cb
use_middle = true
min_value = 0
max_value = 20
alpha = 1
height = 3
" >> ${ini_file}
done
echo "[spacer]
height = 0.2

[DLE]
file = ${pathForAnnot}/DLE.bed
title = DLE
color = #990000
display = collapsed
labels = false
border_color = none

[VISTA]
file = ${pathForAnnot}/VISTA.bed
title = VISTA
color = #0c0055
border_color = none
labels = false
overlay_previous = yes
" >> ${ini_file}
out_file=$(basename $ini_file .ini).pdf
pgt --tracks ${ini_file} --region chr15:94,600,000-96,538,000 -o $out_file

ini_file="fig7A.ini"

echo "[scalebar]
where = top
x_center = 31200000
size = 100000
file_type = scalebar

[spacer]
height = 0.1

[Hi-C HH20 limb]
file = ${pathForHiC}toGEO/HH20_HL_FL_20kb.cool
min_value = 0
max_value = 0.03
depth = 1000000
" > ${ini_file}
for param in default 120kb; do
  echo "[insulation score $param]
file = ${pathForHiC}tads/HH20_HL_FL_20kb.${param}_tad_score.bm
height = 3
" >> ${ini_file}
  if [ $param = default ]; then
    echo "type = lines
" >> ${ini_file}
  else
    echo "type = line
file_type = bedgraph
" >> ${ini_file}
  fi
  echo "[tads $param]
file = ${pathForHiC}toGEO/HH20_HL_FL_20kb.${param}_domains.bed
display = interleaved
labels = false
color = black
" >> ${ini_file}
done
echo "[spacer]

[genes bed9]
file = ${pathForAnnot}genes_galGal6.bed9
title = genes ensembl 101 (galGal6)
color = bed_rgb
display = interleaved
file_type = bed
arrowhead_included = true
height = 1.5
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr1:30,000,000-31,262,500 -o ${ini_file/ini/pdf}

ini_file="figS2.ini"
echo "[scalebar]
where = top
x_center = 96800000
size = 250000
file_type = scalebar

[spacer]

[genes bed9]
file = ${pathForAnnot}genes.bed9
title = genes ensembl
color = bed_rgb
display = interleaved
file_type = bed
arrowhead_included = true
height = 1.5

[spacer]
" > ${ini_file}
for sample in ES_chr15_5kb E14_Cortex_chr15_5kb E12_DL_40kb; do
  for param in default 120kb; do
    echo "[tads $sample $param]
file = ${pathForHiC}toGEO/${sample}.${param}_domains.bed
color = black
title = tads $param
display = interleaved
labels = false
alpha = 0.8
height = 0.8
border_color = none
" >> ${ini_file}
    if [ $param = "default" ]; then
      echo "color = #cdcdcd
" >> ${ini_file}
    else
      echo "color = #66cd99
" >> ${ini_file}
    fi
    echo "[spacer]
height = 0.2
" >> ${ini_file}
  done
  echo "[Hi-C $sample]
file = ${pathForHiC}toGEO/${sample}.cool
title = $sample
min_value = 0
depth = 1000000
" >> ${ini_file}
  if [[ $sample = *"5kb" ]]; then
    echo "max_value = 0.01
" >> ${ini_file}
  else
    echo "max_value = 0.05
" >> ${ini_file}
  fi
  echo "[insulation score $sample]
file = ${pathForHiC}tads/${sample}.default_tad_score.bm
title = tad score $sample
height = 3
type = lines
" >> ${ini_file}
done

pgt --tracks ${ini_file} --region chr15:94200000-97000000 -o ${ini_file/ini/pdf} --dpi 500

ini_file="figS3AB.ini"
echo "[scalebar]
where = top
x_center = 96050000
size = 50000
file_type = scalebar

[tads]
file = ${pathForHiC}toGEO/ES_chr15_5kb.120kb_domains.bed
color = #66cc99
border_color = none
display = interleaved
labels = false

[spacer]
height = 0.2
" > ${ini_file}
cat Dbx2_4C.ini >> ${ini_file}
echo "[spacer]
height = 0.2

[genes]
file = ${pathForAnnot}genes.gtf
height = 1
color = #0c0c78
prefered_name = gene_name
merge_transcripts = true
style = UCSC
all_labels_inside = true

[spacer]
height = 0.2
" >> ${ini_file}
cat acethylation_chip.ini >> ${ini_file}
for target in Hoxa13 Hoxd13; do
  echo "[${target} E11 DFL]
file = ${pathForChIP}toGEO/DFL_E11_${target}_macs_default_q30_norm.bw
title = DFL_E11_${target}
color = #003366
summary_method = max
number_of_bins = 1500
min_value = 0
max_value = 3.1
height = 3
" >> ${ini_file}
done
echo "[spacer]
height = 0.2

[H3K27me3 DFL]
file = ${pathForChIP}toGEO/PFL_E12_WT_H3K27me3_rep1_macs_200_q30_norm.bw
title = E12 H3K27me3
color = #346733
min_value = 0
number_of_bins = 1500
alpha = 0.8
height = 3

[H3K27me3 PFL]
file = ${pathForChIP}toGEO/DFL_E12_WT_H3K27me3_rep2_macs_200_q30_norm.bw
color = #3298cb
number_of_bins = 1500
alpha = 0.8
overlay_previous = share-y
show_data_range = false

[H3K27me3 DFL both]
file = ${pathForChIP}toGEO/PFL_E12_WT_H3K27me3_rep1_macs_200_q30_norm.bw
second_file = ${pathForChIP}toGEO/DFL_E12_WT_H3K27me3_rep2_macs_200_q30_norm.bw
color = #0c0c78
overlay_previous = share-y
operation = min(file, second_file)
number_of_bins = 1500
alpha = 1

[spacer]
height = 0.2

[DLE]
file = ${pathForAnnot}/DLE.bed
title = DLE
color = #990000
display = collapsed
labels = false
border_color = none

[VISTA]
file = ${pathForAnnot}/VISTA.bed
title = VISTA
color = #0c0055
border_color = none
labels = false
overlay_previous = yes
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr15:95209000-96097000 -o ${ini_file/ini/pdf}

ini_file="figS3D.ini"
echo "[scalebar]
where = top
x_center = 96050000
size = 50000
file_type = scalebar

[tads]
file = ${pathForHiC}toGEO/ES_chr15_5kb.120kb_domains.bed
color = #66cc99
border_color = none
display = interleaved
labels = false

[spacer]
height = 0.2
" > $ini_file
cat Dbx2_4C.ini >> ${ini_file}
for viewpoint in Nell2 Ano6; do
  echo "[spacer]
height = 0.2

[4C ${viewpoint} DFL]
file = ${pathFor4Cbedgraph}${viewpoint}_E12_DFL.bedGraph.gz
title = 4C ${viewpoint}
color = #3298cb
use_middle = true
min_value = 0
max_value = 30
alpha = 1
height = 3
" >> ${ini_file}
done
echo "
[spacer]
height = 0.2

[genes gtf]
file = ${pathForAnnot}genes.gtf
height = 1
color = #0c0c78
prefered_name = gene_name
merge_transcripts = true
style = UCSC
all_labels_inside = true
" >> ${ini_file}
for viewpoint in DLE1 DLE2; do
  echo "[spacer]
height = 0.2

[4C ${viewpoint} DFL]
file = ${pathFor4Cbedgraph}${viewpoint}_E12_DFL.bedGraph.gz
title = 4C ${viewpoint}
color = #3298cb
use_middle = true
min_value = 0
max_value = 20
alpha = 1
height = 3
" >> ${ini_file}
done
echo "[spacer]
height = 0.2

[DLE]
file = ${pathForAnnot}/DLE.bed
title = DLE
color = #990000
display = collapsed
labels = false
border_color = none

[VISTA]
file = ${pathForAnnot}/VISTA.bed
title = VISTA
color = #0c0055
border_color = none
labels = false
overlay_previous = yes
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr15:95209000-96097000 -o ${ini_file/ini/pdf} --plotWidth 29.4


cat figS3D.ini | awk '{gsub("96050000", "95787000", $0); gsub("size = 50000", "size = 10000", $0); print}' > figS3E.ini
pgt --tracks figS3E.ini --region chr15:95520000-95797000 -o figS3E.pdf --plotWidth 15.81
