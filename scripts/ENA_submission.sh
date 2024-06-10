#!/usr/bin/env bash

MOSS=/home/def-ilafores/analysis/boreal_moss
cd $MOSS

mkdir -p ENA_submission/manifests

while read -r sample Host; do 
file="ENA_submission/manifests/${sample}_manifest.txt"
touch ${file}
echo -e "STUDY\tPRJEB76464" >> "$file" 
echo -e "SAMPLE\t${sample}" >> "$file"
echo -e "NAME\tBoreal Moss Microbiome" >> "$file"
echo -e "INSTRUMENT\tIllumina NovaSeq 6000" >> "$file"
echo -e "INSERT_SIZE\t150" >> "$file"
echo -e "LIBRARY_SOURCE\tMETAGENOMIC" >> "$file"
echo -e "LIBRARY_SELECTION\tRANDOM" >> "$file"
echo -e "LIBRARY_STRATEGY\tWGS" >> "$file"
for fastq in preproc/"$sample"/*.fastq; do
	echo -e "FASTQ\t${fastq}" >> "$file"
done
done < <(tail -n +2 comptype.tsv)
