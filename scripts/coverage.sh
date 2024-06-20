#!/bin/bash

#SBATCH --job-name=assembly_coverage
#SBATCH -D /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/boreal_moss
#SBATCH -o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/boreal_moss/logs/coverage-%A_%a.slurm.out
#SBATCH --time=48:00:00
#SBATCH --mem=125G
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -A def-ilafores

ml mugqic/bwa mugqic/samtools/1.19.2 mugqic/bedtools/2.30.0 mugqic/seqkit/2.5.0
module list

# parse paths to required files
export SAM_NUM=$(awk "NR==$SLURM_ARRAY_TASK_ID" "${1}")
IFS=$'\t' read -r SAM_ID FQ1 FQ2 ASSEMBLY <<< "$SAM_NUM" # array it up
export SAM_ID FQ1 FQ2 ASSEMBLY
OUT=$(dirname $ASSEMBLY)

if [ ! -d "${OUT}" ]; then
	echo "bin coverage directory not found!"
	exit 1
fi

# Check if script output file is already available, and whether it's empty
if [ -f "$OUT"/average_coverage.txt ]; then 
if [ $(wc -l < "$OUT"/average_coverage.txt) -gt 0 ]; then
	echo "Average coverage already computed!"
	exit 1
else echo "Empty coverage file found, proceeding"
fi
fi

# BWT index
echo "indexing assembly ${SAM_ID}..."
if [[ ! -f "${ASSEMBLY}.bwt" ]]; then
	bwa index "$ASSEMBLY"
else echo "Index found! Skipping."
fi

# Sorting fastqs for bwa mem
echo "Sorting fastqs..."
FQ1_s=${FQ1/.fastq.gz/_sorted.fastq.gz}
FQ2_s=${FQ2/.fastq.gz/_sorted.fastq.gz}
if [[ -f "${FQ1_s}" && -f "${FQ2_s}" ]]; then 
	echo "FASTQs already sorted. Skipping."
else
	seqkit sort -n -j 48 -o "$FQ1_s" "$FQ1"
	seqkit sort -n -j 48 -o "$FQ2_s" "$FQ2"
fi

# Alignemnt
echo "Aligning reads to the assembly..."
if [[ ! -f "$OUT/aligned_reads.sam" ]]; then
	bwa mem -t 48 "$ASSEMBLY" "$FQ1_s" "$FQ2_s" > $OUT/aligned_reads.sam
else echo "Reads already aligned. Skipping..."
fi

#Conversions
echo "Converting and sorting SAM file..."
samtools view -@ 48 -S -b $OUT/aligned_reads.sam -o $OUT/aligned_reads.bam
samtools sort -@ 48 $OUT/aligned_reads.bam -o $OUT/sorted_reads.bam
samtools index $OUT/sorted_reads.bam

#Coverage
echo "Calculating coverage..."
bedtools genomecov -ibam $OUT/sorted_reads.bam > $OUT/assembly_coverage.txt

echo "Calculate average coverage..."
samtools depth $OUT/sorted_reads.bam > $OUT/coverage_per_base.txt
awk '{sum+=$3} END {print sum/NR}' $OUT/coverage_per_base.txt > $OUT/average_coverage.txt

echo "Done!"

# for file in $(find coassembly/assembly/*/coverage -type f -name 'average_coverage.txt'); do mv $file ${file/.txt/_1.txt}; awk -F' ' '{print $3}' ${file/.txt/_1.txt} > $file; done