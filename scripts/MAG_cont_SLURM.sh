#!/bin/bash -l

#SBATCH --mail-type=END,FAIL
#SBATCH -D /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/boreal_moss/
#SBATCH -o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/boreal_moss/MAG_analysis/contamination_check/logs/bowtie-%A_%a.slurm.out
#SBATCH --time=96:00:00
#SBATCH -N 1
#SBATCH -A def-ilafores
#SBATCH -J bowtie2

# fetch and index reference genomes 
ml mugqic/bowtie2

export MAG_DIR="${1}"
export MAG_LIST="${2}"
export SAMPLE_LIST="${3}"
export HOST_IDX="${4}"
export MAG=$(awk "NR==$SLURM_ARRAY_TASK_ID" ${MAG_LIST} | cut -d, -f1)
export SOURMASH="singularity exec -e -B $ANCHOR/home:$ANCHOR/home $ANCHOR/home/def-ilafores/programs/ILL_pipelines/containers/sourmash.4.7.0.sif sourmash"
export OUT=$PWD/MAG_analysis/contamination_check/${MAG}

# Index reference genome if necessary
idx_count=$(ls "$MAG_DIR/${MAG}_idx"/*.bt2 2> /dev/null | wc -l)
if [[ "$idx_count" -ne 6 ]]; then
	echo "Missing genome index files! Exiting."
	exit 1
fi

# initiate result file:
:> $OUT/contamination_by_sample.txt

while read -r SAMPLE FQ_P1 FQ_P2; do
	
	# Align sample reads
	mkdir -p $OUT/${SAMPLE}

if [[ ! -f $OUT/$SAMPLE.aligned_paired.1.fastq.gz ]] || [[ ! -f $OUT/$SAMPLE.aligned_paired.2.fastq.gz ]]; then
	echo "Aligning $SAMPLE reads to $MAG MAG..."	
	bowtie2 --very-sensitive --no-unal -x $MAG_DIR/${MAG}_idx/${MAG} -1 $FQ_P1* -2 $FQ_P2* \
	--al-conc-gz $OUT/$SAMPLE.aligned_paired.%.fastq.gz \
	--threads $SLURM_NTASKS -S $OUT/$SAMPLE.sam 

else echo "Alignments of $SAMPLE reads to $MAG found, skipping!"
fi

if [[ ! -f $OUT/$SAMPLE.aligned.${MAG}.1.fastq.gz ]] || [[ ! -f $OUT/$SAMPLE.aligned.${MAG}.2.fastq.gz ]]; then
	
	echo "Aligning $MAG reads to host genome..."
	bowtie2 --very-sensitive --no-unal -x $HOST_IDX -1 $OUT/$SAMPLE.aligned_paired.1.fastq.gz -2 $OUT/$SAMPLE.aligned_paired.2.fastq.gz \
 	--al-conc-gz $OUT/$SAMPLE.aligned.${MAG}.%.fastq.gz --threads $SLURM_NTASKS -S $OUT/${SAMPLE}_${MAG}.sam
	
else echo "Alignments of $MAG reads to host genome found, skipping..."
fi

#Number of reads mapping to MAG
numReads_MAG=$(wc -l < $(zcat $OUT/$SAMPLE.aligned_paired.1.fastq.gz))

#of which number of reads also mapping to host genome
numReads_HOSt=$(wc -l < $(zcat $OUT/$SAMPLE.aligned.${MAG}.1.fastq.gz))

# Compute the ratio using awk
ratio=$(awk "BEGIN {print $numReads_MAG / $numReads_HOSt}")

echo -e "$MAG\t$SAMPLE\t$ratio" >> $OUT/contamination_by_sample.txt

done < $SAMPLE_LIST
