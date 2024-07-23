#!/bin/bash -l

#SBATCH --mail-type=END,FAIL
#SBATCH -D /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/boreal_moss/
#SBATCH -o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/boreal_moss/MAG_analysis/contamination_check/logs/bowtie-%A_%a.slurm.out
#SBATCH --time=48:00:00
#SBATCH --mem=30G
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -A def-ilafores
#SBATCH -J bowtie2

# fetch and index reference genomes 
module load bowtie2/2.5.1

export MAG_DIR="${1}"
export MAG_NAMES="${2}"
export MAG=$(awk "NR==$SLURM_ARRAY_TASK_ID" ${MAG_NAMES} | cut -d, -f1)
export OUT=MAG_analysis/contamination_check

# Index reference genome if necessary
idx_count=$(ls "$MAG_DIR/${MAG}_idx"/*.bt2 2> /dev/null | wc -l)
if [[ "$idx_count" -ne 6 ]]; then
	echo "Missing genome index files! Exiting."
	exit 1
fi

# Align sample reads
bowtie2 --very-sensitive -x $MAG_DIR/${MAG}_idx/${MAG} -1 $FQ_P1 -2 $FQ_P2 -U $FQ_U1,$FQ_U2 \
	--al-gz $OUT/$MAG.aligned_unpaired.fastq.gz --al-conc-gz $OUT/$MAG.aligned_paired.%.fastq.gz \
		--threads $SLURM_NTASKS -S $OUT/$MAG.aligned.sam