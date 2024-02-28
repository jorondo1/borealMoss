#!/bin/bash -l

#SBATCH --mail-type=END,FAIL
#SBATCH -D /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/boreal_moss/SM_abund
#SBATCH -o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/boreal_moss/SM_abund/logs/sourmash-%A_%a.slurm.out
#SBATCH --time=24:00:00
#SBATCH --mem=250G
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -A def-ilafores
#SBATCH -J sourmash

tmp=$PWD/temp
mkdir -p $tmp

newgrp def-ilafores
export ILAFORES=/nfs3_ib/nfs-ip34/home/def-ilafores
export PARENT_DIR=${ILAFORES}/analysis/boreal_moss
export ILL_PIPELINES=${ILAFORES}/programs/ILL_pipelines
export sourmash="singularity exec --writable-tmpfs -e -B $tmp:$tmp -B /nfs3_ib/nfs-ip34/fast:/nfs3_ib/nfs-ip34/fast -B /nfs3_ib/nfs-ip34/home:/nfs3_ib/nfs-ip34/home ${tmp}/sourmash.4.7.0.sif sourmash"

echo "loading env"
module load StdEnv/2020 apptainer/1.1.5

# copy container
cp $ILL_PIPELINES/containers/sourmash.4.7.0.sif $tmp/

export __sample_line=$(cat $PARENT_DIR/preproc/clean_samples.tsv | awk "NR==$SLURM_ARRAY_TASK_ID")
export sample=$(echo -e "$__sample_line" | cut -f1)
export signature="$PARENT_DIR/genome_sketches/moss_samples/${sample}.sig"
export __fastq_dir="$PARENT_DIR/preproc/$sample"

if [[ ! -f $signature ]]; then
$sourmash sketch dna -p k=31,scaled=1000,abund --merge ${sample} -o $signature \
	${__fastq_dir}/${sample}_paired_1.fastq ${__fastq_dir}/${sample}_paired_2.fastq \
	${__fastq_dir}/${sample}_unmatched_1.fastq ${__fastq_dir}/${sample}_unmatched_2.fastq
fi

if [[ ! -f $PARENT_DIR/SM_abund/${sample}_custom_gather.csv ]]; then
$sourmash gather $signature \
	$PARENT_DIR/genome_sketches/index_custom_k31.sbt.zip \
	-o $PARENT_DIR/SM_abund/${sample}_custom_gather.csv
fi

if [[ ! -f $PARENT_DIR/SM_abund/${sample}_genbank_default_gather.csv ]]; then
$sourmash gather $signature \
	/nfs3_ib/nfs-ip34/fast/def-ilafores/sourmash_db/genbank-2022.03*k31.zip \
	-o $PARENT_DIR/SM_abund/${sample}_genbank_default_gather.csv
fi

if [[ ! -f $PARENT_DIR/SM_abund/${sample}_gtdb_gather.csv ]]; then
$sourmash gather $signature \
	/nfs3_ib/nfs-ip34/fast/def-ilafores/sourmash_db/gtdb-rs214-reps.k31.zip \
	-o $PARENT_DIR/SM_abund/${sample}_gtdb_gather.csv
fi

if [[ ! -f $PARENT_DIR/SM_abund/${sample}_hg38_gather.csv ]]; then
$sourmash gather $signature \
        /nfs3_ib/nfs-ip34/fast/def-ilafores/sourmash_db/gtdb-rs214-reps.k31.zip \
        -o $PARENT_DIR/SM_abund/${sample}_gtdb_gather.csv
fi
