# Index reference genome
for MAG in $(find MAG_analysis/novel_species/genomes -type f -name '*.fa'); do
	MAG_NAME=$(basename ${MAG%.fa})
	OUT_DIR=$(dirname $(realpath ${MAG}))/${MAG_NAME}_idx
	mkdir -p $OUT_DIR
	bowtie2-build --threads 48 ${MAG} $OUT_DIR/${MAG_NAME}
done

# Align sample reads
mkdir -p MAG_analysis/contamination_check/logs
export NOVEL="$MOSS"/MAG_analysis/novel_species
export MAG_LIST="$NOVEL"/novel_MAGs.txt 
export SAMPLE_LIST="$MOSS"/cat_reads/cat_samples.tsv
export NUM_MAGS=$(wc $MAG_LIST | awk '{print $1}')
export REF_GENOMES="$ANCHOR"/fast/def-ilafores/host_genomes
export POLCOM_IDX="$REF_GENOMES"/PolCom_index/PolCom_genome
export SOURMASH="singularity exec -e -B $ANCHOR/home:$ANCHOR/home -B $ANCHOR/fast:$ANCHOR/fast $ANCHOR/home/def-ilafores/programs/ILL_pipelines/containers/sourmash.4.7.0.sif sourmash"

# if [[ ! -f $POLCOM_SIG ]]; then
# 	$SOURMASH sketch dna -p scaled=1000,k=31 --name-from-first \
# 	$REF_GENOMES/_genome_fasta/GCA_950295325.1_cbPolComm4.1_genomic.fna -o $POLCOM_SIG
# fi
sbatch --array=1-"${NUM_MAGS}" -n 6 --mem=2G \
	$PWD/scripts/MAG_cont_SLURM.sh $NOVEL/genomes $MAG_LIST $SAMPLE_LIST $POLCOM_IDX
