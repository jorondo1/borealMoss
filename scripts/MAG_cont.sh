# fetch and index reference genomes 
module load bowtie2/2.5.1

# Index reference genome
for MAG in $(find MAG_analysis/novel_species/genomes -type f -name '*.fa'); do
	MAG_NAME=$(basename ${MAG%.fa})
	OUT_DIR=$(dirname $(realpath ${MAG}))/${MAG_NAME}_idx
	mkdir -p $OUT_DIR
	bowtie2-build --threads 48 ${MAG} $OUT_DIR/${MAG_NAME}
done

# Align sample reads
mkdir -p MAG_analysis/contamination_check