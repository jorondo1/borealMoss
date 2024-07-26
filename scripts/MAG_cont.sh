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


## code snippet at the end of moss_contamination.R used to explore results

### Use MUMMER to compare the suspicious Buchnera MAG to a Buchnera reference genome
mkdir -p Buchnera && cd $_

# Fetch genomes
GET_NCBI="wget -qO- https://ftp.ncbi.nlm.nih.gov/genomes/all"
$GET_NCBI/GCA/016/903/675/GCA_016903675.1_ASM1690367v1/GCA_016903675.1_ASM1690367v1_genomic.fna.gz | gzip -d > StiOce_genomic.fna
$GET_NCBI/GCF/003/099/975/GCF_003099975.1_ASM309997v1/GCF_003099975.1_ASM309997v1_genomic.fna.gz | gzip -d > BucAph_genomic.fna
$GET_NCBI/GCA/950/295/325/GCA_950295325.1_cbPolComm4.1/GCA_950295325.1_cbPolComm4.1_genomic.fna.gz | gzip -d > PolComm4.1_genomic.fna
$GET_NCBI/GCF/016/027/095/GCF_016027095.1_ASM1602709v1/GCF_016027095.1_ASM1602709v1_genomic.fna.gz | gzip -d > SphPau_genomic.fna
$GET_NCBI/GCF/000/512/205/GCF_000512205.2_ASM51220v2/GCF_000512205.2_ASM51220v2_genomic.fna.gz | gzip -d > SphSan_genomic.fna
$GET_NCBI/GCF/010/450/875/GCF_010450875.1_ASM1045087v1/GCF_010450875.1_ASM1045087v1_genomic.fna.gz | gzip -d > SphIns_genomic.fna

# Align bacteria to P. commune
ml mummer

nucmer -p BucAph -t 24 --maxmatch PolComm4.1_genomic.fna BucAph_genomic.fna
nucmer -p StiOce -t 24 --maxmatch PolComm4.1_genomic.fna StiOce_genomic.fna
nucmer -p SphIns -t 24 --maxmatch PolComm4.1_genomic.fna SphIns_genomic.fna
nucmer -p SphPau -t 24 --maxmatch PolComm4.1_genomic.fna SphPau_genomic.fna
nucmer -p SphSan -t 24 --maxmatch PolComm4.1_genomic.fna SphSan_genomic.fna

dnadiff -p SphIns -d SphIns.delta
dnadiff -p SphPau -d SphPau.delta
dnadiff -p SphSan -d SphSan.delta 
find . -maxdepth 1 -type f ! -name '*.fna' -print0 | xargs -0 mv -t tmp/

export MAG_DIR=$MOSS/MAG_analysis/novel_species/genomes
for MAG in $(find $MAG_DIR -type f -name '*.fa'); do
	name=$(basename ${MAG%.fa})
	mkdir -p tmp && cd $_
	nucmer -p $name -t 48 --maxmatch ../PolComm4.1_genomic.fna $MAG
	dnadiff -p $name -d $name.delta
	cd ..
done

grep AlignedBases tmp/*.report > AlignedBases.txt

dnadiff -p StiOce -d StiOce.delta 
dnadiff -p BucMAG -d BucMAG.delta 
