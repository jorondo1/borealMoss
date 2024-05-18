newgrp def-ilafores
ILAFORES=/home/def-ilafores
PARENT_DIR=${ILAFORES}/analysis/boreal_moss
ILL_PIPELINES=${ILAFORES}/programs/ILL_pipelines
SOURMASH="/home/def-ilafores/programs/ILL_pipelines/containers/sourmash.4.7.0.sif sourmash"
MAG_DIR="$PARENT_DIR/MAG_analysis"
cd $PARENT_DIR
source /home/ronj2303/functions.sh 
	
###############
### CLEANING ###
###############
while read -r SAM A B; do
	SAM1=$(find $(realpath ../20220825_boreal_moss/data/) -type f -name "*${SAM}_R1.fastq.gz")
	SAM2=$(find $(realpath ../20220825_boreal_moss/data/) -type f -name "*${SAM}_R2.fastq.gz")
	echo -e "${SAM}\t${SAM1}\t${SAM2}" >> raw_samples.tsv
done < metadata.txt

mkdir -p $PARENT_DIR/preproc
export HOST=/cvmfs/datahub.genap.ca/vhost34/def-ilafores/host_genomes
bash $ILL_PIPELINES/generateslurm_preprocess.kneaddata.sh \
--sample_tsv ${PARENT_DIR}/raw_samples.tsv \
--out ${PARENT_DIR}/preproc \
--db "${HOST}/Pschreberi_index/Pschreberi ${HOST}/Ppatens_index/Ppatens" \
--slurm_email "ronj2303@usherbrooke.ca" \
--trimmomatic_options "SLIDINGWINDOW:4:30 MINLEN:50" \
--slurm_walltime "12:00:00"

find preproc/S-* -name '*paired_1*.fastq' | parallel -j 72 'echo -n {}" "; grep -c "^@" {}' > sequence_counts.txt
find ../20220825_boreal_moss/data -name '*.S-*R1.fastq.gz' | parallel -j 72 "echo -n {} ' '; zcat {} | grep -c '^@'" > raw_sequence_counts.txt


#######################
### REGULAR ASSEMBLY ###
#######################
:> preproc/clean_samples.tsv
while read SAM A1 A2; do
	SAM1=$(find preproc/${SAM} -type f -name '*paired_1.fastq' -exec realpath {} \;)
	SAM2=$(find preproc/${SAM} -type f -name '*paired_2.fastq' -exec realpath {} \;)
	echo -e "${SAM}\t${SAM1}\t${SAM2}" >> preproc/clean_samples.tsv
done < metadata.txt # NEED DIFERENT FILE FOR CAT

bash $ILL_PIPELINES/generateslurm_assembly_bin_refinement.metawrap.sh \
--sample_tsv ${PARENT_DIR}/preproc/clean_samples.tsv \
--out ${PARENT_DIR}/assembly \
--slurm_email "ronj2303@usherbrooke.ca" \
--refinement_min_compl 50 --refinement_max_cont 10 \
--slurm_walltime 72:00:00

###################
### CO ASSEMBLY ###
##################
# sbatch --array=1-50 /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/boreal_moss/coassembly/assembly_bin_refinement.metawrap.slurm.sh
cat_samples() {
    local metadata_file="$1"

    # Extract unique combinations of COM and MS from metadata
    cut -f2,3 "$metadata_file" | sort -u > unique_combinations.txt

    # Loop through the unique combinations and concatenate the files
    while read -r COM MS; do
        # Create an array to store the filenames matching the combination
        file_array_1=()
        file_array_2=()
        # Loop through SAM values to find files matching the combination
        while read -r SAM; do
            file_array_1+=( "preproc/$SAM/${SAM}_paired_1.fastq" )
            file_array_2+=( "preproc/$SAM/${SAM}_paired_2.fastq" )
        done < <(grep -P "$COM\t$MS" "$metadata_file" | cut -f1)

        # Concatenate the matching files to create the output file
        cat "${file_array_1[@]}" > cat_reads/"${COM}_${MS}_1.fastq"
        cat "${file_array_2[@]}" > cat_reads/"${COM}_${MS}_2.fastq"
    done < unique_combinations.txt
}
mkdir -p cat_reads
cat_samples metadata.txt

sed 's/\t/_/g' unique_combinations.txt > cat_samples.txt
tsv=${PARENT_DIR}/cat_reads/cat_samples.tsv
:> $tsv
while read SAM; do
	SAM1=$(find cat_reads -type f -name "${SAM}_1.fastq" -exec realpath {} \;)
	SAM2=$(find cat_reads -type f -name "${SAM}_2.fastq" -exec realpath {} \;)
	echo -e "${SAM}\t${SAM1}\t${SAM2}" >> $tsv
done < cat_samples.txt

# NODES:
bash $ILL_PIPELINES/generateslurm_assembly_bin_refinement.metawrap.sh \
--sample_tsv $tsv \
--out ${PARENT_DIR}/coassembly \
--slurm_email "ronj2303@usherbrooke.ca" \
--refinement_min_compl 50 --refinement_max_cont 10 \
--slurm_walltime 72:00:00 --slurm_threads 48 --slurm_mem 240G

### OR LOOP THROUgH locally on ip34
mkdir -p coassembly/logs
while read sample; do
	log_file="coassembly/logs/${sample}_out.log"
	bash /home/def-ilafores/programs/ILL_pipelines/scripts/fastq_to_bins.metawrap.sh \
	-o /home/def-ilafores/analysis/boreal_moss/coassembly \
	-tmp /fast/def-ilafores/temp/${sample} \
	-t 48 -m 250G \
	-s ${sample} \
	-fq1 /home/def-ilafores/analysis/boreal_moss/cat_reads/${sample}_1.fastq \
	-fq2 /home/def-ilafores/analysis/boreal_moss/cat_reads/${sample}_2.fastq \
	--metaspades --megahit \
	--metabat2 --maxbin2 --concoct --run-checkm \
	--refinement_min_compl 50 --refinement_max_cont 10 > "$log_file" 2>&1
done < cat_left.txt

#mkdir -p MAG_analysis/all_bins && cp coassembly/bin_refinement/*/metawrap_50_10_bins/*.fa $_

####################
### DEREPLICATION ###
###################

# Gather/rename all bins
DREP_OUT=${PARENT_DIR}/MAG_analysis/drep_genomes
mkdir -p $DREP_OUT ${PARENT_DIR}/MAG_analysis/all_bins

for bin in $(find ${PARENT_DIR}/coassembly/bin_refinement/*/metawrap_50_10_bins/ -type f); do
	SAM=$(echo $bin | sed 's/.*bin_refinement\///' | sed 's/\/metawrap.*//')
	BIN=$(basename $bin)
	cp $bin MAG_analysis/all_bins/${BIN/#bin/${SAM}_bin}
done

bash $ILL_PIPELINES/scripts/dereplicate_bins.dRep.sh \
	-o $DREP_OUT \
	-t 32 \
	-bin_path_regex "MAG_analysis/all_bins/*.fa" \
	-comp 50 \
	-con 10

#####################
### CLEAN UP GUNC ###
#####################

module load apptainer
mkdir -p $DREP_OUT/gunc
singularity exec -e -B /home:/home -B /fast:/fast \
	/home/def-ilafores/programs/gunc_1.0.5.sif \
	gunc run -t 12 -d $DREP_OUT/dereplicated_genomes \
	--db_file /fast/def-ilafores/gunc_db/gunc_db_progenomes2.1.dmnd \
	-o $DREP_OUT/gunc

# Remove contaminated genomes
cat $DREP_OUT/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv | \
	grep -v "n_genes_mapped" | \
	awk '{if($8 > 0.45 && $9 > 0.05 && $12 > 0.5)print$1}' > $DREP_OUT/gunc/gunc_contaminated.txt

for R in $(cat $DREP_OUT/gunc/gunc_contaminated.txt); do 
rm $DREP_OUT/dereplicated_genomes/${R}.fa 
done

# Print assembly stats (custom function in /home/ronj2303/functions.sh)
assembly_stats

##############
# Annotation
##############
# mash_db=/cvmfs/datahub.genap.ca/vhost34/def-ilafores/GTDB/release207_v2
# mkdir -p coassembly/gtdb_tmp; tmp=$(realpath gtdb_tmp)
# singularity exec --writable-tmpfs -e \
# --env GTDBTK_DATA_PATH=$gtdb_db \
# -B $tmp:$tmp \
# -B $gtdb_db:$gtdb_db \
# -e $ILL_PIPELINES/containers/gtdbtk.2.2.5.sif \
# 	gtdbtk classify_wf --cpus 32 --genome_dir  $DREP_OUT/dereplicated_genomes\
# 	--out_dir MAG_analysis/Annotation --mash_db $gtdb_db --extension fa
#
# The whole annotation pipeline:
bash $ILL_PIPELINES/scripts/annotate_bins.sh \
	-o MAG_analysis/Annotation \
	-t 72 -drep $DREP_OUT/dereplicated_genomes

###################################
# Novelty : mash distance + ANI ####
#####################################
mkdir -p $MAG_DIR/mash_dist

# Build sketch
cd $ILAFORES/EnsemblBacteria57
cp $PARENT_DIR/Hannah_MAGs/*.gz .
find -type f -name '*.fa.gz' | sed 's/\.\///' > $MAG_DIR/mash_dist/genomes_list.txt
cd $MAG_DIR/mash_dist

module load StdEnv/2020  gcc/9.3.0 mash/2.3
novel_genomes #from custom functions in /home/ronj2303/functions.sh
cd ../..

## Adding a layer : skani validation of genome demarcation
# Setup : https://github.com/bluenote-1577/skani/wiki/Tutorial:-setting-up-the-GTDB-genome-database-to-search-against
skani=/home/def-ilafores/programs/skani/skani
cd /fast/def-ilafores/GTDB
find gtdb_genomes_reps_r214/ -name '*.fna.gz' > gtdb_file_names.txt
"$skani" sketch -l gtdb_file_names.txt -o gtdb_skani_database_ani -t 72

# Compute ANI
cd $MAG_DIR/novel_species
find $PWD/genomes/ -type f -name '*.fa' > novel_species_mash.txt
"$skani" search -d /fast/def-ilafores/GTDB/gtdb_skani_database_ani -o ANI_results.txt -t 72 --ql novel_species_mash.txt 

# Compute additionnal MAG statistics
module load bbmap/38.86
statswrapper.sh MAG_analysis/novel_species/genomes/*.fa > MAG_analysis/novel_species/genomes/genome.stats

#############
## SOURMASH #
#############

# list genomes to sketch (full path)
mkdir -p genome_sketches && cd $_
mkdir -p moss_MAGs ensembl_bact

:> moss_MAGs.txt
find ../Hannah_MAGs -type f -name '*.f*a*' -exec realpath {} \; >> moss_MAGs.txt
find $MAG_DIR/all_bins/ -type f -name '*.f*a*' -print0 | xargs -0 realpath | 
	grep -wf <(cat $MAG_DIR/novel_species/mash_dist/novel_MAGs.txt | cut -d ',' -f1) >> moss_MAGs.txt

# sketch novel genomes
ml apptainer; tmp=$PWD
singularity exec --writable-tmpfs -e -B /home:/home -B $tmp:$tmp \ 
	$SOURMASH sketch dna -p scaled=1000,k=31,abund \
	--name-from-first --from-file moss_MAGs.txt \
		--output-dir moss_MAGs

ml apptainer; tmp=$PWD
singularity exec --writable-tmpfs -e -B /home:/home -B $tmp:$tmp -B /fast:/fast \
	$SOURMASH index /fast/def-ilafores/sourmash_db/genbank-2022.03-index-k31 \
	/fast/def-ilafores/sourmash_db/genbank-2022.03*k31.zip

#check signature names
# singularity exec --writable-tmpfs -e -B /home:/home -B $tmp:$tmp \
# 	$SOURMASH sketch dna -p scaled=1000,k=31,abund \
# 	MAG_analysis/all_bins/Brown_AD.bin.9.fa --output-dir .
#
# singularity exec --writable-tmpfs -e -B /home:/home -B $tmp:$tmp \
# 	$SOURMASH sig describe ./Brown_AD.bin.9.fa.sig
#

mkdir moss_MAGs_renamed/
for file in $(find moss_MAGs -type f -name '*.sig'); do
	new_name=$(basename $file)
	singularity exec --writable-tmpfs -e -B /home:/home -B $tmp:$tmp \
	$SOURMASH sig rename $file "${new_name%.fa}" -o moss_MAGs_renamed/${new_name}
done


#################################
### CONTAINMENT & ABUNDANCE #####
#################################
# cd genome_sketches
singularity exec --writable-tmpfs -e -B /home:/home -B $tmp:$tmp -B /fast:/fast \
	$SOURMASH index index_custom_k31 /fast/def-ilafores/sourmash_db/gtdb-rs214-reps.k31.zip moss_MAGs_renamed/*

# Comparing with Genbank
singularity exec --writable-tmpfs -e -B /home:/home -B $tmp:$tmp -B /fast:/fast \
	$SOURMASH index index_genbank_k31 /fast/def-ilafores/sourmash_db/genbank-2022.03*k31.zip moss_MAGs_renamed/*
cd ..

# Gather on both databases
export gather=$PWD/SM_abund
mkdir -p genome_sketches/moss_samples $gather/logs
sbatch --array=1-140 $PWD/scripts/sourmash_compare_custom.sh

# If some didn,t work, find the missing array numbers you need to rerun :
redo=$(grep -vnf <(ls SM_abund/*genbank_default* | sed 's/.genbank_default_gather.csv//' | sed 's/SM_abund\///') preproc/clean_samples.tsv | cut -d':' -f1 |  paste -sd,)
sbatch --array="$redo" $PWD/scripts/sourmash_compare_custom.sh

# HOST CONTAMINATION
for file in $(find genome_sketches/moss_samples -type f -name '*sig' | grep POLCOM); do
sample=$(basename $file)
singularity exec --writable-tmpfs -e -B /home:/home -B /fast:/fast -B $tmp:$tmp \
	$SOURMASH gather $file /fast/def-ilafores/host_genomes/sourmash_sketches/GCA_950295325.1_cbPolComm4.1_genomic.fna.sig \
	-o SM_abund/${sample%.sig}_POLCOM_gather.csv
done

fix_gtdb # GTDB taxonomy has commas in it, this messes up the column recognition we need below. Fixing it:
eval_cont SM_abund # Compare containment before and after adding MAGs

# Assign taxonomy
cat /fast/def-ilafores/sourmash_db/gtdb-rs214.lineages.csv \
	MAG_analysis/mash_dist/novel_MAGs.txt > genome_sketches/custom_lineages.csv
# $sourmash tax prepare --taxonomy $SM_DB/gtdb-rs214.taxonomy.sqldb $SM_DB/gtdb-rs214.lineages.csv.gz -o tax.db

# Subset gtdb taxonomy with species found only
head -n1 /fast/def-ilafores/sourmash_db/gtdb-rs214.lineages.csv > SM_abund/gtdb_taxonomy_subset.csv
cat SM_abund/*.csv | cut -d',' -f10 | sed 's/\ .*//' | sed 's/\"//' | sort -u | \
	grep -wf - /fast/def-ilafores/sourmash_db/gtdb-rs214.lineages.csv >> SM_abund/gtdb_taxonomy_subset.csv


#######################################
#### FUNCTIONAL ANNOTATION #############
#######################################

# Run HUMANN on clean reads to get community functional potential profile
# bash /nfs3_ib/nfs-ip34$ILL_PIPELINES/generateslurm_functionnal_profile.humann.sh \
# 	--sample_tsv /nfs3_ib/nfs-ip34$PWD/preproc/preprocessed_reads.sample.tsv \
# 	--slurm_threads 24 --slurm_mem 30G --slurm_walltime 72:00:00 \
# 	--out /nfs3_ib/nfs-ip34$PWD/HUMANN

#### We need all genomes found (not just the MAGs) to be annotated
### We download them, then we'll annotate them all together
readarray -t ids < <(cat SM_abund/*custom_gather.csv | awk -F\" '{print $2}' | awk -F' ' '{print $1}' | sort -u)
mkdir -p Annotations/translated_genomes && cd $_

# Download CDS of all genomes found by Sourmash in GTDB (via Genbank)
for id in "${ids[@]}"; do
	curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${id}/download?include_annotation_type=PROT_FASTA&filename=${id}.zip" -H "Accept: application/zip"
    unzip -o ${id}.zip
	mv ncbi_dataset/data/${id}/*protein.faa ./${id}.faa
	rm ${id}.zip
done
rm -r README.md ncbi_dataset
cd ../..

# Add our novel genomes to these
cp Annotations/metawrap_out/bin_translated_genes/*faa Annotations/translated_genomes/

find $PWD/Annotations/translated_genomes -type f -name '*.faa' -exec realpath {} \; > Annotations/translated_genomes_list.txt

###### MISSING ANNOTATIONS
# Not all genomes will have coding sequences in ncbi. Use prodigal to predict them
# List missing:
readarray -t ids_missing < <(cat SM_abund/*custom_gather.csv | awk -F\" '{print $2}' | awk -F' ' '{print $1}' | sort -u | \
    grep -vf <(head -n 1 Annotations/microbeannotator_out/metabolic_summary__module_completeness.tab | sed 's/\.faa\.ko//g' | tr '\t' '\n'))
mkdir -p Annotations/missing_cds && cd $_

for id in "${ids_missing[@]}"; do
	curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${id}/download?include_annotation_type=GENOME_FASTA&filename=${id}.zip" -H "Accept: application/zip"
    unzip -o ${id}.zip
	mv ncbi_dataset/data/${id}/*genomic.fna ./${id}.fna
	rm ${id}.zip
done
rm -r README.md ncbi_dataset
cd ../..

# Predict coding sequence using Prodigal
module load prodigal/2.6.3
for fna in $(find . -type f -name '*.fna' -print0 | xargs -0 -I {} basename {}); do
	prodigal -q -i ${fna} -a ${fna%.fna}.faa
done

find $PWD/Annotations/missing_cds -type f -name '*.faa' -exec realpath {} \; >> Annotations/translated_genomes_list.txt
# find $PWD/Annotations/missing_cds -type f -name '*.faa' -exec realpath {} \; > Annotations/translated_genomes_list_missing.txt

# Run MicrobeAnnotator
## Make sure /fast/def-ilafores/MicrobeAnnotator_DB is there (otherwise it's in /home/def-ilafores/ref_dbs/MicrobeAnnotator_DB)
singularity exec --writable-tmpfs -e \
--env MPLCONFIGDIR=$tmp \
-B /home:/home -B /cvmfs:/cvmfs \
$ILL_PIPELINES/containers/microbeannotator.2.0.5.sif \
microbeannotator --method diamond --processes 12 --threads 4 --refine \
-l $PWD/Annotations/translated_genomes_list_missing.txt \
-d /home/def-ilafores/ref_dbs/MicrobeAnnotator_DB \
-o $PWD/Annotations/microbeannotator_out_missing

#Some didn't work, here's how I listed them:
grep -vf <(find Annotations/microbeannotator_out/annotation_results/ -type f -name '*.ko' -print0 | xargs -0 -I {} basename {} | sed 's/\.faa\.ko//' | grep -v ".bin.") <(find Annotations/found_genomes/ -type f -exec realpath {} \;) | grep -v "fa" > Annotations/try_again.tsv
head Annotations/try_again.tsv

#### Finding fungi??
singularity exec --writable-tmpfs -e -B /home:/home -B $tmp:$tmp -B /fast:/fast \
	$SOURMASH gather genome_sketches/moss_samples/S-18-POLJUN-G.sig /fast/def-ilafores/sourmash_db/genbank-2022.03-fungi-k31.zip

##############################
### FULL PHYLOGENY ################
##############################

# Download relevant genomes
readarray -t ids < <(cut -d',' -f1 SM_abund/gtdb_taxonomy_subset.csv | grep GC)
mkdir -p Phylogeny_fast/genomes && cd $_

# Download all genomes found by Sourmash in GTDB (via Genbank)
for id in "${ids[@]}"; do
	curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${id}/download?include_annotation_type=GENOME_FASTA&filename=${id}.zip" -H "Accept: application/zip"
    unzip -o ${id}.zip
	mv ncbi_dataset/data/${id}/*genomic.fna ./${id}.fna
	rm ${id}.zip
done
rm -r README.md ncbi_dataset
cp /home/def-ilafores/analysis/boreal_moss/MAG_analysis/novel_species/genomes/* .
for file in *.fa; do
    mv "$file" "$(echo $file | sed 's/\.fa$/.fna/')"
done
cd ../..


# Run Phylophlan
# Directories might be wrong, needs REVISION !!!!!
ml apptainer

singularity exec --writable-tmpfs -e -B /fast:/fast \
-B /home:/home \
/home/def-ilafores/programs/ILL_pipelines/containers/phylophlan.3.0.3.sif \
phylophlan --nproc 48 \
--verbose --genome_extension .fna \
-t a -i /home/def-ilafores/analysis/boreal_moss/Annotations/found_genomes \
-o /fast/def-ilafores/temp/phylophlan_acc_full \
--diversity high --accurate \
--databases_folder /fast/def-ilafores/phylophlan_db -d phylophlan \
-f /home/def-ilafores/analysis/boreal_moss/Phylogeny_accurate/phylophlan_nodes.cfg 
# To rerun if crashed: find Phylogeny/Phylophlan -type f -name "*.bkp" -size 0 -exec rm {} \;

# We -output to /fast/ to avoid some performance problems, lets copy results
mkdir -p Phylogeny_acc_full
cp /fast/def-ilafores/temp/phylophlan_acc_full/*.* Phylogeny_acc_full

###################################
## What the fuck is in our reads ##
###################################

mkdir -p DarkMatter/fasta DarkMatter/blast/nt
module load gcc/9.3.0 blast+/2.12.0 seqtk/1.3
samples=$(find ./preproc -maxdepth 1 -type d -name 'S-*' -exec basename {} \;)

# Random select subset of reads
for sample in $samples; do
seqtk sample -s100 preproc/${sample}/${sample}_paired_1.fastq 100000 | \
awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' > DarkMatter/fasta/${sample}_subset.fa
done

# Blasting against ncbi 
for sample in $samples; do
blastn -query DarkMatter/fasta/${sample}_subset.fa \
	-db /cvmfs/bio.data.computecanada.ca/content/databases/Core/blast_dbs/2022_03_23/nt \
	-evalue 0.01 -qcov_hsp_perc 75 -word_size 20 -max_target_seqs 5 -num_threads 48 \
	-out DarkMatter/blast/nt/${sample}_subset.blastout \
	-outfmt "6 staxids slen qstart qend sstart send evalue bitscore score length pident nident mismatch qseqid qlen sseqid"
done

# For now we don't need much info
mkdir -p DarkMatter/blast/nt/staxid_slen
for file in $(find DarkMatter/blast/nt -name '*blastout'); do
	name=$(basename $file)
	cut -f1,4 $file > DarkMatter/blast/nt/staxid_slen/${name}
done

tar -zcvf blast.tar.gz $(find DarkMatter/blast/nt/staxid_slen -name "*.blastout")
