MOSS=/home/def-ilafores/analysis/boreal_moss
cd $MOSS
module load java

## TODO 
## add missing coassemblies

####################################
### Submitting raw (clean) samples #
####################################

mkdir -p ENA_submission/run_manifests

while IFS=$'\t' read -r id sample; do 
file="ENA_submission/run_manifests/${id}_manifest.txt"
fq1=$(find $MOSS -mindepth 2 -maxdepth 3 -type f -name "${sample}_1.fastq.gz")
fq2=$(find $MOSS -mindepth 2 -maxdepth 3 -type f -name "${sample}_2.fastq.gz")
echo -e "STUDY\tPRJEB76464" > "$file"
echo -e "SAMPLE\t${id}" >> "$file"
echo -e "NAME\tBoreal Moss Microbiome ${sample}" >> "$file"
echo -e "INSTRUMENT\tIllumina NovaSeq 6000" >> "$file"
echo -e "INSERT_SIZE\t150" >> "$file"
echo -e "LIBRARY_SOURCE\tMETAGENOMIC" >> "$file"
echo -e "LIBRARY_SELECTION\tRANDOM" >> "$file"
echo -e "LIBRARY_STRATEGY\tWGS" >> "$file"
echo -e "FASTQ\t$fq1" >> "$file"
echo -e "FASTQ\t$fq2" >> "$file"

done < <(tail -n +2 ENA_submission/samples_report.tsv | awk -F'\t' '{print $1 "\t" $2}')

for man in $(find ENA_submission/run_manifests -type f -name '*manifest.txt'); do
java -jar /home/def-ilafores/programs/webin-cli-7.2.1.jar \
	-context reads -userName Webin-67053 -password LRCGsnth==1 \
	-submit -ascp -manifest "$man"
done

# RETRIEVE RUN REPORT and save as ENA_submission/run_report.tsv

###################################
### Submitting primary assemblies #
###################################

paste <(cut -f1,7,8 loctype.tsv) <(cut -f1,5,6 comptype.tsv) > unique_combinations.txt

# These are co-assemblies, therefore I am supplying the original samples they are from in the REF_RUN field, which can be removed for regular assemblies.
assemblies=$MOSS/coassembly/assembly
# Requires a reference file mapping which samples were concatenated
# here it's ./unique_combinations.txt

# Also requires the accessions-to-sample name file
mkdir -p ENA_submission/assembly_manifests

while IFS=$'\t' read -r id sample; do 

# sample reports contain all samples, we're only interested in coassembly names
if ! grep -qw "$sample" cat_reads/cat_samples.tsv; then # using a list of coassembly names
	continue # skip this iteration
fi

file="ENA_submission/assembly_manifests/${id}_manifest.txt"
COM=${sample%_*} # compartment
MS=${sample#*_} # microsite

# find samples from this compartment & microsite
sam_list=($(awk -v com="$COM" -v ms="$MS" '$2 == com && $3 == ms' unique_combinations.txt | cut -f1))

# extract their corresponding ENA id
sam_ID=()
for samples in "${sam_list[@]}"; do
	sam_ID+=($(grep -w $samples ENA_submission/run_report.tsv | cut -f1))
done
echo ${sam_ID[@]}

echo -e "STUDY\tPRJEB76464" > "$file"
echo -e "SAMPLE\t${id}" >> "$file" 
echo -e "ASSEMBLYNAME\tBoreal Moss Metagenome Co-Assembly ${COM}_${MS}" >> "$file" 
echo -e "ASSEMBLY_TYPE\tprimary metagenome" >> "$file" 
echo -e "COVERAGE\t$(cat $assemblies/$sample/coverage/average_coverage.txt)" >> "$file" 
echo -e "PROGRAM\tSPAdes3.15.4,MEGAHIT1.2.9" >> "$file" 
echo -e "PLATFORM\tIllumina Novaseq 6000" >> "$file" 
echo -e "DESCRIPTION\tCoassembly of multiple samples by host species, gametophyte section and microsite." >> "$file" 
echo -e "FASTA\t$assemblies/$sample/$sample.fa.gz" >> "$file"
echo -e "RUN_REF\t$(echo ${sam_ID[@]} | tr ' ', ',')" >> "$file"

done < <(tail -n +2 ENA_submission/samples_report.tsv | awk -F'\t' '{print $1 "\t" $2}')

for man in $(find ENA_submission/assembly_manifests -type f -name '*manifest.txt'); do
java -jar /home/def-ilafores/programs/webin-cli-7.2.1.jar \
	-context genome -userName Webin-67053 -password LRCGsnth==1 \
	-submit -manifest "$man" #THIS IS A TEST 
done


##################################
### Submitting MAGs [Incomplete] #
##################################
# UNTESTED #!!!
MAGs=${MOSS}/MAG_analysis/all_bins

mkdir -p ENA_submission/mags_manifests

while IFS=$'\t' read -r id mag; do 

file="ENA_submission/mags_manifests/${id}_manifest.txt"
COM=${mag%_*} # compartment
segment=${mag#*_}; MS=${mag%%.*} # microsite
NUM=${segment##*.} # bin number

# find samples from this compartment & microsite
sam_list=($(awk -v com="$COM" -v ms="$MS" '$2 == com && $3 == ms' unique_combinations.txt | cut -f1))
# extract their corresponding ENA id
sam_ID=()
for samples in "${sam_list[@]}"; do
	sam_ID+=($(grep -l $samples ENA_submission/run_manifests/*.txt | xargs grep SAMPLE | cut -f2))
done

echo -e "STUDY\tPRJEB76464" > ${file}
echo -e "SAMPLE\t${id}" >> "$file" 
echo -e "ASSEMBLYNAME\tBoreal Moss Metagenome Coassembled MAG ${COM}_${MS}_${NUM}" >> "$file" 
echo -e "ASSEMBLY_TYPE\tMetagenome-Assembled Genome (MAG)" >> "$file" 
echo -e "COVERAGE\t$(cat ${MAGs}/${mag}_coverage/average_coverage.txt)" >> "$file" 
echo -e "PROGRAM\tSPAdes3.15.4,MEGAHIT1.2.9" >> "$file" 
echo -e "PLATFORM\tIllumina Novaseq 6000" >> "$file" 
echo -e "DESCRIPTION\tCoassembly of multiple samples by host species, gametophyte section and microsite. Binned using Concoct, Maxbin2 and Metabat2, refined using metawrap's bin_refinement module." >> "$file" 
echo -e "FASTA\t$MAGs/${mag}.fa.gz" >> "$file"
echo -e "RUN_REF\t$(echo ${sam_ID[@]} | tr ' ', ',')" >> "$file"

# find samples from this compartment & microsite
done <(tail -n +2 ENA_submission/accessions_bins.tsv | awk -F'\t' '{print $2 "\t" $3}')

for man in $(find ENA_submission/mags_manifests -type f -name '*manifest.txt'); do
java -jar /home/def-ilafores/programs/webin-cli-7.2.1.jar \
	-context genome -userName Webin-67053 -password LRCGsnth==1 \
	-submit -manifest "$man" #THIS IS A TEST 
done
