MOSS=/home/def-ilafores/analysis/boreal_moss
assemblies=$MOSS/coassembly/assembly
bins=${MOSS}/MAG_analysis/all_bins
cd $MOSS
module load java

####################################
### Submitting raw (clean) samples #
####################################

mkdir -p ENA_submission/run_manifests

while IFS=$'\t' read -r id sample; do 
file="ENA_submission/run_manifests/${id}_manifest.txt"
> ${file}
echo -e "STUDY\tPRJEB76464" >> "$file" 
echo -e "SAMPLE\t${id}" >> "$file"
echo -e "NAME\tBoreal Moss Microbiome ${sample}" >> "$file"
echo -e "INSTRUMENT\tIllumina NovaSeq 6000" >> "$file"
echo -e "INSERT_SIZE\t150" >> "$file"
echo -e "LIBRARY_SOURCE\tMETAGENOMIC" >> "$file"
echo -e "LIBRARY_SELECTION\tRANDOM" >> "$file"
echo -e "LIBRARY_STRATEGY\tWGS" >> "$file"
echo -e "FASTQ\t$MOSS"/preproc/$sample/${sample}_1.fastq.gz >> "$file"
echo -e "FASTQ\t$MOSS"/preproc/$sample/${sample}_2.fastq.gz >> "$file"

done < <(tail -n +2 ENA_submission/samples_report.tsv | awk -F'\t' '{print $1 "\t" $2}')

for man in $(find ENA_submission/run_manifests -type f -name '*manifest.txt'); do
java -jar /home/def-ilafores/programs/webin-cli-7.2.1.jar \
	-context reads -userName Webin-67053 -password LRCGsnth==1 \
	-submit -manifest "$man"
done

###################################
### Submitting primary assemblies #
###################################

# These are co-assemblies, therefore I am supplying the original samples they are from in the REF_RUN field, which can be removed for regular assemblies.

# Requires a reference file mapping which samples were concatenated
# here it's ./unique_combinations.txt

# Also requires the accessions-to-sample name file

for dir in $(ls "$assemblies"); do
	gzip -c $assemblies/$dir/final_assembly.fasta > $assemblies/$dir/${dir}.fa.gz
done

mkdir -p ENA_submission/assembly_manifests

while IFS=$'\t' read -r id sample; do 

if [[ ! -d "$assemblies/$sample" ]]; then
	continue #skip iteration if assembly doesn't exist
fi

file="ENA_submission/assembly_manifests/${id}_manifest.txt"
COM=${sample%_*} # compartment
MS=${sample#*_} # microsite
# find samples from this compartment & microsite
sam_list=($(cut -f1,5,6 comptype.tsv | awk -v com="$COM" -v ms="$MS" '$2 == com && $3 == ms' | cut -f1))
echo ${sam_ID[@]}

# extract their corresponding ENA id
sam_ID=()
for samples in "${sam_list[@]}"; do
	sam_ID+=($(grep -l $samples ENA_submission/run_manifests/*.txt | xargs grep SAMPLE | cut -f2))
done
# compute coverage 

> ${file}
echo -e "STUDY\tPRJEB76464" >> "$file"
echo -e "SAMPLE\t${id}" >> "$file" ### MULTIPLE ?!
echo -e "ASSEMBLYNAME\tBoreal Moss Metagenome Assembly" >> "$file" 
echo -e "ASSEMBLY_TYPE\tprimary metagenome" >> "$file" 
echo -e "COVERAGE\t$(cat $assemblies/$sample/coverage/average_coverage.txt)" >> "$file" 
echo -e "PROGRAM\tSPAdes3.15.4,MEGAHIT1.2.9" >> "$file" 
echo -e "PLATFORM\tIllumina Novaseq 6000" >> "$file" 
echo -e "DESCRIPTION\tCoassembly of multiple samples by host species, gametophyte section and microsite." >> "$file" 
echo -e "FASTA\t$assemblies/$sample/$sample.fa.gz" >> "$file"
echo -e "RUN_REF\t$(echo ${sam_ID[@]} | tr ' ', ',')" >> "$file"

done < <(tail -n +2 ENA_submission/accessions_coass.tsv | awk -F'\t' '{print $2 "\t" $3}')

for man in $(find ENA_submission/assembly_manifests -type f -name '*manifest.txt'); do
java -jar /home/def-ilafores/programs/webin-cli-7.2.1.jar \
	-context genome -userName Webin-67053 -password LRCGsnth==1 \
	-submit -test -manifest "$man" #THIS IS A TEST 
done

##################################
### Submitting binned assemblies #
##################################

for bin in $(find ${bins} -name '*.fa'); do
	gzip -k $bin
done 

mkdir -p ENA_submission/bins_manifests



