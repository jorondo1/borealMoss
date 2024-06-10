MOSS=/home/def-ilafores/analysis/boreal_moss
cd $MOSS

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
for fastq in $(find "$MOSS"/preproc/"$sample" -type f -name '*.fastq.gz'); do
	echo -e "FASTQ\t${fastq}" >> "$file"
done
done < <(tail -n +2 ENA_submission/samples_report.tsv | awk -F'\t' '{print $1 "\t" $2}')

for man in $(find ENA_submission/run_manifests -type f -name '*manifest.txt'); do
java -jar /home/def-ilafores/programs/webin-cli-7.2.1.jar \
	-context reads -userName Webin-67053 -password LRCGsnth==1 \
	-submit -manifest "$man"
done

###################################
### Submitting primary assemblies #
###################################

mkdir -p ENA_submission/assembly_manifests

assemblies=$MOSS/coassembly/assembly
for dir in $(ls "$assemblies"); do
	gzip -c $assemblies/$dir/final_assembly.fasta > $assemblies/$dir/${dir}_coassembly.fa.gz

while read -r sample Host; do 
file="ENA_submission/assembly_manifests/${sample}_manifest.txt"
> ${file}
echo -e "STUDY\tPRJEB76464" >> "$file" 
echo -e "SAMPLE\t${sample}" >> "$file" ### MULTIPLE ?!
echo -e "ASSEMBLYNAME\tBoreal Moss Metagenome Assembly"
echo -e "ASSEMBLY_TYPE\tprimary metagenome"
echo -e "COVERAGE\t"
echo -e "PROGRAM\tSPAdes3.15.4;MEGAHIT1.2.9"
echo -e "FASTQ\t$(find coassembly/assembly/$assembly -type f -name '*.fastq.gz')" >> "$file"
###### INCOMPLETE
done < <(tail -n +2 ENA_submission/samples_report.tsv)

for man in $(find ENA_submission/run_manifests -type f '*manifest.txt'); do
java -jar /home/def-ilafores/programs/webin-cli-7.2.1.jar \
	-context reads -userName Webin-67053 -password LRCGsnth==1 \
	-submit -manifest "$man"
done
