grep root $PWD/coassembly/binning/*/*.stats | grep -v "unbinned|binner" | \
	sed 's/.stats:/\//g' | awk -F'\t' '{print $1 "\t" $8}' > EukCC/all_root_bins.txt


# # catch bins that are not in the final set:
# for bin in $(ls MAG_analysis/all_bins/*.fa); do
# 	name=$(basename $bin)
# 	if [[ ! -f "MAG_analysis/drep_genomes/dereplicated_genomes/${name}" ]]; then
# 	cp $bin EukCC/candidate_bins/$name
# 	fi
# done

#EukRep -i coassembly/bin_refinement/Brown_AE/concoct_bins/bin.63.fa -o Euk.Brown_AE.bin.63.fa
export EUKCC="singularity exec  -e -B $ANCHOR:$ANCHOR $ILAFORES/programs/eukcc_latest.sif eukcc"
ml apptainer
$EUKCC folder EukCC/candidate_bins --out EukCC/out --threads 48 --db "$EUKCC2_DB"

$EUKCC single GCA_036851105.1_ASM3685110v1_genomic.fna --out EukCC/single_test/ \
	--DNA --threads 48 --db $EUKCC2_DB


