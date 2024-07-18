# catch bins that are not in the final set:
for bin in $(ls MAG_analysis/all_bins/*.fa); do 
	name=$(basename $bin)
	if [[ ! -f "MAG_analysis/drep_genomes/dereplicated_genomes/${name}" ]]; then 
	cp $bin EukCC/candidate_bins/$name
	fi
done

#EukRep -i coassembly/bin_refinement/Brown_AE/concoct_bins/bin.63.fa -o Euk.Brown_AE.bin.63.fa
ml apptainer
singularity exec  -e -B $ANCHOR:$ANCHOR /nfs3_ib/nfs-ip34/home/def-ilafores/ref_dbs/eukccdb/eukcc_latest.sif \
eukcc folder EukCC/candidate_bins --out EukCC/out --threads 48 --db "$EUKCC2_DB"
