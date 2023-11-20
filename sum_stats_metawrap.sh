#!/usr/bin/env bash -e

cd /home/def-ilafores/analysis/20220825_boreal_moss
for file in $(find denovo_assembly/bin_refinement -type f -name '*metawrap_50_10_bins.stats'); do
	echo $file
	SAM=$(echo $file | sed 's/.*bin_refinement\///' | sed 's/\/metawrap.*//')
	NUM=$(cat $file | tail -n +2 | wc -l)
	SIZE=$(find preprocess/"$SAM" -type f -name '*_paired_1.fastq' -exec wc -l {} + | awk '{total += $1} END {print total}')
	echo -e "${SAM}\t${NUM}\t${SIZE}" >> num_bins_per_sample_112.txt
done
