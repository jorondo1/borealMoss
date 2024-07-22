#!/bin/bash -l

#SBATCH --mail-type=END,FAIL
#SBATCH -D /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/boreal_moss
#SBATCH -o /nfs3_ib/nfs-ip34/home/def-ilafores/analysis/boreal_moss/EukCC/logs/eukCC-%A_%a.slurm.out
#SBATCH --time=4:00:00
#SBATCH --mem=30G
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -A def-ilafores
#SBATCH -J EukCC

if [ $# -lt 3 ]; then
  echo "Error: Missing arguments. Please provide three positional arguments."
  echo "Usage: script.sh BIN_LIST BINNING_DIRECTORY/ PARENT_OUTPUT_DIRECTORY"
  exit 1
fi

export BIN_LIST="${1}"
export BIN_FILE=$(awk "NR==$SLURM_ARRAY_TASK_ID" ${BIN_LIST} | cut -f1).fa
export BINNER=$(awk "NR==$SLURM_ARRAY_TASK_ID" ${BIN_LIST} | cut -f2)
export BINNING_DIR="${2}"
export SAMPLE=$(echo $BIN_FILE | sed "s|${BINNING_DIR}||" | sed 's,/.*,,')
export OUT_DIR=${3}/${SAMPLE}/${BINNER}
export EUKCC="singularity exec -e -B $ILAFORES:$ILAFORES -B $ANCHOR/fast:$ANCHOR/fast $ILAFORES/programs/eukcc.sif eukcc"

mkdir -p ${OUT_DIR}

# Run EukCC2
ml singularity
$EUKCC single ${BIN_FILE} --out ${OUT_DIR} --DNA --threads $SLURM_NTASKS --db $EUKCC2_DB

dos2unix ${OUT_DIR}/eukcc.log # output has a windows encoding, wtf

# add clade to output (using logs)
CLADE=$(grep -h 'Genome belongs to clade:' ${OUT_DIR}/eukcc.log |  awk -F 'clade: ' '{print $2}' | awk '{print $1}')

if [[ ! -f ${3}/full_results.csv ]]; then
	echo -e "$(head -n 1 ${OUT_DIR}/eukcc.log)\tClade" > ${3}/full_results.csv
fi

# Append results
echo -e "$(tail -n 1 ${OUT_DIR}/eukcc.csv)\t${CLADE}" >> ${3}/full_results.csv

echo "Done !"