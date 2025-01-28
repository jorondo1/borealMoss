## Bioinformatic scripts for _Boreal moss-microbe interactions are revealed through metagenome assembly of novel bacterial species_

Published in Nature Scientific Reports : [Ishak et al, 2024](https://www.nature.com/articles/s41598-024-73045-z)

The `moss_pipeline.sh` script contains the main steps for metagenome processing, assembly, and profiling. The pipeline in not seamless and must be run step by step. Each steps rely on modules of a pipeline implemented by our team [available here](https://github.com/jflucier/ILL_pipelines) and are scripted for SLURM workload manager. Several steps use singularity-contained programs to process data in parallel, but the underlying scripts ( [here](https://github.com/jflucier/ILL_pipelines/tree/main/scripts)) can be run on a sample-wise basis.

The `moss_pipeline.sh` script was run on Université de Sherbrooke's [Mammouth-mp2 clusters](https://docs.alliancecan.ca/wiki/Mp2/en), whereas R scripts were run locally on a Macbook Pro M2 through RStudio. Raw (clean) sequencing reads, assemblies, and MAGs were deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number [PRJEB76464](https://www.ebi.ac.uk/ena/browser/view/PRJEB76464). Intermediate data, i.e. the outputs of `Sourmash gather`, `BLASTn`, `PhyloPhlAn`, `GTDB-Tk`, `CheckM` and `MicrobeAnnotator`, as well as community composition tables (as `phyloseq` objects) required for the R script analyses, are published within this repository.

## Pipeline summary :
![alt text](https://github.com/jorondo1/borealMoss/blob/main/out/Boreal_Moss_WF.jpg)

## Moss artwork by [Isabel Ramirez](https://www.ilustrobiologia.com/portfolio) 
![moss_artwork_isabel_ramirez](https://github.com/jorondo1/borealMoss/blob/main/out/Bryophytes.jpg)
