## Bioinformatic scripts for _Boreal moss-microbe interactions are revealed through metagenome assembly of novel bacterial species_

Published in Nature Scientific Reports : [Ishak et al, 2024](https://www.nature.com/articles/s41598-024-73045-z)
Raw data available at ENA under project accession PRJEB76464

The `moss_pipeline.sh` script contains the main steps for metagenome processing, assembly, and profiling. The pipeline in not seamless and I highly recommend to run it step by step. These steps rely on modules of a pipeline implemented by our team [available here](https://github.com/jflucier/ILL_pipelines) and are scripted to work with the SLURM workload manager. Several steps use singularity-contained programs to process data in parallel, but the underlying scripts (in the [ILL_pipelines repository](https://github.com/jflucier/ILL_pipelines)) can be run on a sample-wise basis.

The `moss_pipeline.sh` script was run on the Mammouth-mp2 clusters, whereas R scripts were run locally on a Macbook Pro M2. Raw data (metagenomes and assemblies) required for the bash script will be available on ENA (BioProject `PRJEB76464`). Intermediate data, i.e. the `Sourmash gather`, `BLASTn`, `PhyloPhlAn`, `GTDB-Tk`, `CheckM` and `MicrobeAnnotator` outputs, as well as community composition tables (as `phyloseq` objects) required for the R script analyses, are published within this repository. 

Pipeline summary :
![alt text](https://github.com/jorondo1/borealMoss/blob/main/out/Boreal_Moss_WF.jpg)
