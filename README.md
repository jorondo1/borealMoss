# borealMoss
Bioinformatic scripts for _Boreal moss-microbe interactions are revealed through metagenome assembly of novel bacterial species_
Preprint: https://doi.org/10.1101/2023.04.06.535926

All bash (.sh) scripts were run on the Mammouth-mp2 clusters, whereas R scripts were run locally on a Macbook Pro M2. The moss_pipeline.sh script contains the main steps for metagenome processing, assembly, and profiling. These steps rely on modules of a pipeline implemented by our team [available here](https://github.com/jflucier/ILL_pipelines) and processed using the SLURM workload manager. Several steps use singularity-contained programs to process data in parallel, but the underlying scripts (in the [ILL_pipelines repository](https://github.com/jflucier/ILL_pipelines)) can be run on a sample-wise basis.

