#!/bin/bash

# Specify our project
#SBATCH --account=acad

# Request resources
# 2 processor or task
# 2 gigabytes of memory
# 45 minutes of walltime
#SBATCH --nodes=1
#SBATCH --mem=2GB
#SBATCH --time=4:00:00

# Exporting none of the login environment
#SBATCH --export=None

#Place the slurm output file in a dir
#SBATCH --output=/hpcfs/users/a1844642/project/test_seda_dna/logs/slurm-%j.out

# Specify the partition
#SBATCH --partition=batch

module load Nextflow/25.10.2

export SINGULARITY_CACHEDIR="/hpcfs/users/a1844642/project/seda_dna_processing/logs/cache_dir"
export SINGULARITY_LIBRARYDIR="/hpcfs/users/a1844642/project/seda_dna_processing/logs/library_dir"
export NXF_SINGULARITY_CACHEDIR="/hpcfs/users/a1844642/project/seda_dna_processing/logs/nxf_cachedir"
export NXF_SINGULARITY_LIBRARYDIR="/hpcfs/users/a1844642/project/seda_dna_processing/logs/nxf_librarydir"

nextflow run /hpcfs/users/a1844642/project/seda_dna_processing/main.nf \
    -profile singularity \
    --outdir '../output' \
    --input '../samplesheet.csv' \
    -w '../logs/work'