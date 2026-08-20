#!/bin/bash

# Specify our project
#SBATCH --account=acad

# Request resources
# 2 processor or task
# 2 gigabytes of memory
# 45 minutes of walltime
#SBATCH --nodes=1
#SBATCH --mem=2GB
#SBATCH --time=20:00:00

# Exporting none of the login environment
#SBATCH --export=None

#Place the slurm output file in a dir
#SBATCH --output=/hpcfs/users/a1844642/project/test_seda_dna/logs/slurm-%j.out
#SBATCH --error=/hpcfs/users/a1844642/project/test_seda_dna/logs/slurm-%j.err

# Specify the partition
#SBATCH --partition=batch

module load Nextflow/25.10.2

export TMPDIR="/hpcfs/users/a1844642/project/test_seda_dna/logs/tmp_dir"
export SINGULARITY_CACHEDIR="/hpcfs/users/a1844642/project/test_seda_dna/logs/cache_dir"
export SINGULARITY_LIBRARYDIR="/hpcfs/users/a1844642/project/test_seda_dna/logs/library_dir"
export NXF_SINGULARITY_CACHEDIR="/hpcfs/users/a1844642/project/test_seda_dna/logs/nxf_cachedir"
export NXF_SINGULARITY_LIBRARYDIR="/hpcfs/users/a1844642/project/test_seda_dna/logs/nxf_librarydir"

nextflow -log '/hpcfs/users/a1844642/project/test_seda_dna/logs/main.nextflow.log' \
    run /hpcfs/users/a1844642/project/seda_dna_processing \
    -resume \
    -profile conda \
    --outdir '/hpcfs/users/a1844642/project/test_seda_dna/output' \
    --input '/hpcfs/users/a1844642/project/test_seda_dna/samplesheet_vilma_data.csv' \
    -c '/hpcfs/users/a1844642/project/seda_dna_processing/nextflow.config' \
    -w '/hpcfs/users/a1844642/project/test_seda_dna/logs/work'

#--use_adna_trim \
#--adna_trim_path /hpcfs/groups/acad_users/shyrav/resources/adna/adna-trim \