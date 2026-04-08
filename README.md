# seda_dna_processing  
The workflows are in seda_dna_processing:
- ./main.nf: The very main pipeline to run
- ./workflows/seda_dna_processing.nf: The seda_dna pipeline with tools
- ./subworkflows/local/utils_nfcore_seda_dna_processing_pipeline/input_check.nf: The input check script
- ./modules/local/samplesheet_check.nf
- ./modules/nf-core/: Recent nf-core tools  

The running scripts and outputs are in test_seda_dna:
- ./scripts/run.sh: The script to start the workflows
- ./logs: Including nextflow logs and tmp dir
- ./output: Output files
