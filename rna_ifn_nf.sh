#!/bin/bash -l
#SBATCH -A naiss2023-22-1332
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 24:00:00

 module load bioinfo-tools Nextflow nf-core nf-core-pipelines
export NXF_HOME=/proj/naiss2023-23-112/private/eneritz/bulk_rna_ifng_CZ/

 cd /proj/naiss2023-23-112/private/eneritz/bulk_rna_ifng_CZ

nextflow run nf-core/rnaseq  --outdir outs  -profile uppmax  --input /proj/naiss2023-23-112/private/eneritz/bulk_rna_ifng_CZ/sample_sheet.csv --genome GRCm38   --project naiss2023-22-1332

