#!/bin/bash -l
#SBATCH -A naiss2023-22-1332
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 24:00:00

 module load bioinfo-tools Nextflow nf-core nf-core-pipelines
export NXF_HOME=/crex/proj/uppstore2017150/private/eneritz/SICILIAN/nCT_MK_pipeline/atac_mm10_CZBH_2024/

 cd /crex/proj/uppstore2017150/private/eneritz/SICILIAN/nCT_MK_pipeline/atac_mm10_CZBH_2024/

nextflow run -resume nf-core/atacseq  --outdir outs  -profile uppmax  --input /crex/proj/uppstore2017150/private/eneritz/SICILIAN/nCT_MK_pipeline/atac_mm10_CZBH_2024/sample_sheet.csv --genome GRCm38   --read_len
gth 50 --project naiss2023-22-1332

