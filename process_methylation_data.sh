#!/bin/bash
## Perform DNA methylation array pre-processing
source ~/.bashrc
conda activate methylation_processing

## Perform sample QC and filtering
#Rscript 1_sampleQC.R --help
Rscript 1_sampleQC.R \
        --idat_dir "/vault/veeramah/stress_aging/EPIC144_pilot/RawData/" \
        --sample_sheet_dir "/vault/veeramah/stress_aging/EPIC144_pilot/SampleSpecificFiles/" \
        --sample_sheet_csv "kVeeramah_EPIC144_SampleSheet.csv" \
        --sample_manifest "/vault/veeramah/stress_aging/sample_manifests/21788_21785_2_ext_ak_manifest_kreid_1__1_.csv" \
        --drop_sex_mismatches FALSE \
        --output "post_step1"

## Perform probe QC and filtering
#Rscript 2_probeQC.R --help
Rscript 2_probeQC.R \
        --input "post_step1.RData" \
        --cross_reactive_probes "/vault/veeramah/people/shyamie/methylation/meth_datasets/Pidsley_etal_2016/unique_cross-reactive_probes.RData" \
        --exclude_sex_chrs TRUE \
        --output "post_step2"

## Perform correction, normalization, and optionally, cell count and age estimation
#Rscript 3_normalization.R --help
Rscript 3_normalization.R \
        --input "post_step2.RData" \
        --estimate_WBC_props TRUE \
        --estimate_age TRUE \
        --set_poor_CpGs_to_NA TRUE \
        --output "pilot_144_blood"
