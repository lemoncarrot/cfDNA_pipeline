#!/bin/bash

# conda activate /mnt/ix1/Projects_lite/kchen24/envs/bgzip_env 
# ^^ has installs specified by wgbstools github

PATS_PATH="/mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A00_Sample/b2_meth_pats"
OUT_DIR="/mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A01_Samples/d_fragmentomics/deconv_results"

uxm deconv ${PATS_PATH}/*pat.gz -o ${OUT_DIR}/uxm_output.csv