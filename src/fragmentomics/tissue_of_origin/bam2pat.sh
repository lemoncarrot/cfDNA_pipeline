#!/bin/bash

#wgbstools init_genome GRCh38 -f --fasta_path /mnt/ix1/Resources/GenomeRef/Homo_sapiens/Ensembl/GRCh38_no_alt/Sequence/WholeGenomeFasta/hs38_naa.fna

# conda activate /mnt/ix1/Projects_lite/kchen24/envs/bgzip_env

BAM_DIR="/mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A00_Sample/01_sorted_bam"
OUT_DIR="/mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A00_Sample/b2_meth_pats"

mkdir -p ${OUT_DIR}

wgbstools bam2pat ${BAM_DIR}/*.bam --genome GRCh38 --out_dir ${OUT_DIR}  --nanopore