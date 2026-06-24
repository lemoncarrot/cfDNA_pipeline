#!/bin/bash

nohup nextflow run bams_to_bedgraph.nf \
    --parent_dir /mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A00_Sample \
    --unsorted_bam_dir /mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A00_Sample/00_renamed_bam \
    --ref /mnt/ix1/Resources/GenomeRef/Homo_sapiens/Ensembl/GRCh38_no_alt/Sequence/WholeGenomeFasta/hs38_naa.fna \
    > pipeline.log 2>&1 &

echo "Pipeline running with PID $!"