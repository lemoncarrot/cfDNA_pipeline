#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.parent_dir = "/path/to/parent_dir"
params.unsorted_bam_dir = "/path"
params.ref = "/path/to/ref"

params.samtools = "/mnt/ix1/Resources/tools/samtools/v1.19.2/bin/samtools"
params.modkit = "/mnt/ix1/Resources/tools/modkit/v0.6.1/modkit"

process SORT_INDEX {
    tag "${sample_id}"

    publishDir "${params.parent_dir}/01_sorted_bam", mode: 'copy'

    input:
    tuple val(sample_id), path(unsorted_bam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.csi"), emit: sorted_bam

    script:
    """
    ${params.samtools} sort \
        --write-index \
        "$unsorted_bam" \
        -o ${sample_id}_sorted.bam
    """
}

process RUN_FLAGSTAT {
    tag "${sample_id}"

    publishDir "${params.parent_dir}/02_flagstat", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam), path(csi)

    output:
    tuple val(sample_id), path("${sample_id}_flagstat.txt"), emit: flagstat

    script:
    """
    ${params.samtools} flagstat "${sorted_bam}" > ${sample_id}_flagstat.txt
    """
}

process MAKE_BEDS {
    tag "${sample_id}"

    publishDir "${params.parent_dir}/03_bedmethyl", mode: 'copy', pattern: "*.bed"
    publishDir "${params.parent_dir}/03_bedmethyl/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample_id), path(sorted_bam), path(csi)

    output:
    tuple val(sample_id), path("${sample_id}.5mC5hmC.bed"), emit: bed
    path "${sample_id}.modkit.log", emit: log

    script:
    """
    ${params.modkit} pileup \
        "${sorted_bam}" \
        "${sample_id}.5mC5hmC.bed" \
        --cpg \
        --combine-strands \
        --ref "${params.ref}" \
        --modified-bases 5mC 5hmC \
        -t ${task.cpus} \
        --log-filepath "${sample_id}.modkit.log"
    """
}

process MAKE_BEDGRAPHS {
    tag "${sample_id}"

    publishDir "${params.parent_dir}/03_bedmethyl/5hmCG_bedgraph", mode: 'copy', pattern: "*.5hmC.sorted.bedgraph"
    publishDir "${params.parent_dir}/03_bedmethyl/5mCG_bedgraph", mode: 'copy', pattern: "*.5mC.sorted.bedgraph"

    input:
    tuple val(sample_id), path(bed)

    output:
    tuple val(sample_id), path("${sample_id}.5hmC.sorted.bedgraph"), emit: hmC_bedgraph
    tuple val(sample_id), path("${sample_id}.5mC.sorted.bedgraph"), emit: mC_bedgraph

    script:
    """
    awk -F"\t" -v mod="h" '
    \$4==mod && \$5>0 {
        print \$1"\t"\$2"\t"\$3"\t"\$11
    }' "${bed}" | sort -k1,1 -k2,2n > "${sample_id}.5hmC.sorted.bedgraph"

    awk -F"\t" -v mod="m" '
    \$4==mod && \$5>0 {
        print \$1"\t"\$2"\t"\$3"\t"\$11
    }' "${bed}" | sort -k1,1 -k2,2n > "${sample_id}.5mC.sorted.bedgraph"
    """
}

workflow {
    bam_ch = Channel
        .fromPath("${params.unsorted_bam_dir}/*.bam", glob: true)
        .map { bam -> tuple(bam.baseName, bam) }

    SORT_INDEX(bam_ch)
    RUN_FLAGSTAT(SORT_INDEX.out.sorted_bam)
    MAKE_BEDS(SORT_INDEX.out.sorted_bam)
    MAKE_BEDGRAPHS(MAKE_BEDS.out.bed)
}