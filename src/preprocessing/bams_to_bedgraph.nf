#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------------------------------------------------------------------------
// Inputs. Only --unsorted_bam_dir and --ref must be supplied; everything else
// is derived from them or auto-resolved relative to this pipeline file.
//
//   nextflow run bams_to_bedgraph.nf -profile conda \
//       --unsorted_bam_dir <parent>/A00_Samples/00_unsorted_bam \
//       --ref /path/to/reference.fasta
// ---------------------------------------------------------------------------

// Folder containing the input *.bam files (required).
params.unsorted_bam_dir = null
// Reference genome FASTA the BAMs were aligned to (required by modkit pileup).
params.ref              = null
// Root for numbered output folders (01_sorted_bam, 02_flagstat, ...).
// Defaults to the parent of the input BAM folder.
params.parent_dir       = params.unsorted_bam_dir ? file(params.unsorted_bam_dir).parent.toString() : null
// Where QC plots are published. Defaults to <parent_dir>/qc.
params.qc_out_dir       = params.parent_dir ? "${params.parent_dir}/qc" : null
// Analysis scripts shipped in the repo. Auto-resolved relative to this .nf:
// src/preprocessing/ -> ../../tools. Override only if running tools elsewhere.
params.tools_path       = "${projectDir}/../../tools"

if( !params.unsorted_bam_dir )
    error "Missing --unsorted_bam_dir: folder containing the input *.bam files"
if( !params.ref )
    error "Missing --ref: reference genome FASTA the BAMs were aligned to"

def qc_script = file("${params.tools_path}/alignment_metrics.py")
if( !qc_script.exists() )
    error "Could not resolve QC script at ${qc_script}. Pass --tools_path <repo>/tools"

process SORT_INDEX {
    tag "${sample_id}"

    conda 'bioconda::samtools=1.19.2'

    publishDir "${params.parent_dir}/01_sorted_bam", mode: 'copy'

    input:
    tuple val(sample_id), path(unsorted_bam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.csi"), emit: sorted_bam

    script:
    """
    samtools sort \
        --write-index \
        "$unsorted_bam" \
        -o ${sample_id}_sorted.bam
    """
}

process RUN_FLAGSTAT {
    tag "${sample_id}"

    conda 'bioconda::samtools=1.19.2'

    publishDir "${params.parent_dir}/02_flagstat", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam), path(csi)

    output:
    tuple val(sample_id), path("${sample_id}_flagstat.txt"), emit: flagstat

    script:
    """
    samtools flagstat "${sorted_bam}" > ${sample_id}_flagstat.txt
    """
}

process MAKE_BEDS {
    tag "${sample_id}"

    conda 'bioconda::modkit=0.6.1'

    publishDir "${params.parent_dir}/03_bedmethyl", mode: 'copy', pattern: "*.bed"
    publishDir "${params.parent_dir}/03_bedmethyl/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample_id), path(sorted_bam), path(csi)

    output:
    tuple val(sample_id), path("${sample_id}.5mC5hmC.bed"), emit: bed
    path "${sample_id}.modkit.log", emit: log

    script:
    """
    modkit pileup \
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

process RUN_QC {
    tag "qc"

    conda 'conda-forge::pandas conda-forge::numpy conda-forge::matplotlib bioconda::samtools=1.19.2'

    publishDir "${params.qc_out_dir}", mode: 'copy', pattern: "*.png"

    input:
    path qc_script
    path "5mCG_bedgraph/*"
    path "5hmCG_bedgraph/*"

    output:
    path "*.png", emit: qc_plots

    script:
    // bedgraphs are staged into per-modification subdirs so the script's
    // suffix-based directory listing stays clean; flagstat/sorted-bam paths
    // feed the optional alignment-metrics section and point at publishDir.
    // PNGs are written to the task dir and copied out by publishDir.
    """
    python ${qc_script} \
        --flagstat-folder-path ${params.parent_dir}/02_flagstat \
        --sorted-bam-path ${params.parent_dir}/01_sorted_bam \
        --bedgraph-path 5mCG_bedgraph \
        --h-bedgraph-path 5hmCG_bedgraph \
        --output-dir .
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

    // gate QC on all samples' bedgraphs so it runs once at the end
    RUN_QC(
        qc_script,
        MAKE_BEDGRAPHS.out.mC_bedgraph.map { it[1] }.collect(),
        MAKE_BEDGRAPHS.out.hmC_bedgraph.map { it[1] }.collect()
    )
}