# cfDNA pipeline overview

### Setting up the environment

1. Clone the repo and enter it:

```bash
git clone https://github.com/lemoncarrot/cfDNA_pipeline
cd cfDNA_pipeline
```

2. Build the conda environment with `mamba` and activate it:

```bash
mamba env create -f tools/cfdna.yml
mamba activate cfdna
```

The pipeline calls the scripts under `tools/` by path, so nothing needs to be added to your `$PATH`. Optionally, `pip install -e .` makes `tools/` importable as a package, but it is not required to run the pipeline.

## 1. Processing initial BAMs

Start with your input BAMs laid out as:

```
<parent>/
  A00_Sample/
    00_unsorted_bam/    <- input *.bam files go here
```

Run the preprocessing pipeline (BAMs -> sorted/indexed -> flagstat -> bedmethyl -> bedgraphs -> QC). Only the input BAM folder and the reference genome the BAMs were aligned to need to be supplied; output folders are derived from the BAM folder's parent and the QC tool is resolved from the repo automatically.

```bash
nextflow run src/preprocessing/bams_to_bedgraph.nf -profile conda \
    --unsorted_bam_dir <parent>/A00_Sample/00_unsorted_bam \
    --ref /path/to/reference.fasta
```

To run in the background and capture logs:

```bash
nohup nextflow run src/preprocessing/bams_to_bedgraph.nf -profile conda \
    --unsorted_bam_dir <parent>/A00_Sample/00_unsorted_bam \
    --ref /path/to/reference.fasta \
    > pipeline.log 2>&1 &
```

Outputs land next to the input folder under `<parent>/A00_Sample/`: `01_sorted_bam/`, `02_flagstat/`, `03_bedmethyl/` (with `5mCG_bedgraph/` and `5hmCG_bedgraph/`), and QC plots in `qc/`.

## 2. Methylation analysis

Takes the per-sample bedgraphs from step 1 (`<parent>/A00_Sample/03_bedmethyl/`) and runs CpG-level differential methylation between cases and controls, for both 5hmCG (`a01`) and 5mCG (`a02`). Each branch: prep bedgraphs -> build a case/control table -> union across samples -> t-test -> FDR correction -> figures.

Intended run (case samples are identified by `--CASE_LABEL` appearing in the filename; everything else is treated as control):

```bash
nextflow run src/methylation_analysis/full_analysis.nf -profile conda \
    --DATA_DIR <parent>/A00_Sample \
    --PARENT_RESULTS_DIR <parent>/A01_Results \
    --REFS_DIR <parent>/refs \
    --CASE_LABEL case
```

Results publish under `PARENT_RESULTS_DIR/` (`a01_diff_cpg_5hmCG/`, `a02_diff_cpg_5mCG/`).

## 3. Fragmentomics analysis
