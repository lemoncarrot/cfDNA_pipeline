import argparse
import os
import re
import subprocess
from io import StringIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

GENOME_CPG_COUNT = 28991704

# this can be chained onto the end of bams_to_bedgraph.nf

"""
run using:

python alignment_metrics.py --flagstat-folder-path /mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A00_Sample/02_flagstat --sorted-bam-path /mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A00_Sample/01_sorted_bam --bedgraph-path /mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A00_Sample/03_bedmethyl/5mCG_bedgraph --h-bedgraph-path /mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A00_Sample/03_bedmethyl/5hmCG_bedgraph --output-dir /mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A01_Results/qc
"""

p = argparse.ArgumentParser()
p.add_argument('--flagstat-folder-path', required=True)
p.add_argument('--sorted-bam-path', required=True)
p.add_argument('--bedgraph-path', required=True)
p.add_argument('--h-bedgraph-path', required=True)
p.add_argument('--output-dir', required=True) 
args = p.parse_args()

# sample names must match between flagstat and sorted bams

# folder = args.flagstat_folder_path
# 
# files = os.listdir(folder)
# 
# def get_primary_mapped(flagstat_path):
#     with open(flagstat_path) as f:
#         for line in f:
#             if "primary mapped" in line:
#                 count = int(line.split()[0])
#                 m = re.search(r"\(([\d.]+)%", line)
#                 percent = float(m.group(1)) if m else None
#                 return count, percent
#     return None, None
# 
# def get_mean_depth(bam_path):
#     res = subprocess.run(["samtools", "coverage", bam_path],
#                          capture_output=True, text=True, check=True)
#     cov = pd.read_csv(StringIO(res.stdout), sep="\t")
#     lengths = cov["endpos"] - cov["startpos"] + 1
#     return float((cov["meandepth"] * lengths).sum() / lengths.sum())
# 
# 
# rows = []
# for file in files:
#     if not file.endswith("_flagstat.txt"):
#         raise ValueError("SOMething wrong here!!")
#     sample_id = file.replace("_flagstat.txt", "")
#     full_path = os.path.join(folder, file)
#     
#     count, percent = get_primary_mapped(full_path)
#     bam_path = os.path.join(args.sorted_bam_path, f"{sample_id}_sorted.bam")
#     coverage = get_mean_depth(bam_path) if os.path.exists(bam_path) else None
#     if coverage is None:
#         raise ValueError("Something is critically wrong here!")
#     
#     rows.append({"sample_id": sample_id,
#                  "primary_mapped": count,
#                  "mapped_pct": percent,
#                  "mean_depth": coverage})
#     
# qc = pd.DataFrame(rows).set_index("sample_id")
# print(qc.head())
# 
# # using this dataframe, plot metrics stratified by sample group
# 
# qc["group"] = np.where(qc.index.str.contains("case"), "case", "control")
# print(qc["group"].value_counts())
# 
# no_outliers_qc = qc[(qc["primary_mapped"] < 20_000_000) & (qc["mapped_pct"] > 60) & (qc["mean_depth"] < 1.5)]
# 
# metrics = ["mapped_pct", "primary_mapped", "mean_depth"]
# labels  = ["Primary mapped %", "Primary mapped reads", "Mean depth (X)"]
# order   = ["control", "case"]
# 
# fig, axes = plt.subplots(3, 1, figsize=(6, 11))
# 
# for ax, metric, label in zip(axes, metrics, labels):
#     groups = [qc.loc[qc["group"] == g, metric].dropna() for g in order]
#     ax.boxplot(groups, labels=order, showfliers=True, widths=0.5)
# 
#     # jitter the individual samples over the box — small n, you want to see points
#     for i, vals in enumerate(groups, start=1):
#         x = np.random.normal(i, 0.04, size=len(vals))
#         ax.scatter(x, vals, s=18, alpha=0.6, color="black", zorder=3)
# 
#     ax.set_ylabel(label)
#     ax.set_title(label, fontsize=10)
# 
# axes[-1].set_xlabel("group")
# fig.tight_layout()
# fig.savefig(os.path.join(args.output_dir, "qc_boxplots.png"), dpi=300)
# 
# # now for no outliers
# 
# fig, axes = plt.subplots(3, 1, figsize=(6, 11))
# 
# for ax, metric, label in zip(axes, metrics, labels):
#     groups = [no_outliers_qc.loc[no_outliers_qc["group"] == g, metric].dropna() for g in order]
#     ax.boxplot(groups, labels=order, showfliers=True, widths=0.5)
# 
#     # jitter the individual samples over the box — small n, you want to see points
#     for i, vals in enumerate(groups, start=1):
#         x = np.random.normal(i, 0.04, size=len(vals))
#         ax.scatter(x, vals, s=18, alpha=0.6, color="black", zorder=3)
# 
#     ax.set_ylabel(label)
#     ax.set_title(label, fontsize=10)
# 
# axes[-1].set_xlabel("group")
# fig.tight_layout()
# fig.savefig(os.path.join(args.output_dir, "no_outliers_qc_boxplots.png"), dpi=300)


# looking at total sites covered and methylation percentage per sample
files = os.listdir(args.bedgraph_path)

rows = []
for file in files:
    if not file.endswith("5mC.sorted.bedgraph"):
        raise ValueError("Something has gone terribly wrong!")
    
    sample_id = file.replace(".5mC.sorted.bedgraph", "")
    full_path = os.path.join(args.bedgraph_path, file)
    
    temp = subprocess.run(["wc", "-l", full_path], capture_output=True, text=True)
    total_cpg_count = int(temp.stdout.split()[0])
    
    methylated_count = 0
    with open(full_path) as f:
        for line in f:
            fields = line.split("\t")
            if len(fields) >=4 and float(fields[3]) > 0:
                methylated_count += 1
                
    methylated_percent = methylated_count / total_cpg_count
                
    rows.append({"sample_id": sample_id,
                  "total_cpg_count": total_cpg_count,
                  "methylated_percentage": methylated_percent})
    
cpg_qc = pd.DataFrame(rows).set_index("sample_id")
print(cpg_qc.head())

# same as above, but for 5hmC
files = os.listdir(args.h_bedgraph_path)

rows = []
for file in files:
    if not file.endswith("5hmC.sorted.bedgraph"):
        raise ValueError("Something has gone terribly wrong!")

    sample_id = file.replace(".5hmC.sorted.bedgraph", "")
    full_path = os.path.join(args.h_bedgraph_path, file)

    temp = subprocess.run(["wc", "-l", full_path], capture_output=True, text=True)
    total_cpg_count = int(temp.stdout.split()[0])

    methylated_count = 0
    with open(full_path) as f:
        for line in f:
            fields = line.split("\t")
            if len(fields) >=4 and float(fields[3]) > 0:
                methylated_count += 1

    methylated_percent = methylated_count / total_cpg_count

    rows.append({"sample_id": sample_id,
                  "total_cpg_count": total_cpg_count,
                  "methylated_percentage": methylated_percent})

hcpg_qc = pd.DataFrame(rows).set_index("sample_id")
print(hcpg_qc.head())

# combined plot: total CpGs covered, 5mC %, 5hmC %
cpg_qc["group"] = np.where(cpg_qc.index.str.contains("case"), "case", "control")
hcpg_qc["group"] = np.where(hcpg_qc.index.str.contains("case"), "case", "control")
print(cpg_qc["group"].value_counts())

order = ["control", "case"]
# (dataframe, column, panel label) -- total CpGs from the 5mC table (same
# covered sites as 5hmC), then the two modification percentages
panels = [
    (cpg_qc,  "total_cpg_count",       "CpG sites covered"),
    (cpg_qc,  "methylated_percentage", "Fraction of genomic CpGs methylated (5mC)"),
    (hcpg_qc, "methylated_percentage", "Fraction of genomic CpGs hydroxymethylated (5hmC)"),
]

fig, axes = plt.subplots(3, 1, figsize=(6, 11))

for ax, (df, metric, label) in zip(axes, panels):
    groups = [df.loc[df["group"] == g, metric].dropna() for g in order]
    ax.boxplot(groups, labels=order, showfliers=True, widths=0.5)

    # jitter the individual samples over the box — small n, you want to see points
    for i, vals in enumerate(groups, start=1):
        x = np.random.normal(i, 0.04, size=len(vals))
        ax.scatter(x, vals, s=18, alpha=0.6, color="black", zorder=3)

    ax.set_ylabel(label)
    ax.set_title(label, fontsize=10)

axes[-1].set_xlabel("group")
fig.tight_layout()
fig.savefig(os.path.join(args.output_dir, "cpg_methylation_qc_boxplots.png"), dpi=300)