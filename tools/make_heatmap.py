#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from logsetup import setup_logging
import datetime

p = argparse.ArgumentParser(description="Differential methylation heatmap")
p.add_argument('--matrix', required=True, help="union matrix (tsv or bedgraph)")
p.add_argument('--output-dir', required=True)
p.add_argument('--pvalues', help="corrected_pvalues.txt (cpg mode)")
p.add_argument('--log-dir', help="printing some information about significance")
args = p.parse_args()
os.makedirs(args.output_dir, exist_ok=True)



# logging

logger = setup_logging(args.log_dir, name="heatmap")
working_dir = os.getcwd()

logger.info("=== run start: %s ===", datetime.now().isoformat(timespec="seconds"))

# linear readable fashion

if args.matrix.endswith(".tsv") or args.matrix.endswith(".bedgraph"):
    df = pd.read_csv(args.matrix, sep='\t', header=0)
else:
    raise ValueError("your input matrix extension is wrong!")
    
p_values = pd.read_csv(args.pvalues, sep='\t')
    
# filter for only significant columns

significant = p_values[p_values["FDR"] < 0.05]

print("%d significant features", df.shape[0])

if significant.shape[0] < 2:
    raise ValueError("no significance, exiting without creating heatmap")

key_cols = list(df.columns[1:3])
df = df.merge(significant[key_cols], on=key_cols, how='inner')

logger.info(df.shape[0], " significant features")

# gene and promoter have chr, cpg has chrom
sample_cols = [c for c in df.columns if c not in ["chrom", "chr", "start", "end", "gene"]]
data = df[sample_cols].to_numpy(dtype=float)

# data = np.nan_to_num(data, nan=0.0)

cmap = plt.cm.get_cmap('coolwarm').copy()
cmap.set_bad('lightgray')

fig = plt.figure(figsize=(14, 8))
ax = fig.add_subplot(111)
im = ax.imshow(data, aspect='auto', cmap=cmap,
               interpolation='nearest')
ax.set_xticks(np.arange(len(sample_cols)))
# ax.set_xticklabels(sample_cols, rotation=90, fontsize=8) # don't show 
ax.set(xlabel=None)

ax.set_ylabel("features")
cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
fig.savefig(os.path.join(args.output_dir, "temp_heatmap"), dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"Saved {os.path.join(args.output_dir, 'temp_heatmap')}")