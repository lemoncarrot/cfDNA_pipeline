#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(logger)
  library(ComplexHeatmap)
  library(circlize)
})

# usage:
# Rscript make_heatmap.R --matrix /mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A01_Results/a02_diff_cpg_5mCG/all.union.sorted.bedgraph --output-dir /mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A01_Results/a02_diff_cpg_5mCG/figs --pvalues /mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A01_Results/a02_diff_cpg_5mCG/corrected_pvalues.txt
# Rscript make_heatmap.R --matrix /mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A01_Results/a01_diff_cpg_5mCG/all.union.sorted.bedgraph --output-dir /mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A01_Results/a01_diff_cpg_5mCG/figs --pvalues /mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A01_Results/a01_diff_cpg_5hmCG/corrected_pvalues.txt

# ---- args ----
option_list <- list(
  make_option("--matrix", type = "character", help = "union matrix (tsv or bedgraph)"),
  make_option("--output-dir", type = "character"),
  make_option("--pvalues", type = "character", default = NULL, help = "corrected_pvalues.txt (cpg mode)")
)
opt <- parse_args(OptionParser(option_list = option_list))

# for testing:
# opt$matrix = "/mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A01_Results/a02_diff_cpg_5mCG/all.union.sorted.bedgraph"
# opt$output_dir = "/mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A01_Results/a02_diff_cpg_5mCG"
# opt$pvalues = "/mnt/ix1/Projects_lite/kchen24/cfDNA_CRC/A01_Results/a02_diff_cpg_5mCG/corrected_pvalues.txt"

if (is.null(opt$matrix) || is.null(opt$`output-dir`)) {
  stop("--matrix and --output-dir are required")
}
dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)

log_info("=== run start: {format(Sys.time(), '%Y-%m-%dT%H:%M:%S')} ===")

# ---- read matrix ----
if (grepl("\\.(tsv|bedgraph)$", opt$matrix)) {
  df <- fread(opt$matrix, sep = "\t", header = TRUE, data.table = FALSE)
} else {
  stop("your input matrix extension is wrong!")
}

# ---- read pvalues + filter significant ----
p_values <- fread(opt$pvalues, sep = "\t", header = TRUE, data.table = FALSE)
filtered_nas <- p_values %>% filter(!is.na(p_values$FDR))
significant <- filtered_nas[filtered_nas$FDR < 0.05, ]
log_info("{nrow(significant)} significant features in pvalue table")

# however, also could be potential significant features that have high/low expression in case but no coverage in control and vice versa
# think about how to define this here

if (nrow(significant) < 2) {
  stop("no significance, exiting without creating heatmap")
}

# ---- keep matrix rows at significant positions ----
key_cols <- colnames(df)[1:3]
df <- inner_join(df, significant, by = key_cols)
log_info("{nrow(df)} matrix rows after restricting to significant positions")

# ---- build numeric matrix ----
meta_cols <- c("chrom", "chr", "start", "end", "gene") # i think "chr" naming happens for gene/promoter aggregation
sample_cols <- setdiff(colnames(df), meta_cols)
data <- as.matrix(df[, sample_cols])
storage.mode(data) <- "double" # equivalent of dtype=float

# matrix has methylation information but also stuff from p-values
# make new column "hyper" or "hypo" based on mean_case and mean_control
data <- as.data.frame(data)
data$group <- ifelse(data$mean_case > data$mean_control, "hyper", "hypo") # everything already significant positions

make_heatmap <- function(X, min_control_low, min_control_high, min_case_low, min_case_high, top_n) {
  # X has sample_cols but extra t-statistic, pval, mean_case, mean_control, n_case, n_control, FDR, group columns from p-vals

  log_info("Running heatmap with coverage filters — min_control_low: {min_control_low}, min_control_high: {min_control_high}, min_case_low: {min_case_low}, min_case_high: {min_case_high}")

  meta_cols <- c(
    "t-statistic", "pval", "mean_case", "mean_control",
    "n_case", "n_control", "FDR", "group"
  )
  sample_cols <- setdiff(colnames(X), meta_cols)

  is_case_col <- grepl("case", sample_cols, ignore.case = TRUE)
  is_ctrl_col <- grepl("ctrl", sample_cols, ignore.case = TRUE)

  col_group <- ifelse(is_case_col, "case", "ctrl")
  col_group <- factor(col_group, levels = c("ctrl", "case"))

  # heatmap should have 4 sections (hypo first, then hyper), (in each category, first samples passing min_control_low and min_case_high, then samples passing
  # min_control_high, min_case_low), this is to make sure we see the pattern of cpgs that are sparse in case / sparse in control

  X <- X %>%
    mutate(subgroup = case_when(
      n_control >= min_control_low & n_case >= min_case_high ~ "A",
      n_control >= min_control_high & n_case >= min_case_low ~ "B",
      TRUE ~ NA_character_
    ))


  size_before <- dim(X)[1]
  print("Dropping all rows that do not pass the coverage filter!")
  X <- X %>% filter(!is.na(subgroup))
  size_after <- dim(X)[1]
  print(paste0("dropped ", size_before - size_after, " cpg for not passing the coverage filter!"))

  X <- X %>%
    mutate(
      section = factor(
        paste(group, subgroup, sep = "_"),
        levels = c("hypo_A", "hypo_B", "hyper_A", "hyper_B")
      )
    )

  section_counts <- table(X$section)
  log_info("Rows per section: {paste(names(section_counts), section_counts, sep = '=', collapse = ', ')}")


  # only plot top_n per group
  X_plot <- X %>%
    group_by(section) %>%
    arrange(pval) %>%
    slice_head(n = top_n) %>%
    ungroup()


  # reorder rows
  X_plot <- X_plot %>% arrange(section)
  X_plot <- X_plot %>% mutate(group_label = factor(group, levels = c("hypo", "hyper")))

  # using beautiful complexheatmap package here, annotating the rows appropriately as specified above, cluster within the groups if possible
  # rows should be ordered exactly as in X (already control, then case, add a small divider between them), "ctrl" and "case" within column names to specify

  # filter down to plotting matrix)
  mat <- as.matrix(X_plot[, sample_cols])

  # annotation
  col_anno <- HeatmapAnnotation(
    group = col_group,
    col = list(group = c(ctrl = "#4575B4", case = "#D73027")),
    annotation_name_side = "left"
  )

  print(which(mat < 0, arr.ind = TRUE))

  row_anno <- rowAnnotation(
    section = X_plot$group_label,
    col = list(section = c(
      hypo  = "#74ADD1",
      hyper = "#F46D43"
    )),
    annotation_name_side = "top"
  )

  # set color to handle na values
  mat_range <- range(mat, na.rm = TRUE, finite = TRUE)
  col_fun <- colorRamp2(
    c(mat_range[1], mean(mat_range), mat_range[2]),
    c("blue", "white", "red")
  )

  ht <- Heatmap(
    mat,
    name = "Methylation Heatmap",
    na_col = "grey90",
    col = col_fun,
    top_annotation = col_anno,
    left_annotation = row_anno,
    row_split = X_plot$group_label,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE, # keeps the 4 section blocks in our specified order (no reordering of slices)
    cluster_columns = FALSE, # preserve X's column order exactly
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_names_gp = grid::gpar(fontsize = 8),
    row_title_rot = 0,
    border = TRUE,
    column_gap = unit(1, "mm")
  )

  out_path <- file.path(opt$`output-dir`, "methylation_heatmap.pdf")
  pdf(out_path, width = 10, height = 8)
  draw(ht)
  dev.off()
  log_info("heatmap saved to {out_path}")

  # after plotting, return the df used to plot
  return(X)
}

temp <- make_heatmap(data, 2, 4, 2, 8, 100)

# sort by p-value, then greatest absolute diff in mean_case and mean_control
# take top_n to plot


# now have data heatmap to work with, add some metadata columns to make it pretty
# - stratify hyper and hypo (mean_case and mean_control columns in p_values), stratify by passing separate threshold filters (n_case and n_control columns in p_values)
# can have heatmap with topn, filtered first for significance, then filtered by difference of means, p-value, etc.

# - annotated hyper & hypo, low coverage in what category, etc.

# after heatmap, use the same filtering conditions to create another volcano
