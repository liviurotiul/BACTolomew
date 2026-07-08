#!/usr/bin/env Rscript

tree_file <- "/home/bals/Documents/test-Contexar/trees/SNP_trr_2_fixed.nwk"
gpa_file  <- "/home/bals/Documents/test-Contexar3/panaroo/out/gene_presence_absence.csv"
output    <- "~/Documents/dev/BACTolomew/figure_phandango_fullpage.pdf"

# ───────────────────────── diagnostics ─────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(purrr)
  library(lubridate)
  library(stringr)
  library(ggplot2)
})

if (!requireNamespace("ggtree", quietly = TRUE)) {
  stop("ggtree package is required but not installed. Please install it via BiocManager::install('ggtree').")
}
library(ggtree)
library(ape)

cat("Using R version:", paste(R.version$major, R.version$minor, sep = "."), "\n")
cat("Loaded ggtree version:", as.character(packageVersion("ggtree")), "\n")

cat("Tree file: ", tree_file, "\n", sep = "")
cat("GPA file : ", gpa_file, "\n", sep = "")
cat("Output   : ", output, "\n\n", sep = "")

# expand ~ in output path and ensure directory exists
out_path <- path.expand(output)
out_dir  <- dirname(out_path)

if (!file.exists(tree_file)) stop("Tree file does not exist: ", tree_file)
if (!file.exists(gpa_file))  stop("GPA file does not exist: ", gpa_file)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) stop("Could not create output directory: ", out_dir)
}

# ---- load tree ----
tree <- read.tree(tree_file)
tip_order <- tree$tip.label
cat("N tips in tree: ", length(tip_order), "\n", sep = "")

# ---- load Panaroo gene_presence_absence.csv ----
gpa <- read.csv(gpa_file,
                header = TRUE,
                check.names = FALSE,
                stringsAsFactors = FALSE)
cat("N rows in GPA (genes+meta): ", nrow(gpa), "\n", sep = "")
cat("N cols in GPA: ", ncol(gpa), "\n", sep = "")

# gene name column
gene_col <- if ("Gene" %in% names(gpa)) "Gene" else names(gpa)[1]
genes <- gpa[[gene_col]]
cat("Gene column used: ", gene_col, "\n", sep = "")
cat("N genes (rows): ", length(genes), "\n", sep = "")

# ---- match samples between GPA and tree ----
sample_cols <- intersect(colnames(gpa), tip_order)
cat("N sample columns matching tree tips: ", length(sample_cols), "\n", sep = "")
if (length(sample_cols) == 0) {
  cat("Tree tips (first 10):\n")
  print(head(tip_order, 10))
  cat("GPA column names (first 20):\n")
  print(head(colnames(gpa), 20))
  stop("No overlapping sample IDs between tree tip labels and GPA columns.")
}

presence_raw <- gpa[, sample_cols, drop = FALSE]   # genes × samples
presence_chr <- as.matrix(presence_raw)

# binary presence/absence: non-empty → 1, empty/NA → 0
presence_bin <- ifelse(is.na(presence_chr) | presence_chr == "", 0L, 1L)

# samples × genes, initial GPA sample order
gpa_mat <- t(presence_bin)
rownames(gpa_mat) <- sample_cols
colnames(gpa_mat) <- genes

# ---- reorder to tree order ----
common_samples <- intersect(tip_order, rownames(gpa_mat))
cat("N common samples tree∩GPA: ", length(common_samples), "\n", sep = "")
if (length(common_samples) == 0) {
  stop("No common samples between tree tips and GPA matrix after processing.")
}

# drop any tips not in GPA and re-root tree on common samples
tree <- keep.tip(tree, common_samples)
common_samples <- tree$tip.label

# reorder matrix rows to tree tip order
gpa_mat <- gpa_mat[common_samples, , drop = FALSE]

n_samples <- nrow(gpa_mat)
n_genes   <- ncol(gpa_mat)
cat("Final matrix: ", n_samples, " samples × ", n_genes, " genes\n", sep = "")

# ---- build ggtree + gheatmap plot ----

# cladogram: branch.length = "none"
p_tree <- ggtree(tree, layout = "rectangular", branch.length = "none")

# FIX: Removed the manual mapping assignment (p_tree@mapping <- ...)

# gheatmap expects rows = taxa (samples), cols = genes
# we already have that in gpa_mat with rownames matching tip labels
p <- gheatmap(
  p_tree,
  gpa_mat,
  width          = 0.8,   # fraction of total width used by heatmap
  colnames       = FALSE, # don't print gene names (too many)
  color          = NA     # no tile borders
) +
  scale_fill_gradient(
    low  = "white",
    high = "black",
    name = "Presence"
  ) +
  theme(
    legend.position = "none",
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.title.x    = element_text(size = 10),
    plot.margin     = margin(5, 5, 5, 5)
  ) +
  xlab("Gene")  # acts as overall x-label

# ---- save to PDF ----
ggsave(
  filename = out_path,
  plot     = p,
  device   = "pdf",
  width    = 14,
  height   = 8,
  dpi      = 300
)

cat("\nSaved PDF to: ", normalizePath(out_path), "\n", sep = "")
