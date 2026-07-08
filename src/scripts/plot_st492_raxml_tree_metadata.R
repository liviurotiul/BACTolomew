#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(ggtree)
  library(readr)
  library(stringr)
  library(tibble)
})

run_dir <- "/home/bals/Documents/dev/BACTolomew/runs/sample_78eaa9126850/run_78eaa9126850"
tree_path <- file.path(run_dir, "gubbins", "core_gene.final_tree.tre")
metadata_path <- "/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.csv"
fastgear_path <- file.path(run_dir, "fastgear", "sample_clades.tsv")
recent_recombination_path <- file.path(run_dir, "fastgear", "output", "recombinations_recent.txt")
output_dir <- "/Data/ST492/figures2026"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#132033"
branch_col <- "#111827"

normalise_country <- function(x) {
  y <- as.character(x) %>% str_squish() %>% na_if("")
  case_when(
    is.na(y) | str_to_lower(y) %in% c("unknown", "na", "n/a", "not available") ~ "Unknown",
    y %in% c("United Kingdom", "UK", "Scotland", "UK (Scotland)") ~ "UK (Scotland)",
    y %in% c("United States", "United States of America", "U.S.A.", "US") ~ "USA",
    TRUE ~ y
  )
}

format_complex_sort_key <- function(x) {
  y <- suppressWarnings(as.integer(round(as.numeric(x))))
  if_else(is.na(y), "ZZZ_missing", sprintf("CT%05d", y))
}

read_recent_recombination <- function(path) {
  if (!file.exists(path)) {
    return(tibble(
      Sample = character(),
      recent_recombination_segments = integer(),
      recent_recombination_bp = numeric(),
      recent_recombination_high_conf_segments = integer()
    ))
  }

  read_table(
    path,
    skip = 1,
    col_types = cols(
      Start = col_double(),
      End = col_double(),
      DonorLineage = col_character(),
      RecipientStrain = col_character(),
      `log(BF)` = col_double(),
      StrainName = col_character()
    )
  ) %>%
    transmute(
      Sample = StrainName,
      segment_bp = pmax(0, End - Start + 1),
      log_bf = `log(BF)`
    ) %>%
    filter(!is.na(Sample), Sample != "") %>%
    group_by(Sample) %>%
    summarise(
      recent_recombination_segments = n(),
      recent_recombination_bp = sum(segment_bp, na.rm = TRUE),
      recent_recombination_high_conf_segments = sum(log_bf >= 10, na.rm = TRUE),
      .groups = "drop"
    )
}

metadata <- read_csv(
  metadata_path,
  col_types = cols(
    .default = col_guess(),
    Sample = col_character(),
    country = col_character(),
    complex_type = col_double()
  )
) %>%
  transmute(
    Sample,
    subtype = suppressWarnings(as.integer(subtype)),
    country = normalise_country(country),
    complex_order = format_complex_sort_key(complex_type)
  )

fastgear <- read_tsv(
  fastgear_path,
  col_types = cols(
    Sample = col_character(),
    FastGearLineage = col_integer(),
    FastGearCluster = col_integer(),
    FastGearStrainIndex = col_integer(),
    FastGearMethod = col_character()
  )
) %>%
  transmute(
    Sample,
    fastgear_lineage = if_else(is.na(FastGearLineage), "Lineage NA", paste0("Lineage ", FastGearLineage)),
    fastgear_cluster = if_else(is.na(FastGearCluster), "Cluster NA", paste0("Cluster ", FastGearCluster))
  )

recent_recombination <- read_recent_recombination(recent_recombination_path)

tree <- read.tree(tree_path)
reference_dropped <- "Reference" %in% tree$tip.label
if (reference_dropped) {
  tree <- ape::drop.tip(tree, "Reference")
}

lineage1_outgroup <- fastgear %>%
  filter(fastgear_lineage == "Lineage 1") %>%
  pull(Sample) %>%
  intersect(tree$tip.label)
rooted_on_lineage1 <- FALSE
if (length(lineage1_outgroup) > 1 && ape::is.monophyletic(tree, lineage1_outgroup)) {
  tree <- ape::root(tree, outgroup = lineage1_outgroup, resolve.root = TRUE)
  rooted_on_lineage1 <- TRUE
}

annotation <- tibble(Sample = tree$tip.label) %>%
  left_join(metadata, by = "Sample") %>%
  left_join(fastgear, by = "Sample") %>%
  left_join(recent_recombination, by = "Sample") %>%
  mutate(
    country = coalesce(country, "Unknown"),
    subtype = if_else(is.na(subtype), -1L, subtype),
    subtype_label = if_else(subtype == -1L, "ST NA", paste0("ST", subtype)),
    is_st492 = subtype == 492L,
    complex_order = coalesce(complex_order, "ZZZ_missing"),
    fastgear_lineage = coalesce(fastgear_lineage, "Lineage NA"),
    fastgear_cluster = coalesce(fastgear_cluster, "Cluster NA"),
    recent_recombination_segments = coalesce(recent_recombination_segments, 0L),
    recent_recombination_bp = coalesce(recent_recombination_bp, 0),
    recent_recombination_high_conf_segments = coalesce(recent_recombination_high_conf_segments, 0L)
  )

# Rotate only branch display order. Branch lengths and topology remain from RAxML.
# FastGEAR lineage is the primary grouping and FastGEAR cluster is secondary.
desired_order <- annotation %>%
  arrange(fastgear_lineage, fastgear_cluster, subtype_label, country, complex_order, Sample) %>%
  pull(Sample)
tree <- ape::rotateConstr(ladderize(tree, right = FALSE), desired_order)

# ---- Truncate the divergent Lineage 1 outgroup stem for display ----
# The Gubbins tree is scaled in (recombination-corrected) substitutions and the
# Lineage 1 outgroup hangs off a single enormous stem branch that would compress
# the entire ingroup into an unreadable sliver. Shorten only that one branch for
# display (drawn with an axis break) and keep every other branch length intact.
ingroup_tips <- setdiff(tree$tip.label, lineage1_outgroup)
outgroup_stem_truncated <- FALSE
stem_node <- NA_integer_
stem_true_length <- NA_real_
if (length(ingroup_tips) > 1 && length(lineage1_outgroup) >= 1) {
  stem_node <- ape::getMRCA(tree, ingroup_tips)
  stem_edge <- which(tree$edge[, 2] == stem_node)
  if (length(stem_edge) == 1) {
    stem_true_length <- tree$edge.length[stem_edge]
    tmp <- tree
    tmp$edge.length[stem_edge] <- 0
    span_without_stem <- max(node.depth.edgelength(tmp)[seq_len(ape::Ntip(tmp))])
    display_stem <- span_without_stem * 0.13
    if (is.finite(display_stem) && display_stem > 0 && stem_true_length > display_stem) {
      tree$edge.length[stem_edge] <- display_stem
      outgroup_stem_truncated <- TRUE
    }
  }
}
subtitle_text <- if (outgroup_stem_truncated) {
  sprintf("Gubbins recombination-corrected phylogeny  ·  Lineage 1 outgroup branch truncated (//, ~%.1fk SNPs)",
          stem_true_length / 1000)
} else {
  "Gubbins recombination-corrected phylogeny; branch lengths preserved"
}

annotation <- tibble(Sample = tree$tip.label) %>%
  left_join(annotation, by = "Sample")

lineage_colours <- c(
  "Lineage 1" = "#0ea5e9",
  "Lineage 2" = "#b11242",
  "Lineage NA" = "#b8c0cc"
)

cluster_colours <- c(
  "Cluster 1" = "#38bdf8",
  "Cluster 2" = "#fbbf24",
  "Cluster 3" = "#f97316",
  "Cluster 4" = "#b11242",
  "Cluster NA" = "#b8c0cc"
)

base_tree <- ggtree(tree, linewidth = 0.28, colour = branch_col)
tree_data <- base_tree$data

tip_data <- tree_data %>%
  filter(.data$isTip) %>%
  transmute(Sample = .data$label, tip_x = .data$x, y = .data$y) %>%
  left_join(annotation, by = "Sample") %>%
  mutate(
    lineage_colour = unname(lineage_colours[fastgear_lineage]),
    cluster_colour = unname(cluster_colours[fastgear_cluster])
  )

fastgear_tip_data <- tip_data %>%
  filter(fastgear_lineage != "Lineage NA", fastgear_cluster != "Cluster NA")

make_contiguous_blocks <- function(data, field, colour_field) {
  data %>%
    arrange(y) %>%
    mutate(
      value = .data[[field]],
      colour = .data[[colour_field]],
      run_id = cumsum(value != lag(value, default = first(value)))
    ) %>%
    group_by(run_id, value) %>%
    summarise(
      ymin = min(y) - 0.48,
      ymax = max(y) + 0.48,
      y = mean(range(y)),
      n = n(),
      colour = first(colour),
      samples = list(Sample),
      .groups = "drop"
    ) %>%
    filter(!str_detect(value, " NA$"))
}

lineage_blocks <- make_contiguous_blocks(fastgear_tip_data, "fastgear_lineage", "lineage_colour")
cluster_blocks <- make_contiguous_blocks(fastgear_tip_data, "fastgear_cluster", "cluster_colour")

max_x <- max(tree_data$x, na.rm = TRUE)
tree_width <- max_x
recomb_bar_x <- max_x + tree_width * 0.030
recomb_bar_max_width <- tree_width * 0.115
cluster_track_xmin <- recomb_bar_x + recomb_bar_max_width + tree_width * 0.025
cluster_track_xmax <- cluster_track_xmin + tree_width * 0.235
label_x <- cluster_track_xmin + tree_width * 0.045
plot_right <- max_x + tree_width * 0.42
plot_left <- -tree_width * 0.055
block_right <- plot_right + tree_width * 0.045
top_y <- max(tip_data$y, na.rm = TRUE)
legend_x <- plot_right + tree_width * 0.195
legend_swatch_x <- legend_x + tree_width * 0.018
legend_text_x <- legend_x + tree_width * 0.046

# Axis-break marker placed at the midpoint of the truncated outgroup stem.
break_x <- NA_real_
break_y <- NA_real_
if (isTRUE(outgroup_stem_truncated) && !is.na(stem_node)) {
  stem_child_x <- tree_data$x[match(stem_node, tree_data$node)]
  break_y <- tree_data$y[match(stem_node, tree_data$node)]
  stem_parent <- tree_data$parent[match(stem_node, tree_data$node)]
  stem_parent_x <- tree_data$x[match(stem_parent, tree_data$node)]
  break_x <- (stem_parent_x + stem_child_x) / 2
}

st492_country_points <- tip_data %>%
  filter(is_st492) %>%
  mutate(marker_x = tip_x + tree_width * 0.018)

max_recomb_log <- max(log1p(tip_data$recent_recombination_bp), na.rm = TRUE)
recombination_bars <- tip_data %>%
  filter(recent_recombination_bp > 0) %>%
  mutate(
    recombination_log = log1p(recent_recombination_bp),
    recombination_bar_width = if (max_recomb_log > 0) {
      recombination_log / max_recomb_log * recomb_bar_max_width
    } else {
      0
    },
    recombination_bar_x = recomb_bar_x,
    recombination_bar_xend = recomb_bar_x + recombination_bar_width
  )

cluster_bands <- cluster_blocks %>%
  mutate(
    xmin = vapply(samples, function(sample_names) {
      sample_names <- intersect(sample_names, tree$tip.label)
      if (length(sample_names) <= 1) {
        return(min(tip_data$tip_x[tip_data$Sample %in% sample_names], na.rm = TRUE))
      }
      node <- ape::getMRCA(tree, sample_names)
      if (is.null(node) || is.na(node)) {
        return(min(tip_data$tip_x[tip_data$Sample %in% sample_names], na.rm = TRUE))
      }
      tree_data$x[match(node, tree_data$node)]
    }, numeric(1)) - tree_width * 0.006,
    xmax = cluster_track_xmax
  )

cluster_labels <- cluster_blocks %>%
  mutate(
    label = paste0(value, "  n=", n),
    x = label_x
  )

cluster_legend <- fastgear_tip_data %>%
  count(fastgear_cluster, name = "n") %>%
  arrange(fastgear_cluster) %>%
  mutate(
    label = paste0(fastgear_cluster, " (", n, ")"),
    colour = unname(cluster_colours[fastgear_cluster]),
    y = top_y - 4.6 - row_number() * 3.25
  )

figure <- base_tree +
  geom_rect(
    data = cluster_bands,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = colour),
    inherit.aes = FALSE,
    alpha = 0.76,
    colour = "white",
    linewidth = 0.18
  ) +
  geom_tree(linewidth = 0.30, colour = branch_col) +
  geom_point(
    data = st492_country_points,
    aes(x = marker_x, y = y),
    inherit.aes = FALSE,
    shape = 16,
    size = 1.05,
    colour = "#111827",
    alpha = 0.95
  ) +
  geom_segment(
    data = recombination_bars,
    aes(x = recombination_bar_x, xend = recombination_bar_xend, y = y, yend = y),
    inherit.aes = FALSE,
    colour = "#111827",
    linewidth = 0.28,
    alpha = 0.42,
    lineend = "round"
  ) +
  geom_text(
    data = cluster_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    family = "DejaVu Sans",
    fontface = "bold",
    colour = "#111827",
    size = 2.55,
    hjust = 0
  ) +
  annotate(
    "text",
    x = max_x * 0.08,
    y = top_y + 7,
    label = "FastGEAR lineage and cluster structure in ST2/ST492",
    family = "DejaVu Sans",
    fontface = "bold",
    size = 4.5,
    colour = ink,
    hjust = 0
  ) +
  annotate(
    "text",
    x = max_x * 0.08,
    y = top_y + 3.6,
    label = subtitle_text,
    family = "DejaVu Sans",
    fontface = "bold",
    size = 2.65,
    colour = "#627286",
    hjust = 0
  ) +
  annotate("text", x = legend_x, y = top_y - 2.2, label = "FastGEAR assignment", family = "DejaVu Sans", fontface = "bold", size = 3.25, colour = ink, hjust = 0) +
  annotate("text", x = legend_x, y = top_y - 5.0, label = "Cluster", family = "DejaVu Sans", fontface = "bold", size = 2.45, colour = "#627286", hjust = 0) +
  geom_tile(data = cluster_legend, aes(x = legend_swatch_x, y = y, fill = colour), inherit.aes = FALSE, width = tree_width * 0.026, height = 1.75, colour = "white", linewidth = 0.12, alpha = 0.76) +
  geom_text(data = cluster_legend, aes(x = legend_text_x, y = y, label = label), inherit.aes = FALSE, family = "DejaVu Sans", fontface = "bold", size = 2.35, colour = ink, hjust = 0) +
  annotate("text", x = legend_x, y = top_y - 22.0, label = "Recent recombination", family = "DejaVu Sans", fontface = "bold", size = 2.45, colour = "#627286", hjust = 0) +
  annotate("segment", x = legend_swatch_x - tree_width * 0.013, xend = legend_swatch_x + tree_width * 0.013, y = top_y - 25.0, yend = top_y - 25.0, colour = "#111827", linewidth = 0.60, alpha = 0.55, lineend = "round") +
  annotate("text", x = legend_text_x, y = top_y - 25.0, label = "total bp, log-scaled", family = "DejaVu Sans", fontface = "bold", size = 2.25, colour = ink, hjust = 0) +
  annotate("text", x = legend_x, y = top_y - 31.0, label = "Sequence type", family = "DejaVu Sans", fontface = "bold", size = 2.45, colour = "#627286", hjust = 0) +
  annotate("point", x = legend_swatch_x, y = top_y - 34.0, shape = 16, size = 2.1, colour = "#111827") +
  annotate("text", x = legend_text_x, y = top_y - 34.0, label = "ST492 isolate", family = "DejaVu Sans", fontface = "bold", size = 2.35, colour = ink, hjust = 0) +
  scale_fill_identity() +
  scale_colour_identity() +
  coord_cartesian(xlim = c(plot_left, plot_right + tree_width * 0.52), ylim = c(0, top_y + 10), clip = "off") +
  theme_tree(base_family = "DejaVu Sans") +
  theme(
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white", colour = "white"),
    legend.position = "none",
    plot.margin = margin(8, 14, 8, 8)
  )

# Draw the axis break ("//") across the truncated outgroup stem: erase a slice of
# the branch, then two parallel diagonal ticks to signal the branch is not to scale.
if (isTRUE(outgroup_stem_truncated) && is.finite(break_x)) {
  brk_gap <- tree_width * 0.014
  brk_ox <- tree_width * 0.0075
  brk_sx <- tree_width * 0.0055
  brk_hy <- 2.0
  figure <- figure +
    annotate("segment", x = break_x - brk_gap, xend = break_x + brk_gap,
             y = break_y, yend = break_y, colour = "white", linewidth = 1.1) +
    annotate("segment", x = break_x - brk_ox - brk_sx, xend = break_x - brk_ox + brk_sx,
             y = break_y - brk_hy, yend = break_y + brk_hy,
             colour = branch_col, linewidth = 0.5, lineend = "round") +
    annotate("segment", x = break_x + brk_ox - brk_sx, xend = break_x + brk_ox + brk_sx,
             y = break_y - brk_hy, yend = break_y + brk_hy,
             colour = branch_col, linewidth = 0.5, lineend = "round")
}

png_path <- file.path(output_dir, "st492_raxml_tree_fastgear_blocks.png")
pdf_path <- file.path(output_dir, "st492_raxml_tree_fastgear_blocks.pdf")
csv_path <- file.path(output_dir, "st492_raxml_tree_fastgear_blocks_metadata.csv")

write_csv(
  annotation %>%
    select(
      Sample,
      subtype_label,
      country,
      fastgear_lineage,
      fastgear_cluster,
      complex_order,
      recent_recombination_segments,
      recent_recombination_bp,
      recent_recombination_high_conf_segments
    ),
  csv_path
)

png(png_path, width = 10.8, height = 12.0, units = "in", res = 450, type = "cairo", bg = "white")
print(figure)
dev.off()

cairo_pdf(pdf_path, width = 10.8, height = 12.0, bg = "white")
print(figure)
dev.off()

message("Saved: ", png_path)
message("Saved: ", pdf_path)
message("Saved metadata: ", csv_path)
message("Tree tips plotted: ", length(tree$tip.label))
message("FastGEAR tips highlighted: ", nrow(fastgear_tip_data))
message("Reference tip dropped before plotting: ", reference_dropped)
message("Rooted on FastGEAR Lineage 1 split: ", rooted_on_lineage1)
message("Outgroup stem truncated for display: ", outgroup_stem_truncated,
        if (outgroup_stem_truncated) sprintf(" (true length ~%.0f substitutions)", stem_true_length) else "")
