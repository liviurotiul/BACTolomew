#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
})

metadata_path <- "/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.pre_standardize_2026-06-11.csv"
mlst_path <- "/home/bals/Documents/dev/BACTolomew/runs/sample_78eaa9126850/run_78eaa9126850/mlst/results.csv"
output_dir <- "/Data/ST492/figures2026"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

normalise_broad_source <- function(x) {
  y <- x %>% str_squish() %>% str_to_lower() %>% na_if("")
  case_when(
    is.na(y) | y %in% c("unknown", "not collected", "not available", "na", "n/a") ~ "Not available",
    y == "clinical" ~ "Clinical",
    y == "animal" ~ "Animal",
    y == "wastewater" ~ "Wastewater",
    y == "respiratory" ~ "Respiratory",
    y == "sterile body site" ~ "Sterile site/body fluid",
    TRUE ~ str_to_sentence(y)
  )
}

normalise_detail_source <- function(x) {
  y <- x %>% str_squish() %>% str_to_lower() %>% na_if("")
  case_when(
    is.na(y) | y %in% c("unknown", "not collected", "not available", "na", "n/a") ~ "Not available",
    y == "clinical sample" ~ "Not available/generic",
    y == "wound/skin" ~ "Wound/skin",
    y == "sterile site/body fluid" ~ "Sterile site/body fluid",
    TRUE ~ str_to_sentence(y)
  )
}

make_wedge <- function(row, inner_radius = 0.44, outer_radius = 1.0, steps = 120) {
  theta_outer <- seq(row$theta_start, row$theta_end, length.out = steps)
  theta_inner <- rev(theta_outer)
  offset_x <- cos(row$theta_mid) * row$explode
  offset_y <- sin(row$theta_mid) * row$explode
  tibble(
    segment_id = row$segment_id,
    source = row$source,
    n = row$n,
    x = c(cos(theta_outer) * outer_radius, cos(theta_inner) * inner_radius) + offset_x,
    y = c(sin(theta_outer) * outer_radius, sin(theta_inner) * inner_radius) + offset_y
  )
}

adjust_label_positions <- function(label_df, min_gap = 0.15, ymin = -1.03, ymax = 1.03) {
  adjust_one_side <- function(d) {
    if (nrow(d) == 0) {
      return(d)
    }
    d <- d[order(d$label_y_raw), ]
    y <- d$label_y_raw
    if (length(y) > 1) {
      for (i in 2:length(y)) {
        y[i] <- max(y[i], y[i - 1] + min_gap)
      }
    }
    overflow <- max(y) - ymax
    if (overflow > 0) {
      y <- y - overflow
    }
    underflow <- ymin - min(y)
    if (underflow > 0) {
      y <- y + underflow
    }
    d$label_y <- y
    d
  }

  bind_rows(
    adjust_one_side(label_df %>% filter(side == "left")),
    adjust_one_side(label_df %>% filter(side == "right"))
  )
}

metadata <- read_csv(
  metadata_path,
  col_types = cols(
    .default = col_guess(),
    Sample = col_character(),
    isolation_source = col_character(),
    isolation_source_detail = col_character()
  )
) %>%
  mutate(
    Sample = as.character(Sample),
    subtype = suppressWarnings(as.integer(subtype))
  )

mlst <- read_csv(mlst_path, show_col_types = FALSE) %>%
  transmute(
    Sample = as.character(sample),
    mlst_sequence_type = suppressWarnings(as.integer(sequence_type))
  )

st492 <- metadata %>%
  left_join(mlst, by = "Sample") %>%
  filter(subtype == 492L | mlst_sequence_type == 492L) %>%
  mutate(
    broad_source = normalise_broad_source(isolation_source),
    detail_source = normalise_detail_source(isolation_source_detail)
  )

if (nrow(st492) == 0) {
  stop("No ST492 isolates were found.")
}

total_n <- nrow(st492)

ink <- "#132033"
muted <- "#64748b"
line_colour <- "#64748b"
broad_order <- c("Clinical", "Not available", "Animal", "Wastewater", "Respiratory", "Sterile site/body fluid")
segment_order <- c(
  "Clinical: Blood",
  "Clinical: Respiratory",
  "Clinical: Wound/skin",
  "Clinical: Fecal",
  "Clinical: Not available/generic",
  "Clinical: Urine",
  "Not available",
  "Animal",
  "Wastewater"
)
segment_colours <- c(
  "Clinical: Blood" = "#b11242",
  "Clinical: Respiratory" = "#b11242",
  "Clinical: Wound/skin" = "#b11242",
  "Clinical: Fecal" = "#b11242",
  "Clinical: Not available/generic" = "#b11242",
  "Clinical: Urine" = "#b11242",
  "Not available" = "#f97316",
  "Animal" = "#fbbf24",
  "Wastewater" = "#38bdf8"
)

clinical_segments <- st492 %>%
  filter(broad_source == "Clinical") %>%
  count(detail_source, name = "n") %>%
  mutate(
    broad_source = "Clinical",
    source = paste0("Clinical: ", detail_source),
    label_name = detail_source
  )

other_segments <- st492 %>%
  filter(broad_source != "Clinical") %>%
  count(broad_source, name = "n") %>%
  transmute(
    broad_source = broad_source,
    source = broad_source,
    label_name = broad_source,
    n
  )

segments <- bind_rows(
  clinical_segments %>% select(broad_source, source, label_name, n),
  other_segments
) %>%
  mutate(
    source = factor(source, levels = segment_order),
    fraction = n / total_n
  ) %>%
  arrange(source) %>%
  mutate(
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    pad = if_else(str_starts(as.character(source), "Clinical:"), 0.006, 0.014),
    theta_start_raw = pi / 2 - 2 * pi * ymin,
    theta_end_raw = pi / 2 - 2 * pi * ymax,
    arc_direction = sign(theta_end_raw - theta_start_raw),
    theta_start = theta_start_raw + arc_direction * pad,
    theta_end = theta_end_raw - arc_direction * pad,
    theta_mid = (theta_start + theta_end) / 2,
    explode = 0.0,
    segment_id = paste0("segment_", row_number()),
    percent_label = percent(n / total_n, accuracy = 0.1),
    label = paste0(label_name, "\n", n, " | ", percent_label)
  )

inner_segments <- st492 %>%
  count(broad_source, name = "n") %>%
  mutate(
    broad_source = factor(broad_source, levels = broad_order),
    source = as.character(broad_source),
    label_name = as.character(broad_source),
    fraction = n / total_n
  ) %>%
  arrange(broad_source) %>%
  mutate(
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    pad = 0.006,
    theta_start_raw = pi / 2 - 2 * pi * ymin,
    theta_end_raw = pi / 2 - 2 * pi * ymax,
    arc_direction = sign(theta_end_raw - theta_start_raw),
    theta_start = theta_start_raw + arc_direction * pad,
    theta_end = theta_end_raw - arc_direction * pad,
    theta_mid = (theta_start + theta_end) / 2,
    explode = 0.0,
    segment_id = paste0("inner_", row_number()),
    label = paste0(label_name, "\n", n)
  )

outer_wedges <- bind_rows(lapply(seq_len(nrow(segments)), function(i) {
  make_wedge(segments[i, ], inner_radius = 0.74, outer_radius = 1.18, steps = 140)
})) %>%
  mutate(ring = "outer")

inner_wedges <- bind_rows(lapply(seq_len(nrow(inner_segments)), function(i) {
  make_wedge(inner_segments[i, ], inner_radius = 0.34, outer_radius = 0.71, steps = 140)
})) %>%
  mutate(ring = "inner")

wedges <- bind_rows(inner_wedges, outer_wedges)

label_data <- segments %>%
  filter(broad_source == "Clinical") %>%
  mutate(
    display_label = case_when(
      label_name == "Not available/generic" ~ "Not\navailable",
      label_name == "Wound/skin" ~ "Wound/\nskin",
      label_name == "Respiratory" ~ "Respir-\natory",
      TRUE ~ label_name
    ),
    label = display_label,
    label_radius = case_when(
      label_name == "Wound/skin" ~ 0.970,
      label_name == "Respiratory" ~ 0.940,
      n <= 4 ~ 0.975,
      n <= 6 ~ 0.965,
      TRUE ~ 0.960
    ),
    x = cos(theta_mid) * label_radius,
    y = sin(theta_mid) * label_radius,
    raw_angle = theta_mid * 180 / pi,
    angle = if_else(raw_angle > 90 | raw_angle < -90, raw_angle + 180, raw_angle),
    text_size = case_when(
      n >= 20 ~ 2.45,
      label_name == "Respiratory" ~ 1.60,
      label_name == "Wound/skin" ~ 1.72,
      label_name == "Not available/generic" ~ 1.58,
      n >= 10 ~ 1.90,
      str_length(label_name) >= 14 ~ 1.08,
      TRUE ~ 1.46
    ),
    text_colour = "white"
  )

inner_label_data <- inner_segments %>%
  mutate(
    display_label = case_when(
      label_name == "Not available" ~ "Not\navailable",
      label_name == "Wastewater" ~ "Waste\nwater",
      TRUE ~ label_name
    ),
    label = display_label,
    label_radius = if_else(label_name == "Wastewater", 0.560, 0.525),
    x = cos(theta_mid) * label_radius,
    y = sin(theta_mid) * label_radius,
    raw_angle = theta_mid * 180 / pi,
    angle = if_else(raw_angle > 90 | raw_angle < -90, raw_angle + 180, raw_angle),
    text_size = case_when(
      label_name == "Not available" ~ 1.74,
      label_name == "Wastewater" ~ 1.08,
      n >= 20 ~ 2.18,
      TRUE ~ 1.50
    ),
    text_colour = if_else(source %in% c("Clinical", "Not available"), "white", ink)
  )

percent_label_data <- segments %>%
  mutate(
    label_radius = 1.390,
    x = cos(theta_mid) * label_radius,
    y = sin(theta_mid) * label_radius,
    raw_angle = theta_mid * 180 / pi,
    angle = if_else(raw_angle > 90 | raw_angle < -90, raw_angle + 180, raw_angle),
    label = paste0(percent(n / total_n, accuracy = 0.1), " (", n, ")"),
    text_size = case_when(
      n >= 10 ~ 1.88,
      TRUE ~ 1.55
    )
  )

center_label <- tibble(x = 0, y = 0, label = paste0("ST492\nn = ", total_n))

clinical_total <- segments %>%
  filter(broad_source == "Clinical") %>%
  summarise(n = sum(n), .groups = "drop") %>%
  pull(n)

clinical_angle_start <- segments %>%
  filter(broad_source == "Clinical") %>%
  summarise(theta = max(theta_start_raw), .groups = "drop") %>%
  pull(theta)

clinical_angle_end <- segments %>%
  filter(broad_source == "Clinical") %>%
  summarise(theta = min(theta_end_raw), .groups = "drop") %>%
  pull(theta)

clinical_bracket_radius <- 1.700
clinical_bracket <- tibble(
  theta = seq(clinical_angle_start - 0.020, clinical_angle_end + 0.020, length.out = 260),
  x = cos(theta) * clinical_bracket_radius,
  y = sin(theta) * clinical_bracket_radius
)

clinical_bracket_ticks <- tibble(
  theta = c(clinical_angle_start - 0.020, clinical_angle_end + 0.020),
  x = cos(theta) * clinical_bracket_radius,
  y = sin(theta) * clinical_bracket_radius,
  xend = cos(theta) * (clinical_bracket_radius + 0.050),
  yend = sin(theta) * (clinical_bracket_radius + 0.050)
)

clinical_bracket_label_angle <- (clinical_angle_start + clinical_angle_end) / 2
clinical_bracket_label <- tibble(
  x = cos(clinical_bracket_label_angle) * 1.890,
  y = sin(clinical_bracket_label_angle) * 1.890,
  label = paste0("Clinical\n", percent(clinical_total / total_n, accuracy = 0.1))
)

broad_colours <- c(
  "Clinical" = "#b11242",
  "Not available" = "#f97316",
  "Animal" = "#fbbf24",
  "Wastewater" = "#38bdf8"
)
fill_colours <- c(broad_colours, segment_colours)

figure <- ggplot() +
  geom_polygon(
    data = wedges,
    aes(x, y, group = segment_id, fill = source),
    colour = "white",
    linewidth = 0.55
  ) +
  geom_text(
    data = inner_label_data,
    aes(x = x, y = y, label = label, angle = angle, size = text_size, colour = text_colour),
    family = "DejaVu Sans",
    fontface = "bold",
    lineheight = 0.90,
    hjust = 0.5,
    vjust = 0.5
  ) +
  geom_text(
    data = label_data,
    aes(x = x, y = y, label = label, angle = angle, size = text_size, colour = text_colour),
    family = "DejaVu Sans",
    fontface = "bold",
    lineheight = 0.90,
    hjust = 0.5,
    vjust = 0.5
  ) +
  geom_text(
    data = percent_label_data,
    aes(x = x, y = y, label = label, angle = angle, size = text_size),
    family = "DejaVu Sans",
    fontface = "bold",
    colour = ink,
    hjust = 0.5,
    vjust = 0.5
  ) +
  geom_path(
    data = clinical_bracket,
    aes(x, y),
    colour = alpha("#7f102f", 0.82),
    linewidth = 0.34,
    lineend = "round"
  ) +
  geom_segment(
    data = clinical_bracket_ticks,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = alpha("#7f102f", 0.82),
    linewidth = 0.34,
    lineend = "round"
  ) +
  geom_label(
    data = clinical_bracket_label,
    aes(x = x, y = y, label = label),
    family = "DejaVu Sans",
    fontface = "bold",
    colour = "#7f102f",
    fill = alpha("white", 0.94),
    lineheight = 0.90,
    label.r = unit(0.07, "lines"),
    label.padding = unit(0.10, "lines"),
    linewidth = 0,
    size = 2.05
  ) +
  geom_text(
    data = center_label,
    aes(x, y, label = label),
    family = "DejaVu Sans",
    fontface = "bold",
    size = 2.65,
    lineheight = 0.92,
    colour = ink
  ) +
  scale_fill_manual(
    values = fill_colours,
    breaks = names(broad_colours),
    labels = names(broad_colours),
    guide = guide_legend(
      title = "Source",
      override.aes = list(colour = "white", linewidth = 0.35)
    )
  ) +
  scale_colour_identity() +
  scale_size_identity() +
  coord_equal(xlim = c(-1.42, 2.02), ylim = c(-1.86, 1.86), clip = "off") +
  labs(title = "Isolation source composition of A. baumannii ST492 isolates") +
  theme_void(base_family = "DejaVu Sans") +
  theme(
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.title = element_text(
      family = "DejaVu Sans",
      face = "bold",
      size = 10.5,
      colour = ink,
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    legend.position = "right",
    legend.justification = "center",
    legend.title = element_text(
      family = "DejaVu Sans",
      face = "bold",
      size = 9.2,
      colour = ink
    ),
    legend.text = element_text(
      family = "DejaVu Sans",
      face = "bold",
      size = 8.0,
      colour = ink
    ),
    legend.key.size = unit(0.36, "cm"),
    legend.spacing.y = unit(0.08, "cm"),
    legend.margin = margin(0, 0, 0, 8),
    plot.title.position = "plot",
    plot.margin = margin(8, 14, 8, 0)
  )

png_path <- file.path(output_dir, "st492_isolation_source_donut.png")
pdf_path <- file.path(output_dir, "st492_isolation_source_donut.pdf")
csv_path <- file.path(output_dir, "st492_isolation_source_counts.csv")

figure_width_in <- 7.8
figure_height_in <- 5.0

png(
  filename = png_path,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  res = 450,
  type = "cairo",
  bg = "white"
)
print(figure)
dev.off()

cairo_pdf(
  filename = pdf_path,
  width = figure_width_in,
  height = figure_height_in,
  bg = "white"
)
print(figure)
dev.off()

write_csv(
  st492 %>%
    count(broad_source, detail_source, name = "st492_isolates") %>%
    group_by(broad_source) %>%
    mutate(percent_within_source = percent(st492_isolates / sum(st492_isolates), accuracy = 0.1)) %>%
    ungroup() %>%
    transmute(
      isolation_source = broad_source,
      isolation_source_detail = detail_source,
      st492_isolates,
      percent_within_source
    ) %>%
    arrange(isolation_source, desc(st492_isolates), isolation_source_detail),
  csv_path
)

message("Saved: ", png_path)
message("Saved: ", pdf_path)
message("Saved: ", csv_path)
message("ST492 isolates plotted: ", total_n)
