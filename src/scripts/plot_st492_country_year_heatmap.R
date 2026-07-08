#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
})

metadata_path <- "/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.csv"
mlst_path <- "/home/bals/Documents/dev/BACTolomew/runs/sample_78eaa9126850/run_78eaa9126850/mlst/results.csv"
output_dir <- "/Data/ST492/figures2026"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

display_country <- function(x) {
  y <- x %>% str_squish() %>% na_if("")
  case_when(
    y %in% c("United States", "United States of America", "U.S.A.", "US") ~ "USA",
    y == "Scotland" ~ "UK (Scotland)",
    TRUE ~ y
  )
}

metadata <- read_csv(
  metadata_path,
  col_types = cols(
    .default = col_guess(),
    Sample = col_character(),
    country = col_character(),
    collection_date = col_character()
  )
) %>%
  mutate(
    Sample = as.character(Sample),
    subtype = suppressWarnings(as.integer(subtype)),
    country = display_country(country),
    collection_year = suppressWarnings(as.integer(str_extract(collection_date, "\\d{4}")))
  )

mlst <- read_csv(mlst_path, show_col_types = FALSE) %>%
  transmute(
    Sample = as.character(sample),
    mlst_sequence_type = suppressWarnings(as.integer(sequence_type))
  )

st492_all <- metadata %>%
  left_join(mlst, by = "Sample") %>%
  filter(subtype == 492L | mlst_sequence_type == 492L)

excluded_missing_country <- st492_all %>% filter(is.na(country)) %>% nrow()
excluded_missing_year <- st492_all %>% filter(is.na(collection_year)) %>% nrow()

st492 <- st492_all %>%
  filter(!is.na(country), !is.na(collection_year))

if (nrow(st492) == 0) {
  stop("No ST492 isolates with both country and collection-year metadata were found.")
}

country_totals <- st492 %>%
  count(country, name = "total") %>%
  arrange(total, country)

year_totals <- st492 %>%
  count(collection_year, name = "total")

year_levels <- seq(min(st492$collection_year), max(st492$collection_year))
country_levels <- country_totals$country

country_year_counts <- st492 %>%
  count(country, collection_year, name = "n") %>%
  complete(
    country = country_levels,
    collection_year = year_levels,
    fill = list(n = 0L)
  ) %>%
  mutate(
    country = factor(country, levels = country_levels),
    year_label = factor(as.character(collection_year), levels = as.character(year_levels)),
    count_label = if_else(n > 0L, as.character(n), "")
  )

year_counts <- tibble(collection_year = year_levels) %>%
  left_join(year_totals, by = "collection_year") %>%
  mutate(
    total = replace_na(total, 0L),
    year_label = factor(as.character(collection_year), levels = as.character(year_levels)),
    total_label = if_else(total > 0L, as.character(total), "")
  )

compact_country_label <- function(x) {
  case_when(
    x == "UK (Scotland)" ~ "UK\n(Scotland)",
    x == "Saudi Arabia" ~ "Saudi\nArabia",
    TRUE ~ x
  )
}

pal_low <- "#fee8c8"
pal_mid <- "#f97316"
pal_high <- "#9f1239"
zero_fill <- "#f8fafc"
ink <- "#132033"
muted <- "#64748b"
grid_line <- "#e7edf4"

positive_counts <- country_year_counts$n[country_year_counts$n > 0]
mid_count <- if (length(positive_counts) > 0) median(positive_counts) else 1
max_count <- max(country_year_counts$n)

count_scale <- scale_fill_gradientn(
  colours = c(zero_fill, pal_low, pal_mid, pal_high),
  values = rescale(c(0, 1, mid_count, max_count)),
  limits = c(0, max_count),
  breaks = pretty_breaks(n = 4),
  labels = label_number(accuracy = 1),
  name = "Isolates"
)

base_theme <- theme_minimal(base_family = "DejaVu Sans") +
  theme(
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_text(face = "bold", size = 7.2, colour = ink),
    legend.text = element_text(size = 6.4, colour = muted),
    axis.title = element_blank(),
    axis.text = element_text(face = "bold", colour = ink),
    plot.title = element_text(face = "bold", size = 11, colour = ink, hjust = 0),
    plot.subtitle = element_text(size = 6.8, colour = muted, hjust = 0, margin = margin(t = 1, b = 4)),
    plot.margin = margin(2, 4, 2, 4)
  )

top_totals_plot <- ggplot(year_counts, aes(year_label, total)) +
  geom_col(width = 0.82, fill = pal_mid, colour = "white", linewidth = 0.18) +
  scale_x_discrete(drop = FALSE, expand = expansion(add = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "Yearly totals") +
  base_theme +
  theme(
    plot.title = element_text(size = 7.8, margin = margin(b = 1)),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6.2, colour = muted),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = grid_line, linewidth = 0.18),
    plot.margin = margin(0, 6, 0, 31)
  )

heatmap_plot <- ggplot(country_year_counts, aes(year_label, country, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.35, width = 0.96, height = 0.92) +
  geom_text(
    aes(label = count_label),
    family = "DejaVu Sans",
    fontface = "bold",
    size = 2.15,
    colour = ink
  ) +
  count_scale +
  scale_x_discrete(drop = FALSE, expand = expansion(add = 0)) +
  scale_y_discrete(drop = FALSE, expand = expansion(add = 0), labels = compact_country_label) +
  labs(title = NULL, subtitle = NULL) +
  base_theme +
  theme(
    axis.text.x = element_text(size = 6.35, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 6.35, margin = margin(r = 1)),
    legend.key.height = unit(0.22, "in"),
    legend.key.width = unit(0.12, "in"),
    plot.margin = margin(1, 0, 2, 0)
  )

figure <- heatmap_plot
figure <- figure +
  plot_annotation(
    title = "A. baumannii ST492 over time and geography",
    theme = theme(
      plot.background = element_rect(fill = "white", colour = "white"),
      plot.title.position = "plot",
      plot.title = element_text(
        family = "DejaVu Sans",
        face = "bold",
        size = 10.5,
        colour = ink,
        hjust = 0.5,
        margin = margin(b = 5)
      ),
      plot.margin = margin(6, 4, 5, 0)
    )
  )

png_path <- file.path(output_dir, "st492_country_year_heatmap.png")
pdf_path <- file.path(output_dir, "st492_country_year_heatmap.pdf")
csv_path <- file.path(output_dir, "st492_country_year_counts.csv")

figure_width_in <- 6.18
figure_height_in <- 3.85

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
  country_year_counts %>%
    transmute(
      country = as.character(country),
      collection_year,
      st492_isolates = n
    ) %>%
    filter(st492_isolates > 0) %>%
    arrange(collection_year, country),
  csv_path
)

message("Saved: ", png_path)
message("Saved: ", pdf_path)
message("Saved: ", csv_path)
message("ST492 isolates plotted: ", nrow(st492))
message("Excluded for missing country: ", excluded_missing_country)
message("Excluded for missing collection year: ", excluded_missing_year)
