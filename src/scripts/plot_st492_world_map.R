#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(maps)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
})

metadata_path <- "/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.csv"
mlst_path <- "/home/bals/Documents/dev/BACTolomew/runs/sample_78eaa9126850/run_78eaa9126850/mlst/results.csv"
output_dir <- "/Data/ST492/figures2026"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

normalise_country <- function(x) {
  y <- x %>% str_squish() %>% na_if("")
  case_when(
    y %in% c("United States", "United States of America", "U.S.A.", "US") ~ "USA",
    y %in% c("United Kingdom", "UK (Scotland)", "Scotland") ~ "UK",
    TRUE ~ y
  )
}

display_country <- function(x) {
  y <- x %>% str_squish() %>% na_if("")
  case_when(
    y %in% c("United States", "United States of America", "U.S.A.", "US") ~ "USA",
    y == "Scotland" ~ "UK (Scotland)",
    TRUE ~ y
  )
}

metadata <- read_csv(metadata_path, show_col_types = FALSE) %>%
  mutate(
    Sample = as.character(Sample),
    subtype = suppressWarnings(as.integer(subtype)),
    country = display_country(as.character(country))
  )

mlst <- read_csv(mlst_path, show_col_types = FALSE) %>%
  transmute(
    Sample = as.character(sample),
    mlst_sequence_type = suppressWarnings(as.integer(sequence_type))
  )

st492 <- metadata %>%
  left_join(mlst, by = "Sample") %>%
  filter(subtype == 492L | mlst_sequence_type == 492L) %>%
  filter(!is.na(country))

if (nrow(st492) == 0) {
  stop("No ST492 isolates with country metadata were found.")
}

country_counts <- st492 %>%
  count(country, name = "n") %>%
  mutate(
    map_region = normalise_country(country),
    label = paste0(country, "\n", n)
  ) %>%
  arrange(desc(n), country)

manual_coords <- tribble(
  ~country, ~lon, ~lat,
  "Croatia", 16.40, 45.10,
  "France", 2.20, 46.20,
  "Germany", 10.45, 51.20,
  "Hungary", 19.40, 47.15,
  "Kosovo", 20.90, 42.60,
  "Lebanon", 35.85, 33.85,
  "Montenegro", 19.30, 42.75,
  "Romania", 24.95, 45.95,
  "Saudi Arabia", 45.00, 24.00,
  "UK (Scotland)", -4.20, 56.60,
  "Serbia", 20.80, 44.05,
  "Slovenia", 14.90, 46.15,
  "Switzerland", 8.25, 46.80,
  "USA", -98.50, 39.60
)

world <- map_data("world") %>%
  as_tibble()

map_counts <- country_counts %>%
  group_by(map_region) %>%
  summarise(n = sum(n), countries = paste(country, collapse = ", "), .groups = "drop")

world_hits <- world %>%
  inner_join(map_counts, by = c("region" = "map_region"))

computed_coords <- world_hits %>%
  group_by(region) %>%
  summarise(
    lon = mean(range(long, na.rm = TRUE)),
    lat = mean(range(lat, na.rm = TRUE)),
    .groups = "drop"
  )

country_counts <- country_counts %>%
  left_join(manual_coords, by = "country") %>%
  left_join(computed_coords, by = c("map_region" = "region"), suffix = c("", "_computed")) %>%
  mutate(
    lon = coalesce(lon, lon_computed),
    lat = coalesce(lat, lat_computed),
    n_label = comma(n)
  ) %>%
  select(-lon_computed, -lat_computed)

leader_label_countries <- c(
  "Lebanon", "Kosovo", "Croatia", "Serbia", "Slovenia", "Hungary", "Switzerland",
  "Montenegro"
)

label_positions <- tribble(
  ~country, ~label_lon, ~label_lat,
  "Croatia", 15.65, 43.55,
  "Hungary", 23.20, 49.65,
  "Kosovo", 22.10, 40.65,
  "Lebanon", 39.95, 34.55,
  "Montenegro", 17.85, 41.75,
  "Serbia", 24.45, 43.45,
  "Slovenia", 15.35, 47.55,
  "Switzerland", 7.35, 45.70
)

country_counts <- country_counts %>%
  left_join(label_positions, by = "country") %>%
  mutate(
    label_lon = coalesce(label_lon, lon),
    label_lat = coalesce(label_lat, lat),
    use_leader = country %in% leader_label_countries
  )

pal_low <- "#fee8c8"
pal_mid <- "#f97316"
pal_high <- "#9f1239"
land <- "#e8edf3"
coast <- "#ffffff"
ink <- "#132033"
muted <- "#637083"

base_map_theme <- theme_void(base_family = "DejaVu Sans") +
  theme(
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.title = element_text(
      face = "bold",
      size = 8.6,
      colour = ink,
      hjust = 0,
      margin = margin(b = 1)
    ),
    legend.position = "none",
    plot.margin = margin(1.5, 2, 1.5, 2)
  )

count_fill <- scale_fill_gradientn(
  colours = c(pal_low, pal_mid, pal_high),
  values = rescale(c(1, median(country_counts$n), max(country_counts$n))),
  limits = c(1, max(country_counts$n)),
  breaks = pretty_breaks(n = 4),
  labels = label_number(accuracy = 1),
  name = "ST492 isolates"
)

map_aspect <- function(xlim, ylim) {
  mean_lat <- mean(ylim)
  diff(ylim) / (diff(xlim) * cos(mean_lat * pi / 180))
}

europe_xlim <- c(-12.5, 32.5)
europe_ylim <- c(36.2, 60.2)
north_america_xlim <- c(-126.5, -66.5)
north_america_ylim <- c(23.5, 53.5)
arabia_xlim <- c(31, 58)
arabia_ylim <- c(12.8, 37.2)

europe_width_in <- 6.8
inter_map_gap_in <- 0.10
outer_margin_width_in <- 0.16
outer_margin_height_in <- 0.10
panel_title_height_in <- 0.34
main_title_height_in <- 0.52

europe_aspect <- map_aspect(europe_xlim, europe_ylim)
north_america_aspect <- map_aspect(north_america_xlim, north_america_ylim)
arabia_aspect <- map_aspect(arabia_xlim, arabia_ylim)

lower_row_map_height_in <- (europe_width_in - inter_map_gap_in) /
  ((1 / north_america_aspect) + (1 / arabia_aspect))
north_america_width_in <- lower_row_map_height_in / north_america_aspect
arabia_width_in <- lower_row_map_height_in / arabia_aspect
europe_map_height_in <- europe_width_in * europe_aspect

europe_row_height_in <- europe_map_height_in + panel_title_height_in
lower_row_height_in <- lower_row_map_height_in + panel_title_height_in
figure_width_in <- europe_width_in + (2 * outer_margin_width_in)
figure_height_in <- main_title_height_in + europe_row_height_in + lower_row_height_in + (2 * outer_margin_height_in)

regional_map <- function(title, xlim, ylim, countries = NULL) {
  point_data <- country_counts
  if (!is.null(countries)) {
    point_data <- point_data %>% filter(country %in% countries)
  }
  leader_data <- point_data %>% filter(use_leader)
  ggplot() +
    geom_polygon(
      data = world,
      aes(long, lat, group = group),
      fill = land,
      colour = coast,
      linewidth = 0.14
    ) +
    geom_polygon(
      data = world_hits,
      aes(long, lat, group = group, fill = n),
      colour = "white",
      linewidth = 0.24
    ) +
    geom_segment(
      data = leader_data,
      aes(x = lon, y = lat, xend = label_lon, yend = label_lat),
      inherit.aes = FALSE,
      colour = alpha(ink, 0.62),
      linewidth = 0.32
    ) +
    geom_label(
      data = point_data,
      aes(label_lon, label_lat, label = label),
      inherit.aes = FALSE,
      family = "DejaVu Sans",
      fontface = "bold",
      size = 1.95,
      colour = ink,
      fill = alpha("white", 0.94),
      linewidth = 0.10,
      label.r = unit(0.07, "lines"),
      label.padding = unit(0.09, "lines")
    ) +
    count_fill +
    coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE, clip = "on") +
    labs(title = title) +
    guides(fill = "none") +
    base_map_theme
}

europe_plot <- regional_map(
  "Europe",
  xlim = europe_xlim,
  ylim = europe_ylim,
  countries = c(
    "Croatia", "France", "Germany", "Hungary", "Kosovo", "Montenegro",
    "Romania", "UK (Scotland)", "Serbia", "Slovenia", "Switzerland"
  )
)

north_america_plot <- regional_map(
  "North America",
  xlim = north_america_xlim,
  ylim = north_america_ylim,
  countries = c("USA")
)

arabia_plot <- regional_map(
  "Levant & Arabian Peninsula",
  xlim = arabia_xlim,
  ylim = arabia_ylim,
  countries = c("Saudi Arabia", "Lebanon")
)

regional_row <- north_america_plot | plot_spacer() | arabia_plot
regional_row <- regional_row +
  plot_layout(widths = unit(c(north_america_width_in, inter_map_gap_in, arabia_width_in), "in"))

figure <- europe_plot / regional_row
figure <- figure +
  plot_layout(heights = unit(c(europe_row_height_in, lower_row_height_in), "in"), guides = "collect") +
  plot_annotation(
    title = "Global distribution of Acinetobacter baumannii ST492",
    theme = theme(
      plot.title = element_text(
        family = "DejaVu Sans",
        face = "bold",
        size = 10.5,
        colour = ink,
        margin = margin(b = 2)
      ),
      plot.background = element_rect(fill = "white", colour = "white"),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

png_path <- file.path(output_dir, "st492_global_country_map.png")
pdf_path <- file.path(output_dir, "st492_global_country_map.pdf")
csv_path <- file.path(output_dir, "st492_country_counts.csv")

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

write_csv(country_counts %>% select(country, st492_isolates = n), csv_path)

message("Saved: ", png_path)
message("Saved: ", pdf_path)
message("Saved: ", csv_path)
message(
  "Computed figure geometry: Europe width=", round(europe_width_in, 2), " in; ",
  "North America width=", round(north_america_width_in, 2), " in; ",
  "gap=", round(inter_map_gap_in, 2), " in; ",
  "Arabia width=", round(arabia_width_in, 2), " in; ",
  "canvas=", round(figure_width_in, 2), " x ", round(figure_height_in, 2), " in"
)
