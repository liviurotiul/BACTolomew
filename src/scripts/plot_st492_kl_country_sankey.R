#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
})

metadata_path <- "/Data/ST492/metadata/metadata_ALL_new_merged_with_clusters_with_complex.csv"
output_dir <- "/Data/ST492/figures2026"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#132033"
region_order <- c("USA", "Western Europe", "Eastern Europe", "Levant and Arabian Peninsula")
region_colours <- c(
  "USA" = "#38bdf8",
  "Western Europe" = "#fbbf24",
  "Eastern Europe" = "#b11242",
  "Levant and Arabian Peninsula" = "#f97316"
)

normalise_country <- function(x) {
  y <- x %>% as.character() %>% str_squish() %>% na_if("")
  case_when(
    is.na(y) | str_to_lower(y) %in% c("unknown", "na", "n/a", "not available") ~ NA_character_,
    y %in% c("United Kingdom", "Scotland", "UK", "UK (Scotland)") ~ "UK (Scotland)",
    y %in% c("United States", "United States of America", "U.S.A.", "US") ~ "USA",
    TRUE ~ y
  )
}

assign_region <- function(country) {
  case_when(
    country == "USA" ~ "USA",
    country %in% c("France", "Germany", "Switzerland", "UK (Scotland)") ~ "Western Europe",
    country %in% c("Croatia", "Hungary", "Kosovo", "Montenegro", "Romania", "Serbia", "Slovenia") ~ "Eastern Europe",
    country %in% c("Lebanon", "Saudi Arabia") ~ "Levant and Arabian Peninsula",
    TRUE ~ NA_character_
  )
}

stack_nodes <- function(nodes, group_col = NULL, item_gap = 0.006, group_gap = 0.028) {
  n_items <- nrow(nodes)
  if (n_items == 0) {
    return(nodes)
  }
  group_breaks <- rep(FALSE, n_items)
  if (!is.null(group_col)) {
    groups <- nodes[[group_col]]
    group_breaks <- c(FALSE, groups[-1] != groups[-length(groups)])
  }
  gap_total <- item_gap * (n_items - 1) + group_gap * sum(group_breaks)
  available <- 1 - gap_total
  nodes <- nodes %>%
    mutate(height = total / sum(total) * available)
  y_top <- 1
  ymin <- numeric(n_items)
  ymax <- numeric(n_items)
  for (i in seq_len(n_items)) {
    if (i > 1) {
      y_top <- y_top - item_gap
    }
    if (group_breaks[i]) {
      y_top <- y_top - group_gap
    }
    ymax[i] <- y_top
    ymin[i] <- y_top - nodes$height[i]
    y_top <- ymin[i]
  }
  nodes %>%
    mutate(ymin = ymin, ymax = ymax, y = (ymin + ymax) / 2)
}

bezier <- function(p0, p1, p2, p3, t) {
  ((1 - t)^3 * p0) + (3 * (1 - t)^2 * t * p1) + (3 * (1 - t) * t^2 * p2) + (t^3 * p3)
}

make_ribbon <- function(row, x0 = 0.26, x1 = 0.74, steps = 36) {
  t <- seq(0, 1, length.out = steps)
  cx0 <- x0 + 0.22
  cx1 <- x1 - 0.22
  top_x <- bezier(x0, cx0, cx1, x1, t)
  top_y <- bezier(row$source_ymax, row$source_ymax, row$target_ymax, row$target_ymax, t)
  bottom_x <- rev(bezier(x0, cx0, cx1, x1, t))
  bottom_y <- rev(bezier(row$source_ymin, row$source_ymin, row$target_ymin, row$target_ymin, t))
  tibble(
    flow_id = row$flow_id,
    region_group = row$region_group,
    country = row$country,
    k_locus = row$k_locus,
    n = row$n,
    x = c(top_x, bottom_x),
    y = c(top_y, bottom_y)
  )
}

metadata <- read_csv(
  metadata_path,
  col_types = cols(
    .default = col_guess(),
    Sample = col_character(),
    country = col_character(),
    k_locus = col_character()
  )
) %>%
  mutate(
    subtype = suppressWarnings(as.integer(subtype)),
    country = normalise_country(country),
    k_locus = k_locus %>% str_squish() %>% na_if(""),
    region_group = assign_region(country)
  )

st492 <- metadata %>%
  filter(subtype == 492L, !is.na(country), !is.na(k_locus), !is.na(region_group))

if (nrow(st492) == 0) {
  stop("No ST492 isolates with country and K-locus metadata were found.")
}

flows <- st492 %>%
  count(region_group, country, k_locus, name = "n") %>%
  mutate(region_group = factor(region_group, levels = region_order))

country_nodes <- flows %>%
  group_by(region_group, country) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(region_group, desc(total), country) %>%
  stack_nodes(group_col = "region_group", item_gap = 0.006, group_gap = 0.030)

kl_nodes <- flows %>%
  group_by(k_locus) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total), k_locus) %>%
  stack_nodes(item_gap = 0.010, group_gap = 0)

source_parts <- flows %>%
  left_join(country_nodes %>% select(country, country_total = total, country_ymin = ymin, country_ymax = ymax), by = "country") %>%
  arrange(country, desc(n), k_locus) %>%
  group_by(country) %>%
  mutate(
    source_prev = lag(cumsum(n), default = 0),
    source_cum = cumsum(n),
    source_ymin = country_ymin + (source_prev / country_total) * (country_ymax - country_ymin),
    source_ymax = country_ymin + (source_cum / country_total) * (country_ymax - country_ymin)
  ) %>%
  ungroup()

target_parts <- flows %>%
  left_join(kl_nodes %>% select(k_locus, kl_total = total, kl_ymin = ymin, kl_ymax = ymax), by = "k_locus") %>%
  arrange(k_locus, region_group, country) %>%
  group_by(k_locus) %>%
  mutate(
    target_prev = lag(cumsum(n), default = 0),
    target_cum = cumsum(n),
    target_ymin = kl_ymin + (target_prev / kl_total) * (kl_ymax - kl_ymin),
    target_ymax = kl_ymin + (target_cum / kl_total) * (kl_ymax - kl_ymin)
  ) %>%
  ungroup()

ribbon_data <- source_parts %>%
  select(region_group, country, k_locus, n, source_ymin, source_ymax) %>%
  left_join(target_parts %>% select(country, k_locus, target_ymin, target_ymax), by = c("country", "k_locus")) %>%
  mutate(flow_id = row_number())

ribbons <- bind_rows(lapply(seq_len(nrow(ribbon_data)), function(i) make_ribbon(ribbon_data[i, ])))

region_spans <- country_nodes %>%
  group_by(region_group) %>%
  summarise(ymin = min(ymin), ymax = max(ymax), y = mean(range(c(ymin, ymax))), .groups = "drop")

figure <- ggplot() +
  geom_polygon(
    data = ribbons,
    aes(x = x, y = y, group = flow_id, fill = region_group),
    alpha = 0.68,
    colour = NA
  ) +
  geom_rect(
    data = country_nodes,
    aes(xmin = 0.05, xmax = 0.19, ymin = ymin, ymax = ymax, fill = region_group),
    colour = "white",
    linewidth = 0.34
  ) +
  geom_rect(
    data = kl_nodes,
    aes(xmin = 0.81, xmax = 0.95, ymin = ymin, ymax = ymax),
    fill = "#e8edf3",
    colour = "white",
    linewidth = 0.34
  ) +
  geom_text(
    data = country_nodes,
    aes(x = 0.03, y = y, label = paste0(country, " (", total, ")")),
    family = "DejaVu Sans",
    fontface = "bold",
    colour = ink,
    size = 2.25,
    hjust = 1
  ) +
  geom_text(
    data = kl_nodes,
    aes(x = 0.97, y = y, label = paste0(k_locus, " (", total, ")")),
    family = "DejaVu Sans",
    fontface = "bold",
    colour = ink,
    size = 2.45,
    hjust = 0
  ) +
  annotate("text", x = 0.12, y = 1.055, label = "Countries", family = "DejaVu Sans", fontface = "bold", colour = ink, size = 3.0) +
  annotate("text", x = 0.88, y = 1.055, label = "K loci", family = "DejaVu Sans", fontface = "bold", colour = ink, size = 3.0) +
  scale_fill_manual(values = region_colours, name = "Region") +
  scale_colour_manual(values = region_colours, guide = "none") +
  coord_cartesian(xlim = c(-0.18, 1.20), ylim = c(-0.02, 1.08), clip = "off") +
  labs(title = "Geographic distribution of K-loci in A. baumannii ST492") +
  theme_void(base_family = "DejaVu Sans") +
  theme(
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.title = element_text(face = "bold", size = 13.5, colour = ink, hjust = 0.5, margin = margin(b = 8)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 8.8, colour = ink),
    legend.text = element_text(face = "bold", size = 7.8, colour = ink),
    plot.margin = margin(8, 18, 8, 18)
  )

png_path <- file.path(output_dir, "st492_kl_country_sankey.png")
pdf_path <- file.path(output_dir, "st492_kl_country_sankey.pdf")

png(png_path, width = 7.8, height = 6.2, units = "in", res = 450, type = "cairo", bg = "white")
print(figure)
dev.off()

cairo_pdf(pdf_path, width = 7.8, height = 6.2, bg = "white")
print(figure)
dev.off()

message("Saved: ", png_path)
message("Saved: ", pdf_path)
message("ST492 isolates plotted: ", nrow(st492))
