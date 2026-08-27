#!/usr/bin/env Rscript
### DATA ANALYSIS ###

##### Load libraries #####
library(phyloseq)
library(dplyr)
library(ggplot2)
library(scales)
library(microbiome)

##### Load data #####

args <- commandArgs(trailingOnly = TRUE)

otu_dir   <-as.character(args[1])
output_dir<-as.character(args[2])

tax_file      <- args[3]
meta_id       <- args[4]
metadata_file <- args[5]

k_value       <- args[6]
r_value       <- args[7]

#otu_file <- list.files(otu_dir, pattern = "\\krakenuniq_abundance_matrix.txt$", full.names = TRUE)

file_pattern <- paste0("^krakenuniq_abundance_matrix_k", k_value, "_r", r_value, "\\.txt$")

otu_file <- list.files(
  path=otu_dir,
  pattern=file_pattern,
  full.names = TRUE)

sample_names <- scan( otu_file, what = character(), nlines = 1, quiet = TRUE )

otu_df <- data.table::fread( otu_file, skip = 1, header = FALSE, col.names = c("taxid", sample_names), colClasses = list(character = 1) )

rownames(otu_df) <- otu_df$taxid
print(otu_df)
otu_df <- read.table( otu_file, header = TRUE, row.names = 1, sep = "",check.names = FALSE, comment.char = "", quote = "" )

otu_mat <- as.matrix(otu_df)

otu_ps <- otu_table(otu_mat, taxa_are_rows = TRUE)

tax_table_df <- data.table::fread(tax_file, header=T)
print(tax_table_df)

metadata <- data.table::fread(metadata_file, header=T)

print(metadata)
metadata <- sample_data(metadata)

otu_ids <- rownames(otu_mat)

tax_table_df <- tax_table_df[rownames(tax_table_df) %in% otu_ids, ]
tax_table_df <- tax_table_df[match(otu_ids, rownames(tax_table_df)), ]

common_ids <- intersect(otu_ids, rownames(tax_table_df))
print(common_ids)

otu_ps2 <- prune_taxa(common_ids, otu_ps)
tax_table_df2 <- tax_table_df[common_ids, ]
print(tax_table_df2)

tax_ps <- tax_table(as.matrix(tax_table_df2))
print(tax_ps)

ps <- phyloseq(otu_ps2, tax_ps, metadata)

##### Colours #####

distinct_bright_cols <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#00CED1", "#000000",
  "#1B9E77", "#FF4500", "#7570B3", "#000080", "#66A61E",
  "#E6AB02", "#A6761D", "#666666", "#00BFFF", "#FF1493"
)


pastel_palette <- colorRampPalette(c(
  "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6",
  "#FFFFCC", "#E5D8BD", "#FDDAEC", "#D9D9D9", "#CCE5FF",
  "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCFF", "#E0FFFF"
))

##### Helper function: top 20 overall plot #####

make_tax_plot_top20_overall <- function(df_input, taxon_col, filename) {
  
  df_input <- df_input %>%
    mutate(taxon_plot = as.character(.data[[taxon_col]]))
  
  top20_taxa <- df_input %>%
    group_by(taxon_plot) %>%
    summarise(total_abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_abundance)) %>%
    slice_head(n = 20) %>%
    pull(taxon_plot)
  
  other_taxa <- setdiff(unique(df_input$taxon_plot), top20_taxa)
  
  bright_cols <- distinct_bright_cols[1:length(top20_taxa)]
  names(bright_cols) <- top20_taxa
  
  pastel_cols <- pastel_palette(length(other_taxa))
  names(pastel_cols) <- other_taxa
  
  taxon_cols <- c(bright_cols, pastel_cols)
  
  depth_order <- df_input %>%
    distinct(depth) %>%
    arrange(desc(depth)) %>%
    pull(depth)
  
  df_input <- df_input %>%
    mutate(
      taxon_plot = factor(taxon_plot, levels = c(top20_taxa, other_taxa)),
      depth_plot = factor(depth, levels = depth_order)
    )
  
  read_labels <- df_input %>%
    distinct(depth_plot, total_reads, read_label)
  
  p <- ggplot(df_input, aes(x = depth_plot, y = Abundance, fill = taxon_plot)) +
    geom_col(width = 0.8) +
    geom_text(
      data = read_labels,
      aes(x = depth_plot, y = 1.03, label = read_label),
      inherit.aes = FALSE,
      size = 6,
      fontface = "bold",
      hjust = 0
    ) +
    labs(x = "Depth (cm)", y = "Relative abundance") +
    scale_y_continuous(
      labels = scales::percent_format(),
      limits = c(0, 1.18),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_manual(
      values = taxon_cols,
      breaks = top20_taxa
    ) +
    coord_flip(clip = "off") +
    theme(
      axis.text.y = element_text(size = 20, face = "bold"),
      axis.text.x = element_text(size = 20, face = "bold"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 20, face = "bold"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 20, face = "bold"),
      legend.position = "bottom",
      panel.background = element_blank(),
      plot.margin = margin(10, 120, 10, 10)
    )
  
  ggsave(filename = filename, plot = p, width = 20, height = 20, dpi = 300)
  
  return(p)
}

##### Helper function: optional domain subset + family plot #####

make_domain_family_plot <- function(ps_input, domain_name = NULL, filename) {
  
  if (!inherits(ps_input, "phyloseq")) {
    stop("ps_input must be a phyloseq object.")
  }
  
  if (!"family" %in% rank_names(ps_input)) {
    stop("The taxonomy table does not contain a 'family' rank.")
  }
  
  ps_domain <- ps_input
  
  if (!is.null(domain_name)) {
    if (!"domain" %in% rank_names(ps_input)) {
      stop("The taxonomy table does not contain a 'domain' rank.")
    }
    
    domain_values <- as.character(tax_table(ps_input)[, "domain"])
    
    keep_taxa <- taxa_names(ps_input)[
      !is.na(domain_values) & domain_values == domain_name
    ]
    
    if (length(keep_taxa) == 0) {
      stop(paste0("No taxa found for domain: ", domain_name))
    }
    
    ps_domain <- prune_taxa(keep_taxa, ps_input)
  }
  
  # Replace missing family assignments before agglomerating
  tax_df <- as.data.frame(tax_table(ps_domain), stringsAsFactors = FALSE)
  
  tax_df$family <- as.character(tax_df$family)
  tax_df$family[
    is.na(tax_df$family) |
      tax_df$family == "" |
      tax_df$family == "NA"
  ] <- "Unclassified family"
  
  tax_table(ps_domain) <- tax_table(as.matrix(tax_df))
  
  # Retain taxa with missing ranks as an additional safeguard
  ps_domain_family <- tax_glom(
    ps_domain,
    taxrank = "family",
    NArm = FALSE
  )
  
  if (ntaxa(ps_domain_family) == 0) {
    stop("No taxa remain after family-level aggregation.")
  }
  
  read_counts <- data.frame(
    Sample = sample_names(ps_domain_family),
    total_reads = as.numeric(sample_sums(ps_domain_family))
  ) %>%
    dplyr::mutate(
      read_label = paste0(
        "total reads=",
        scales::comma(total_reads)
      )
    )
  
  ps_domain_family_RA <- microbiome::transform(
    ps_domain_family,
    "compositional"
  )
  
  df_domain <- phyloseq::psmelt(ps_domain_family_RA) %>%
    dplyr::left_join(read_counts, by = "Sample") %>%
    dplyr::mutate(
      family = as.character(family),
      family = dplyr::if_else(
        is.na(family) | family == "",
        "Unclassified family",
        family
      )
    )
  
  make_tax_plot_top20_overall(
    df_input = df_domain,
    taxon_col = "family",
    filename = filename
  )
}

##### Plot 1: Eukaryotes vs Prokaryotes #####

ps.family <- tax_glom(ps, taxrank = "family")

read_counts_all <- data.frame(
  Sample = names(sample_sums(ps.family)),
  total_reads = as.numeric(sample_sums(ps.family))
) %>%
  mutate(read_label = paste0("total reads=", scales::comma(total_reads)))

ps.family_RA <- microbiome::transform(ps.family, "compositional")

df_all <- psmelt(ps.family_RA) %>%
  left_join(read_counts_all, by = "Sample") %>%
  mutate(domain = as.character(domain))

df_domain <- df_all %>%
  mutate(
    domain_group = case_when(
      domain == "d__Eukaryota" ~ "Eukaryotes",
      domain %in% c("d__Bacteria", "d__Archaea") ~ "Prokaryotes",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(Sample, depth, domain_group, total_reads, read_label) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

depth_order_domain <- df_domain %>%
  distinct(depth) %>%
  arrange(desc(depth)) %>%
  pull(depth)

df_domain <- df_domain %>%
  mutate(
    depth_plot = factor(depth, levels = depth_order_domain),
    domain_group = factor(domain_group, levels = c("Eukaryotes", "Prokaryotes", "Other"))
  )

read_labels_domain <- df_domain %>%
  distinct(depth_plot, total_reads, read_label)

domain_cols <- c(
  "Eukaryotes" = "#E41A1C",
  "Prokaryotes" = "#377EB8",
  "Other" = "#D9D9D9"
)

plot1_domain <- ggplot(df_domain, aes(x = depth_plot, y = Abundance, fill = domain_group)) +
  geom_col(width = 0.8) +
  geom_text(
    data = read_labels_domain,
    aes(x = depth_plot, y = 1.03, label = read_label),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold",
    hjust = 0
  ) +
  labs(x = "Depth (cm)", y = "Relative abundance") +
  #scale_x_continuous(breaks=c(0,50,100)) + 
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1.18),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = domain_cols) +
  coord_flip(clip = "off") +
  theme(
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 20, face = "bold"),
    legend.position = "bottom",
    panel.background = element_blank(),
    plot.margin = margin(10, 120, 10, 10)
  )
#

plot1_domain

ggsave(
  filename = file.path(
    output_dir,
    sprintf("plot1_eukaryotes_prokaryotes_k%s_r%s.png", k_value, r_value)
  ),
  plot = plot1_domain,
  width = 20,
  height = 20,
  dpi = 300
)

##### Plot 2: Bacteria + Archaea together #####

ps_prokaryotes <- prune_taxa(
  taxa_names(ps)[tax_table(ps)[, "domain"] %in% c("d__Bacteria", "d__Archaea")],
  ps
)

plot2_prokaryotes <- make_domain_family_plot(
  ps_input = ps_prokaryotes,
  domain_name = NULL,
  filename = file.path(
    output_dir,
    sprintf("plot2_bacteria_archaea_family_top20_overall_k%s_r%s.png", k_value, r_value)
  )
)

plot2_prokaryotes

##### Plot 3: Eukaryotes only #####

plot3_eukaryotes <- make_domain_family_plot(
  ps_input = ps,
  domain_name = "d__Eukaryota",
  filename = filename = file.path(
    output_dir,
    sprintf("plot3_eukaryotes_family_top20_overall_k%s_r%s.png", k_value, r_value)
  )
)

plot3_eukaryotes
