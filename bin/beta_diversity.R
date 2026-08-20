#!/usr/bin/env Rscript

library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)

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

otu_df <- read.delim(
  otu_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

otu_mat <- as.matrix(otu_df)

otu_ps <- otu_table(otu_mat, taxa_are_rows = TRUE)

tax_table_df <- read.delim(
  tax_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  quote = ""
)

metadata_df <- read.delim(
    metadata_file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
)

metadata <- sample_data(metadata_df)

ps <- phyloseq(otu_ps2, tax_ps, metadata)

##### Optional: agglomerate to genus level #####
# Only do this if your taxonomy table has a column called "genus"
ps_genus <- tax_glom(ps, taxrank = "genus", NArm = FALSE)

##### Remove taxa with zero reads #####

ps_genus <- prune_taxa(taxa_sums(ps_genus) >= 1, ps_genus)

##### Convert to relative abundance per sample #####
ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

##### Extract community matrix #####
otu_rel <- as.data.frame(otu_table(ps_rel))

# vegan needs samples as rows and taxa as columns
if (taxa_are_rows(ps_rel)) {
  otu_rel <- t(otu_rel)
}

##### Remove samples with zero reads, just in case #####
otu_rel <- otu_rel[rowSums(otu_rel) > 0, ]

##### Bray-Curtis distance #####
bray_dist <- vegdist(otu_rel, method = "bray")

##### NMDS #####
set.seed(123)
nmds <- metaMDS(bray_dist, k = 3, trymax = 200)

##### Prepare plot data #####
nmds_df <- as.data.frame(nmds$points)
nmds_df$SampleID <- rownames(nmds_df)

##### Add metadata #####
metadata_df <- as.data.frame(sample_data(ps_rel))
metadata_df$SampleID <- rownames(metadata_df)

nmds_df <- left_join(nmds_df, metadata_df, by = "SampleID")

##### Add total sequences per sample #####
seq_counts <- sample_sums(ps_genus)

nmds_df$Sequences <- seq_counts[match(nmds_df$SampleID, names(seq_counts))]

##### Add observed OTUs / genera per sample #####
observed_otus <- estimate_richness(ps_genus, measures = "Observed")

observed_otus$SampleID <- rownames(observed_otus)
observed_otus$SampleID <- gsub("X", "", observed_otus$SampleID)

nmds_df <- left_join(nmds_df, observed_otus, by = "SampleID")

##### Label: depth + sequences + observed OTUs #####
nmds_df$Label <- paste0(
  nmds_df$depth,
  " cm (",
  comma(nmds_df$Sequences),
  " seqs; ",
  nmds_df$Observed,
  " genera)"
)


##### Make plot #####
nmds_plot <- ggplot(nmds_df, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(colour = depth), size = 4) +
  geom_text_repel(
    aes(label = Label),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    segment.color = "grey60",
    force = 2
  ) +
  theme_classic() +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    colour = "Depth (cm)",
    title = "NMDS of metagenomic taxonomic profiles at genus level",
    subtitle = paste0("Stress = ", round(nmds$stress, 3))
  )

nmds_plot

ggsave(
  filename = file.path(
    output_dir,
    sprintf("NMDS_metagenomic_genus_depth_k%s_r%s.png", k_value, r_value)
  ),
  plot = nmds_plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
