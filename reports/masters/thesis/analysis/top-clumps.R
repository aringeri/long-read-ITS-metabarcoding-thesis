library(phyloseq)
library(tidyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(glue)
library(scales)
library(tidytext)
library(microViz)
library(stringr)

source('helpers/dnabarcoder.R')
source('helpers/vsearch.R')
source('helpers/config.R')
source('helpers/manual-taxonomy.R')

samplesheet <- read_samplesheet(config)

experiment_path <- config$experiment_path
sample_depth <- config$sample_depth
repetition <- config$repetition
min_cluster_size <- 0

load_dnabarcoder_classification <- function(
  experiment="../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-09-12",
  sequence_type="nanoclust_abundant",
  reads_per_sample=2000,
  repetition=2,
  min_cluster_size=0
) {

  path <- if (!is.null(min_cluster_size)) {
    glue('{experiment}/dnabarcoder/{sequence_type}/FULL_ITS/{reads_per_sample}/{repetition}/{min_cluster_size}/classify/')
  } else {
    glue('{experiment}/dnabarcoder/{sequence_type}/FULL_ITS/{reads_per_sample}/{repetition}/classify/')
  }

  classFile <- list.files(
    path,
    pattern='*.unite2024ITS_BLAST.classification',
    full.names = T
  )
  if (!is.null(min_cluster_size)) {
    read.csv(classFile[1], sep='\t', row.names = 1) %>%
      tibble::rownames_to_column('OTU') %>%
      separate(., OTU, into=c('read', 'barcode', 'cluster', NA, 'size'), sep=';') %>%
      mutate_at(vars(read, barcode, cluster, size), \(str) gsub('.*=', '', str)) %>%
      filter(cluster != -1) %>%
      arrange(as.numeric(cluster)) %>%
      tibble::column_to_rownames("cluster")
  } else {
    read.csv(classFile[1], sep='\t', row.names = 1) %>%
      tibble::rownames_to_column('OTU') %>%
      separate(., OTU, into=c('read', 'barcode', NA, 'size'), sep=';') %>%
      mutate_at(vars(read, barcode, size), \(str) gsub('.*=', '', str)) %>%
      tibble::column_to_rownames("read")
  }
}

load_otus <- function (
  experiment="../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14",
  reads_per_sample=2000,
  repetition=2,
  min_cluster_size=0
) {
  # print(experiment)
  otu <- read.csv(glue('{experiment}/hdbscan_clustering/FULL_ITS/{reads_per_sample}/{repetition}/{min_cluster_size}/otu_table/otu_table.tsv'), sep = '\t', row.names = 1)
  otu <- otu[order(as.numeric(rownames(otu))), ]
  otu <- otu[rownames(otu) != -1, ]
  otu
}

classifications <- load_dnabarcoder_classification(
  experiment = experiment_path,
  sequence_type = 'nanoclust_abundant',
  reads_per_sample = sample_depth,
  repetition = repetition,
  min_cluster_size = min_cluster_size) %>%
  rename(read_barcode = 'barcode')

otus <- load_otus(
  experiment = experiment_path,
  reads_per_sample = sample_depth,
  repetition = repetition,
  min_cluster_size = min_cluster_size)

otu_tax_and_actual <- merge(classifications, otus, by=0) %>%
  pivot_longer(barcode25:barcode89, names_to = 'barcode', values_to = 'count') %>%
  filter(count > 0) %>%
  rename(otu = 'Row.names') %>%
  left_join(
    (samplesheet %>%
      tibble::rownames_to_column("barcode") %>%
      select(barcode, order_unite, family_unite, genus_unite, species) %>%
      rename(
        actual_order = order_unite,
        actual_family = family_unite,
        actual_genus = genus_unite,
        actual_species = species
      )
    ),
    join_by(barcode)
  )

major_clumps_nc <- otu_tax_and_actual %>%
  mutate(size = as.numeric(size)) %>%
  group_by(otu) %>%
  filter(all(sum(count > 300) > 1)) %>%
  filter(count > 11) %>%
  select(otu, genus, species, count, barcode, actual_genus, actual_species) %>%
  arrange(otu, actual_genus, actual_species)

View(major_clumps_nc)

vsearch_classifications <- load_dnabarcoder_classification(
  experiment = experiment_path,
  sequence_type = 'vsearch',
  reads_per_sample = sample_depth,
  repetition = repetition,
  min_cluster_size = NULL) %>%
  rename(read_barcode = 'barcode')

vsearch_otus <- readRDS(glue('{experiment_path}/phyloseq/FULL_ITS/{sample_depth}/{repetition}/all_samples/all_samples.phyloseq.rds'))

vsearch_merged <- merge(vsearch_classifications, vsearch_otus, by=0) %>%
  pivot_longer(barcode25:barcode89, names_to = 'barcode', values_to = 'count') %>%
  filter(count > 0) %>%
  rename(otu = 'Row.names') %>%
  left_join(
    (samplesheet %>%
      tibble::rownames_to_column("barcode") %>%
      select(barcode, order_unite, family_unite, genus_unite, species) %>%
      rename(
        actual_order = order_unite,
        actual_family = family_unite,
        actual_genus = genus_unite,
        actual_species = species
      )
    ),
    join_by(barcode)
  )

vsearch_merged %>%
  mutate(size = as.numeric(size)) %>%
  group_by(otu) %>%
  filter(all(sum(count > 300) > 1)) %>%
  filter(count > 15) %>%
  select(otu, genus, species, count, barcode, actual_genus, actual_species) %>%
  arrange(otu, actual_genus, actual_species) %>%
  View()

