library(phyloseq)
library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(scales)
library(tidytext)
library(microViz)
library(ggpattern)

source('../helpers/dnabarcoder.R')
source('../helpers/config.R')
source('../helpers/vsearch.R')

samplesheet <- read_samplesheet(config)

nanoclust <- load_nanoclust_phyloseq_2(samplesheet=samplesheet, experiment=config$experiment_path, reads_per_sample=config$sample_depth, repetition=config$repetition) %>%
  tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)

umap <- read.csv(
  glue('{config$experiment_path}/umap_transform/FULL_ITS/{config$sample_depth}/{config$repetition}/all_samples/umap.output.tsv'),
  sep = '\t')

hdbscan <- read.csv(
  glue('{config$experiment_path}/hdbscan_clustering/FULL_ITS/{config$sample_depth}/{config$repetition}/clusters/hdbscan.output.tsv'),
  sep = '\t')

together <- left_join(umap, hdbscan, join_by(read), relationship = 'many-to-one')

# together[together[,"cluster_id"] != -1, ] %>%
together %>%
  tibble() %>%
  tidyr::separate(read, sep=';barcodelabel=', into=c("read id", "barcode")) %>%
  mutate(barcode = sub(';', '', barcode)) %>%
  filter(cluster_id != -1) %>%
  filter(barcode %in% paste0('barcode', c(25, 27, 28, 36))) %>%
  filter(cluster_id != -1) %>%
  ggplot(aes(x=D1, y=D2, colour = factor(paste0(cluster_id, " - ", tax_table(nanoclust)[as.character(cluster_id), 'species'])))) +
  # facet_wrap(~samplesheet[barcode, 'UpdatedName']) +
  facet_wrap(~samplesheet[barcode, 'UpdatedName']) +
  # ggplot(aes(x=D1, y=D2, colour = factor(cluster_id))) +
  geom_point(size=0.2, alpha=0.5) +
  geom_text(
    aes(label=cluster_id),
    data=together %>%
      tibble() %>%
      tidyr::separate(read, sep=';barcodelabel=', into=c("read id", "barcode")) %>%
      mutate(barcode = sub(';', '', barcode)) %>%
      filter(cluster_id != -1) %>%
      filter(barcode %in% paste0('barcode', c(25, 27, 28, 36))) %>%
      group_by(barcode, cluster_id) %>%
      summarise_at(vars(D1, D2), mean),
    size=3, colour='black')
  # geom_label(aes(label=cluster_id))


together %>%
  tibble() %>%
  tidyr::separate(read, sep=';barcodelabel=', into=c("read id", "barcode")) %>%
  mutate(barcode = sub(';', '', barcode)) %>%
  # filter(cluster_id != -1) %>%
  filter(barcode %in% paste0('barcode', 59:62)) %>%
  filter(cluster_id != -1) %>%
  ggplot(aes(x=D1, y=D2, colour = factor(tax_table(nanoclust)[as.character(cluster_id), 'species']))) +
  facet_wrap(~samplesheet[barcode, 'UpdatedName']) +
  # ggplot(aes(x=D1, y=D2, colour = factor(cluster_id))) +
  geom_point(size=0.2, alpha=0.5)