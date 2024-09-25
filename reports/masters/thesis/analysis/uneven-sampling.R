library(phyloseq)
library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(scales)

source('./helpers/config.R')
source('./helpers/vsearch.R')
source('./helpers/dnabarcoder.R')

samplesheet <- read_samplesheet(config)

low_sampled_barcodes <- c("barcode77" = 5, "barcode35" = 10, "barcode39" = 20, "barcode71" = 50, "barcode57" = 100)
min_cluster_sizes <- c(0, 0.005)
sample_depth <- 200
repetition <- 1

combined_data <- tibble()
for (min_cluster_size in min_cluster_sizes) {
  phylo <- load_nanoclust_phyloseq_3(
    samplesheet,
    experiment = '../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-uneven-reps-lo-09-24', reads_per_sample = sample_depth,
    repetition = repetition, min_cluster_size=min_cluster_size)

  data <- otu_table(phylo) %>%
    data.frame() %>%
    tibble::rownames_to_column(var = 'OTU') %>%
    pivot_longer(barcode25:barcode89, names_to = "barcode", values_to = 'count') %>%
    filter(count > 0) %>%
    # filter(barcode %in% names(low_sampled_barcodes)) %>%
    group_by(OTU) %>%
    mutate(
      clusterSize = sum(count),
      # numOTUs = length(OTU),
      sample_count = ifelse(barcode %in% names(low_sampled_barcodes),
                            as.character(low_sampled_barcodes[barcode]),
                            as.character(sample_depth))
    ) %>%
    ungroup()
  data$min_cluster_size <- min_cluster_size

  combined_data <- rbind(combined_data, data)
}

library_size <- 50*sample_depth + sum(low_sampled_barcodes)

combined_data_2 <- combined_data %>%
  mutate(barcode = sub('barcode', '', barcode)) %>%
  mutate( OTU = tidytext::reorder_within(OTU, clusterSize, min_cluster_size) )
  # mutate( OTU = tidytext::reorder_within(OTU, clusterSize, min_cluster_size) )

uneven_nanoclust <- combined_data_2 %>%
  # filter(min_cluster_size %in% c(0, 0.001)) %>%
  ggplot(aes(
    x=OTU,#reorder(OTU, clusterSize),#reorder(barcode,clusterSize),
    y=count,
    # fill=barcode %in% names(low_sampled_barcodes)
    fill=reorder(sample_count, as.numeric(sample_count))
  )) +
  geom_bar(stat='identity', color='grey42') +
  geom_hline(aes(yintercept = min_cluster_size * library_size), linetype='dashed') +
  # geom_text(aes(label=round(min_cluster_size * library_size), x=52, y=round(min_cluster_size * library_size + 100)), data = combined_data %>% group_by(min_cluster_size)) +
  facet_wrap(ncol=2, ~min_cluster_size,
             labeller = as_labeller(\(x) paste0('minimum ', ifelse(x == 0, "2", as.numeric(x)*library_size), " reads (", x, ")")),
             scales='free_x'
  ) +
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1),
    # axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  tidytext::scale_x_reordered() +
  # scale_fill_discrete(name="sample depth") +
  scale_fill_manual(name="sample depth", values = c("#E87D72", "#B39F33", "#53B64C", "#55BCC2", "#6E9CF8", "grey69")) +
  scale_colour_discrete(guide="none")
uneven_nanoclust


vsearch <- load_vsearch_phyloseq(
  samplesheet = samplesheet,
  experiment = '../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-uneven-reps-lo-09-24',
  reads_per_sample = sample_depth,
  repetition = repetition)


data <- vsearch %>%
  # filter_taxa_by_thresh(0.0001) %>%
  otu_table() %>%
  data.frame() %>%
  tibble::rownames_to_column(var = 'OTU') %>%
  pivot_longer(barcode25:barcode89, names_to = "barcode", values_to = 'count') %>%
  filter(count > 0) %>%
  # filter(barcode %in% names(low_sampled_barcodes)) %>%
  group_by(OTU) %>%
  mutate(
    clusterSize = sum(count),
    # numOTUs = length(OTU),
    sample_count = ifelse(barcode %in% names(low_sampled_barcodes),
                          as.character(low_sampled_barcodes[barcode]),
                          as.character(sample_depth))
  ) %>%
  ungroup()

data %>%
  mutate( OTU = tidytext::reorder_within(OTU, clusterSize, min_cluster_size) ) %>%
  # filter(min_cluster_size %in% c(0, 0.001)) %>%
  ggplot(aes(
    x=OTU,#reorder(OTU, clusterSize),#reorder(barcode,clusterSize),
    y=count,
    # fill=barcode %in% names(low_sampled_barcodes)
    fill=reorder(sample_count, as.numeric(sample_count))
  )) +
  geom_bar(stat='identity', color='grey42') +
  # geom_hline(aes(yintercept = min_cluster_size * library_size), linetype='dashed') +
  # geom_text(aes(label=round(min_cluster_size * library_size), x=52, y=round(min_cluster_size * library_size + 100)), data = combined_data %>% group_by(min_cluster_size)) +
  # facet_wrap(ncol=2, ~min_cluster_size,
  #            labeller = as_labeller(\(x) paste0('minimum ', ifelse(x == 0, "2", as.numeric(x)*library_size), " reads (", x, ")")),
  #            scales='free_x'
  # ) +
  theme(
    # text = element_text(size = 10),
    # axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  # tidytext::scale_x_reordered() +
  # scale_fill_discrete(name="sample depth") +
  scale_fill_manual(name="sample depth", values = c("#E87D72", "#B39F33", "#53B64C", "#55BCC2", "#6E9CF8", "grey69")) +
  scale_colour_discrete(guide="none")

