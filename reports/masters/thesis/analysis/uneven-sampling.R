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

low_sampled_barcodes <- c("barcode77" = 0.05, "barcode35" = 0.1, "barcode39" = 0.2, "barcode71" = 0.4, "barcode57" = 0.5)
min_cluster_sizes <- c(0, 0.001, 0.0025, 0.005, 0.0065)

combined_data <- tibble()
for (min_cluster_size in min_cluster_sizes) {
  sample_depth <- config$sample_depth
  phylo <- load_nanoclust_phyloseq_3(
    samplesheet,
    experiment = '../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-uneven-reps2-09-16', reads_per_sample = sample_depth,
    repetition = config$repetition, min_cluster_size=min_cluster_size)

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
                            as.character(low_sampled_barcodes[barcode]*sample_depth),
                            as.character(sample_depth))
    ) %>%
    ungroup()
  data$min_cluster_size <- min_cluster_size

  combined_data <- rbind(combined_data, data)
}

library_size <- 50*config$sample_depth + sum(config$sample_depth*low_sampled_barcodes)

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
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  # tidytext::scale_x_reordered() +
  # scale_fill_discrete(name="sample depth") +
  scale_fill_manual(name="sample depth", values = c("#E87D72", "#B39F33", "#53B64C", "#55BCC2", "#6E9CF8", "grey69")) +
  scale_colour_discrete(guide="none")

# phylo <- load_nanoclust_phyloseq_3(
#   samplesheet,
#   experiment = '../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-uneven-reps2-09-16', reads_per_sample = config$sample_depth,
#   repetition = config$repetition, min_cluster_size=0.0025)
#
#
# data <- otu_table(phylo) %>%
#   data.frame() %>%
#   tibble::rownames_to_column(var = 'OTU') %>%
#   pivot_longer(barcode25:barcode89, names_to = "barcode", values_to = 'count') %>%
#   filter(count > 0) %>%
#   filter(barcode %in% names(low_sampled_barcodes)) %>%
#   group_by(barcode) %>%
#   mutate(
#     clusterSize = sum(count),
#     numOTUs = length(OTU),
#     sample_count = ifelse(barcode %in% names(low_sampled_barcodes),
#                           as.character(low_sampled_barcodes[barcode]*2000),
#                           "2000")
#   ) %>%
#   ungroup()
# data %>%
#   ggplot(aes(
#     x=reorder(barcode,clusterSize),
#     y=count,
#     # fill=barcode %in% names(low_sampled_barcodes)
#     fill=reorder(sample_count, as.numeric(sample_count))
#   )) +
#   geom_bar(stat='identity') +
#   # geom_label(aes(label=numOTUs), data=(data %>% group_by(barcode))) +
#   # geom_vline(xintercept = 0.001 * library_size) +
#   geom_hline(yintercept = min_cluster_sizes * library_size, linetype='dashed') +
#   # scale_y_continuous(breaks = min_cluster_sizes*library_size)
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   scale_fill_discrete(name="sample depth")
# # scale_x_log10() + scale_y_log10()
