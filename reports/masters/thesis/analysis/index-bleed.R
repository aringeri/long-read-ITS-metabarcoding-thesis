library(phyloseq)
library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(scales)

samplesheet <- read.csv('../../../../experiments/66-fungal-isolate-ONT/Rep1_Run2_ITS_Fungal_database_samplesheet.csv', header = TRUE, row.names = 1)
rownames(samplesheet) <- paste0("barcode", rownames(samplesheet))

otu_table <- readRDS('../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-02/phyloseq/FULL_ITS/1000/1/all_samples/all_samples.phyloseq.rds') %>%
  phyloseq(sample_data(samplesheet), taxa_are_rows=TRUE)

barcode43 <- otu_table %>%
  subset_samples(Sample == 'Aspergillus_fumigatus') %>%
  filter_taxa(\(x) sum(x) > 0, prune=TRUE)

filter_taxa_by_thresh <- function(phylo, thresh) {
  t <- sum(sample_sums(phylo))*thresh
  filter_taxa(phylo, \(x) sum(x) > t, prune=TRUE)
}

cowplot::plot_grid(
  otu_table(barcode43) %>% sort(decreasing = TRUE) %>%
    data.frame() %>%
    tibble::rowid_to_column('OTU') %>%
    ggplot(
      aes(x=OTU, y=barcode43)
    ) +
      geom_col(),
      # scale_y_continuous(transform = 'log10')

  filter_taxa_by_thresh(otu_table, 0.001) %>%
    subset_samples(Sample == 'Aspergillus_fumigatus') %>%
    filter_taxa(\(x) sum(x) > 0, prune=TRUE) %>%
    otu_table %>%
    sort(decreasing = TRUE) %>%
    data.frame() %>%
    tibble::rowid_to_column('OTU') %>%
    ggplot(
      aes(x=OTU, y=barcode43)
    ) +
    geom_col()
  ,
  ncol = 1
)

# otu_table %>%
#   filter_taxa_by_thresh( 0.001) %>%
#   otu_table() %>%
#   data.frame() %>%
#   tibble::rowid_to_column('OTU') %>%
#   pivot_longer(cols = barcode25:barcode47, names_to='barcode', values_to='count' ) %>%
#   ggplot(
#     aes(x=OTU, y=count, colour=factor(OTU))
#   ) +
#   geom_col() +
#   facet_grid(~barcode)

# nanoplot
nanoclust_otus <-
  read.csv('../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-08/hdbscan_clustering/FULL_ITS/2000/1/otu_table/otu_table.tsv', sep = '\t', row.names = 1) %>%
  # read.csv('/Users/alex/repos/long-read-ITS-metabarcoding/output/isolate-sub5000/hdbscan_clustering/FULL_ITS/50/1/otu_table/otu_table.tsv', sep = '\t', row.names = 1) %>%
  otu_table(taxa_are_rows = TRUE) %>%
  phyloseq(sample_data(samplesheet))

nanoclust_otus %>%
  otu_table() %>%
  data.frame() %>%
  tibble::rownames_to_column('OTU') %>%
  # tibble::rownames_to_column('cluster_id') %>%
  mutate_at(vars(OTU), function(cluster) {
    cluster[cluster == -1] <- "U"
    cluster
  }) %>%
  pivot_longer(cols = sample_names(nanoclust_otus)[1:30], names_to='barcode', values_to='count' ) %>%
  filter(count > 0) %>%
  # filter(barcode %in% c('barcode41', 'barcode42', 'barcode43', 'barcode28', 'barcode36', 'barcode47')) %>%
  mutate(
    OTU = tidytext::reorder_within(OTU, -count, barcode)
  ) %>%
  ggplot(
      aes(x=OTU, y=count)#, colour=OTU)
    ) +
    geom_col() +
    facet_wrap(~barcode,
               labeller = as_labeller(\(x)  paste0(
                 samplesheet[x, 'Sample']#, ' (', x, ')'
               )),
               ncol = 6,
               scales='free_x') +
    # coord_flip() +
    tidytext::scale_x_reordered() +
    labs(x="Cluster ID", y="Abundance")