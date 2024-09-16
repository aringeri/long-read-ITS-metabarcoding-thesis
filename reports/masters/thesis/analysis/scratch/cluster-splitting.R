library(phyloseq)
library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(scales)
library(tidytext)

samplesheet <- read.csv('../../../../../experiments/66-fungal-isolate-ONT/Rep1_Run2_ITS_Fungal_database_samplesheet.csv', header = TRUE, row.names = 1)
rownames(samplesheet) <- paste0("barcode", rownames(samplesheet))

nanoclust_otus <-
  read.csv('../../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14/hdbscan_clustering/FULL_ITS/2000/1/otu_table/otu_table.tsv', sep = '\t', row.names = 1) %>%
    otu_table(taxa_are_rows = TRUE) %>%
    phyloseq(sample_data(samplesheet))

vsearch_otus <-
  readRDS('../../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14/phyloseq/FULL_ITS/2000/1/all_samples/all_samples.phyloseq.rds') %>%
    phyloseq(sample_data(samplesheet))
taxa_names(vsearch_otus) <- 1:ntaxa(vsearch_otus)

plot_splitting <- function(phylo, n) {
  phylo %>%
    otu_table() %>%
    data.frame() %>%
    tibble::rownames_to_column('OTU') %>%
    mutate_at(vars(OTU), function(cluster) {
      cluster[cluster == -1] <- "U"
      cluster
    }) %>%
    pivot_longer(cols = sample_names(phylo), names_to='barcode', values_to='count' ) %>%
    filter(count > 0) %>%
    filter(barcode %in% sample_names(phylo)[n]) %>%
    mutate(
      OTU = tidytext::reorder_within(OTU, -count, barcode)
    ) %>%
    ggplot(
      aes(x=OTU, y=count)#, colour=OTU)
    ) +
    geom_col() +
    facet_wrap(~barcode,
               labeller = as_labeller(\(x)  paste0(
                 samplesheet[x, 'Sample'], ' (', gsub('barcode','', x), ')'
               )),
               scales='free_x') +
    # coord_flip() +
    tidytext::scale_x_reordered() +
    labs(x="Cluster ID", y="Abundance")
}

first_30_nano <- plot_splitting(nanoclust_otus, 1:30)
second_30_nano <- plot_splitting(nanoclust_otus, 31:nsamples(nanoclust_otus))
ggsave('../images/06-cluster-splitting-nanoclust-1-30.png', first_30_nano)
ggsave('../images/06-cluster-splitting-nanoclust-31-60.png', second_30_nano)


filter_taxa_by_thresh <- function(phylo, thresh) {
  t <- sum(sample_sums(phylo))*thresh
  filter_taxa(phylo, \(x) sum(x) > t, prune=TRUE)
}
plot_splitting(vsearch_otus, c(1, 3:7))
plot_splitting(filter_taxa_by_thresh(vsearch_otus, 0.001), c(1, 3:7))
