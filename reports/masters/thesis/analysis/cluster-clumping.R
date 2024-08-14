library(phyloseq)
library(tidyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(glue)
library(scales)
library(tidytext)

samplesheet <- read.csv('../../../../experiments/66-fungal-isolate-ONT/Rep1_Run2_ITS_Fungal_database_samplesheet.csv', header = TRUE, row.names = 1)
rownames(samplesheet) <- paste0("barcode", rownames(samplesheet))

nanoclust_otus <-
  read.csv('../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-08/hdbscan_clustering/FULL_ITS/2000/1/otu_table/otu_table.tsv', sep = '\t', row.names = 1) %>%
    otu_table(taxa_are_rows = TRUE) %>%
    phyloseq(sample_data(samplesheet))

vsearch_otus <-
  readRDS('../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-08/phyloseq/FULL_ITS/2000/1/all_samples/all_samples.phyloseq.rds') %>%
    phyloseq(sample_data(samplesheet))
taxa_names(vsearch_otus) <- 1:ntaxa(vsearch_otus)

plot_clumping <- function(phylo, n) {
  phylo %>%
    otu_table() %>%
    data.frame() %>%
    tibble::rownames_to_column('OTU') %>%
    filter(OTU %in% rownames(otu_table(phylo))[n]) %>%
    filter(OTU != -1) %>%
    pivot_longer(cols = sample_names(phylo), names_to='barcode', values_to='count' ) %>%
    filter(count > 0) %>%
    mutate(barcode = gsub('barcode', '', barcode)) %>%
    mutate(
      barcode = tidytext::reorder_within(barcode, -count, OTU)
    ) %>%
    # View() %>%
    ggplot(
      aes(x=barcode, y=count)#, colour=OTU)
    ) +
    geom_col() +
    facet_wrap(~factor(OTU, levels = sort(as.integer(rownames(otu_table(phylo))))),
               labeller = as_labeller(\(x)  paste0("Cluster ", x)),
               scales='free_x') +
    # coord_flip() +
    tidytext::scale_x_reordered() +
    labs(x="Barcode ID", y="Abundance") #+
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

first_30_clump_nano <- plot_clumping(nanoclust_otus, 1:30)
second_30_clump_nano <- plot_clumping(nanoclust_otus, 31:ntaxa(nanoclust_otus))
ggsave('images/06-cluster-clumping-nanoclust-1-30.png', first_30_clump_nano)
ggsave('images/06-cluster-clumping-nanoclust-31-60.png', second_30_clump_nano)

plot_clumping(vsearch_otus, c(750, 128, 465, 193, 706, 1536 ))