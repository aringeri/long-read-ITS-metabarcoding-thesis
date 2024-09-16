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

source('../helpers/config.R')
source('../helpers/dnabarcoder.R')
source('../helpers/vsearch.R')

samplesheet <- read_samplesheet(config)

nanoclust <- load_nanoclust_phyloseq_3(
  samplesheet=samplesheet, experiment=config$experiment_path, reads_per_sample=config$sample_depth, repetition=config$repetition,
  sequence_type='nanoclust_abundant', min_cluster_size = '0.005') %>%
  tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)

samplesheet <- read.csv('../../../../../experiments/66-fungal-isolate-ONT/Rep1_Run2_ITS_Fungal_database_samplesheet.csv', header = TRUE, row.names = 1)
rownames(samplesheet) <- paste0("barcode", rownames(samplesheet))

experiment <- "../../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14"


vsearch_otus <- readRDS(glue('{experiment}/phyloseq/FULL_ITS/2000/1/all_samples/all_samples.phyloseq.rds'))
tax_table <- read_dna_barcoder_classification_vsearch(glue('{experiment}/dnabarcoder/vsearch/FULL_ITS/2000/1/classify/all_samples.unite2024ITS_BLAST.classification'))

vsearch_phylo <- phyloseq(vsearch_otus, tax_table, sample_data(samplesheet))


filter_taxa_by_thresh(vsearch_phylo, 0.001)

plot_clumping <- function(phylo, n) {
  tax <- tax_table(phylo)
  phylo %>%
    otu_table() %>%
    data.frame() %>%
    tibble::rownames_to_column('OTU') %>%
    filter(OTU %in% rownames(otu_table(phylo))[n]) %>%
    filter(OTU != -1) %>%
    pivot_longer(cols = sample_names(phylo), names_to='barcode', values_to='count' ) %>%
    filter(count > 0) %>%
    mutate(barcode = barcode) %>%
    mutate(
      barcode = tidytext::reorder_within(barcode, -count, OTU)
    ) %>%
    # View() %>%
    ggplot(
      aes(x=barcode, y=count, fill=barcode)
    ) +
    geom_col() +
    facet_wrap(~factor(OTU,
                       # levels = sort(as.integer(rownames(otu_table(phylo))))
               ),
               labeller = as_labeller(\(x) paste0(tax[as.character(x), 'species'], " (OTU ", x, ")")),
               scales='free_x') +
    # coord_flip() +
    tidytext::scale_x_reordered(
      # labels = \(x) sample_data(phylo)[reorder_func(x), ][['Sample']]
      # labels = \(x) sample_data(phylo)[reorder_func(x), 'Species']
    ) +
    scale_fill_discrete(guide='none') +
    labs(x="Sample", y="Abundance") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

filter_taxa_by_thresh(vsearch_phylo, 0.001)
plot_clumping(filter_taxa_by_thresh(vsearch_phylo, 0.001), 1:30)
plot_clumping(nanoclust, c(57, 49) + 1)

ggsave('../images/06-cluster-clumping-nanoclust-with-tax-ex.png', plot_clumping(nanoclust, c(57, 56, 34, 22) + 1))
