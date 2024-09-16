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

source('../helpers/dnabarcoder.R')

samplesheet <- read.csv('../../../../../experiments/66-fungal-isolate-ONT/Rep1_Run2_ITS_Fungal_database_samplesheet.csv', header = TRUE, row.names = 1)
rownames(samplesheet) <- paste0("barcode", rownames(samplesheet))

experiment <- "../../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14"

# load_nanoclust_phyloseq <- function(samplesheet, experiment) {
#   otu <- read.csv(glue('{experiment}/hdbscan_clustering/FULL_ITS/2000/2/otu_table/otu_table.tsv'), sep = '\t', row.names = 1)
#   otu <- otu[order(as.numeric(rownames(otu))), ]
#   otu <- otu[rownames(otu) != -1, ]
#
#   tax <- read.csv(glue('{experiment}/dnabarcoder/nanoclust/FULL_ITS/2000/2/classify/all_samples.unite2024ITS_BLAST.classification'), sep='\t', row.names = 1) %>%
#     tibble::rownames_to_column('OTU') %>%
#     separate(OTU, into=c('read', 'barcode', 'cluster', NA, 'size'), sep=';') %>%
#     mutate_at(vars(read, barcode, cluster, size), \(str) gsub('.*=', '', str)) %>%
#     filter(cluster != -1) %>%
#     arrange(as.numeric(cluster)) %>%
#     tibble::column_to_rownames("cluster")
#
#   phyloseq(
#     otu_table(otu, taxa_are_rows = TRUE),
#     tax_table(as.matrix(tax[, !(names(tax) %in% c('read', 'barcode', 'size', 'ReferenceID', 'rank', 'score', 'cutoff', 'confidence'))])),
#     sample_data(samplesheet)
#   )
# }
nanoclust <- load_nanoclust_phyloseq(samplesheet, experiment, 2000, 2) %>%
  tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)

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
      labels = \(x) sample_data(phylo)[reorder_func(x), ][['Sample']]
      # labels = \(x) sample_data(phylo)[reorder_func(x), 'Species']
    ) +
    scale_fill_discrete(guide='none') +
    labs(x="Sample", y="Abundance") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

plot_clumping(nanoclust, 1:6)
plot_clumping(nanoclust, c(43, 54, 55) + 1)

ggsave('../images/06-cluster-clumping-nanoclust-with-tax-ex.png', plot_clumping(nanoclust, c(57, 56, 34, 22) + 1))
