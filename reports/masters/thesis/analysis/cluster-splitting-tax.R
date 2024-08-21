library(phyloseq)
library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(scales)
library(tidytext)
library(microViz)

samplesheet <- read.csv('../../../../experiments/66-fungal-isolate-ONT/Rep1_Run2_ITS_Fungal_database_samplesheet.csv', header = TRUE, row.names = 1)
rownames(samplesheet) <- paste0("barcode", rownames(samplesheet))

experiment <- "../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14"

load_nanoclust_phyloseq <- function(samplesheet, experiment) {
  otu <- read.csv(glue('{experiment}/hdbscan_clustering/FULL_ITS/2000/2/otu_table/otu_table.tsv'), sep = '\t', row.names = 1)
  otu <- otu[order(as.numeric(rownames(otu))), ]
  otu <- otu[rownames(otu) != -1, ]

  tax <- read.csv(glue('{experiment}/dnabarcoder/nanoclust/FULL_ITS/2000/2/classify/all_samples.unite2024ITS_BLAST.classification'), sep='\t', row.names = 1) %>%
    tibble::rownames_to_column('OTU') %>%
    separate(OTU, into=c('read', 'barcode', 'cluster', NA, 'size'), sep=';') %>%
    mutate_at(vars(read, barcode, cluster, size), \(str) gsub('.*=', '', str)) %>%
    filter(cluster != -1) %>%
    arrange(as.numeric(cluster)) %>%
    tibble::column_to_rownames("cluster")

  phyloseq(
    otu_table(otu, taxa_are_rows = TRUE),
    tax_table(as.matrix(tax[, !(names(tax) %in% c('read', 'barcode', 'size', 'ReferenceID', 'rank', 'score', 'cutoff', 'confidence'))])),
    sample_data(samplesheet)
  )
}
nanoclust <- load_nanoclust_phyloseq(samplesheet, experiment) %>%
  tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)


plot_splitting <- function(phylo, n) {
  cbind(
    otu_table(phylo),
    tax_table(phylo)[,'species']
  ) %>% data.frame() %>%
    tibble::rownames_to_column('OTU') %>%
    pivot_longer(cols = sample_names(phylo), names_to='barcode', values_to='count' ) %>%
    mutate(count = as.integer(count)) %>%
    filter(count > 0) %>%
    filter(barcode %in% sample_names(phylo)[n]) %>%
    mutate(
      OTU = tidytext::reorder_within(OTU, -count, barcode)
    ) %>%
    ggplot(
      aes(x=OTU, y=count, fill=species)
    ) +
    geom_col() +
    facet_wrap(~barcode,
               labeller = as_labeller(\(x)  paste0(
                 samplesheet[x, 'Sample'], ' (', gsub('barcode','', x), ')'
               )),
               scales='free_x') +
    # coord_flip() +
    tidytext::scale_x_reordered(
      labels = \(x) paste0(tax_table(phylo)[reorder_func(x), 'species'], " (OTU ", reorder_func(x), ")")
    ) +
    scale_fill_discrete(guide='none') +
    labs(x="Cluster", y="Abundance") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

nanoclust_splitting <- plot_splitting(nanoclust, 13:24)
plot_splitting(nanoclust, 32:45)
ggsave('images/06-cluster-splitting-nanoclust-with-tax-6.png', plot_splitting(nanoclust, 1:6))
ggsave('images/06-cluster-splitting-nanoclust-with-tax-ex.png', plot_splitting(nanoclust, 2))

