library(phyloseq)
library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(magrittr)


samplesheet <- read.csv('../../../../experiments/66-fungal-isolate-ONT/Rep1_Run2_ITS_Fungal_database_samplesheet.csv', header = TRUE, row.names = 1)
rownames(samplesheet) <- paste0("barcode", rownames(samplesheet))

otu <- read.csv('../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14/hdbscan_clustering/FULL_ITS/2000/2/otu_table/otu_table.tsv', sep='\t', row.names = 1)
otu <- otu[order(as.numeric(rownames(otu))), ]
otu <- otu[rownames(otu) != -1, ]

tax <- read.csv('../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14/dnabarcoder/nanoclust/FULL_ITS/2000/2/classify/all_samples.unite2024ITS_BLAST.classification', sep='\t', row.names = 1) %>%
  tibble::rownames_to_column('OTU') %>%
  separate(OTU, into=c('read', 'barcode', 'cluster', NA, 'size'), sep=';') %>%
  mutate_at(vars(read, barcode, cluster, size), \(str) gsub('.*=', '', str)) %>%
  filter(cluster != -1) %>%
  tibble::column_to_rownames("cluster")

# tax[, !(names(tax) %in% c('read', 'barcode', 'size', 'ReferenceID', 'rank', 'score', 'cutoff', 'confidence'))] %>% rownames()
tax <- tax[order(as.numeric(rownames(tax))), ]

phylo <- phyloseq(
  otu_table(otu, taxa_are_rows = TRUE),
  tax_table(as.matrix(tax[, !(names(tax) %in% c('read', 'barcode', 'size', 'ReferenceID', 'rank', 'score', 'cutoff', 'confidence'))])),
  sample_data(samplesheet)
)

plot_bar(
  subset_samples(phylo,  grepl('Aspergillus', Sample)),
  fill='species',
  x = 'sample_Sample'
)

plot_bar(
  subset_samples(phylo, grepl("Asteroma", Sample)),
  fill='species',
  x = 'sample_Sample'
)
plot_bar(
  subset_samples(phylo, grepl("Austropucc", Sample)),
  fill='species',
  x = 'sample_Sample'
)

genus <- "pergillus"

subset <- subset_samples(phylo, grepl(genus, Sample)) %>%
  filter_taxa(\(x) sum(x) > 0, prune = TRUE)

plot_bar(
  subset,
  fill='species',
  x = 'sample_Sample'
)

subset %>% otu_table() %>%
  cbind(
    tax_table(subset)
  ) %>%
  View()