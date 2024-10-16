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

source('helpers/dnabarcoder.R')
source('helpers/vsearch.R')
source('helpers/config.R')

soil_full_dataset <- '../../../../experiments/camille-soil-samples/outputs/c-soil-singleton-10-11/'
experiment <- soil_full_dataset

otus <- readRDS(glue('{experiment}/phyloseq/FULL_ITS/ALL_READS/1/all_samples/all_samples.phyloseq.rds'))
classFile <- list.files(
  glue('{experiment}/dnabarcoder/vsearch/FULL_ITS/ALL_READS/1/classify/'),
  pattern='*.unite2024ITS_BLAST.classification',
  full.names = T
)

tax_table <- read_dna_barcoder_classification_vsearch(classFile[1])
soil_full_phylo <- phyloseq(otu_table(otus), tax_table)

filt <- soil_full_phylo #%>% filter_taxa_by_min_otu_size(700)
filt

# TODO filter by sample 0.005
# ?prune_samples()
# sample_sums(soil_full_phylo) * .005
ntaxa(soil_full_phylo)
estimate_richness(soil_full_phylo, measures = 'Observed')
soil_full_phylo

filt %>%
  filter_taxa_by_min_otu_size(50) %>%
  transform_sample_counts(\(x) x/sum(x)) %>%
  # prune_samples('sample_AL4', .) %>%
  tax_glom(taxrank="species") %>%
  plot_bar(fill="species == \"unidentified\"")
prune_taxa(rownames(tax_table(filt)[tax_table(filt)[,'genus'] != 'unidentified', ]), filt)


soil_subset_dataset <- '../../../../experiments/camille-soil-samples/outputs/c-soil-NC-10-13/'
subset_vsearch <- load_vsearch_phyloseq(NULL, experiment = soil_subset_dataset, 4000, 1, F)
subset_vsearch %>%
  filter_taxa_by_min_otu_size(50) %>%
  transform_sample_counts(\(x) x/sum(x)) %>%
  # prune_samples('sample_AL4', .) %>%
  tax_glom(taxrank="kingdom") %>%
  plot_bar(fill="kingdom == \"unidentified\"")

load_nc <- function(experiment, min_cluster_size) {
  otu <- read.csv(glue('{experiment}/hdbscan_clustering/FULL_ITS/4000/1/{min_cluster_size}/otu_table/otu_table.tsv'), sep = '\t', row.names = 1)
  otu <- otu[order(as.numeric(rownames(otu))), ]
  otu <- otu[rownames(otu) != -1, ]

  classFile <- list.files(
    glue('{experiment}/dnabarcoder/nanoclust_abundant/FULL_ITS/4000/1/{min_cluster_size}/classify/'),
    pattern='*.unite2024ITS_BLAST.classification',
    full.names = T
  )
  tax <- read.csv(classFile[1], sep='\t', row.names = 1) %>%
    tibble::rownames_to_column('OTU') %>%
    separate(OTU, into=c('read', 'cluster', 'barcode', NA, 'size'), sep=';') %>%
    mutate_at(vars(read, barcode, cluster, size), \(str) gsub('.*=', '', str)) %>%
    filter(cluster != -1) %>%
    arrange(as.numeric(cluster)) %>%
    tibble::column_to_rownames("cluster")
  tax <- tax_table(as.matrix(tax[, !(names(tax) %in% c('read', 'barcode', 'size', 'ReferenceID', 'rank', 'score', 'cutoff', 'confidence'))]))

  phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax)
}

subset_nc <- load_nc(soil_subset_dataset, 2)

fix <- subset_nc %>% tax_fix(unknowns = 'unidentified')
fix %>%
  transform_sample_counts(\(x) x/sum(x)) %>%
  # prune_samples('sample_AL4', .) %>%
  tax_glom(taxrank="genus") %>%
  plot_bar(fill="species %in% 1:700")

subset_nc %>%
  # transform_sample_counts(\(x) x/sum(x)) %>%
  # prune_samples('sample_AL4', .) %>%
  tax_glom(taxrank="species") %>%
  plot_bar(fill="paste0(kingdom == \"unidentified\", 1)")


tax_count <- cbind(tax_table(fix), taxa_sums(fix))
un_id <- tax_count[tax_count[,'species'] %in% 1:700, ]
un_id %>% as.data.frame() %>% mutate(count = as.numeric(V8)) %>% arrange(count) %>% summarise(x= sum(count))

load_nanoclust_phyloseq_3(NULL, experiment = soil_subset_dataset, sequence_type="nanoclust_abundant", 4000, 1, 2)

