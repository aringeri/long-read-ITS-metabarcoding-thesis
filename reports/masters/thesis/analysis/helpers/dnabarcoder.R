read_dna_barcoder_classification_nanoclust <- function(filepath) {
  tax <- read.csv(filepath, sep='\t', row.names = 1) %>%
    tibble::rownames_to_column('OTU') %>%
    separate(OTU, into=c('read', 'barcode', 'cluster', NA, 'size'), sep=';') %>%
    mutate_at(vars(read, barcode, cluster, size), \(str) gsub('.*=', '', str)) %>%
    filter(cluster != -1) %>%
    arrange(as.numeric(cluster)) %>%
    tibble::column_to_rownames("cluster")

  tax_table(as.matrix(tax[, !(names(tax) %in% c('read', 'barcode', 'size', 'ReferenceID', 'rank', 'score', 'cutoff', 'confidence'))]))
}

read_dna_barcoder_classification_vsearch <- function(filepath) {
  tax <- read.csv(filepath, sep='\t', row.names = 1) %>%
    tibble::rownames_to_column('OTU') %>%
    separate(OTU, into=c('read', 'barcode', NA, 'size'), sep=';') %>%
    mutate_at(vars(read, barcode, size), \(str) gsub('.*=', '', str)) %>%
    tibble::column_to_rownames('read')

  tax_table(as.matrix(tax[, !(names(tax) %in% c('read', 'barcode', 'size', 'ReferenceID', 'rank', 'score', 'cutoff', 'confidence'))]))
}

load_nanoclust_phyloseq <- function(
  samplesheet,
  experiment="../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14",
  reads_per_sample=2000,
  repetition=2)
{
  otu <- read.csv(glue('{experiment}/hdbscan_clustering/FULL_ITS/{reads_per_sample}/{repetition}/otu_table/otu_table.tsv'), sep = '\t', row.names = 1)
  otu <- otu[order(as.numeric(rownames(otu))), ]
  otu <- otu[rownames(otu) != -1, ]

  tax <- read_dna_barcoder_classification_nanoclust(glue('{experiment}/dnabarcoder/nanoclust/FULL_ITS/{reads_per_sample}/{repetition}/classify/all_samples.unite2024ITS_BLAST.classification'))

  phylo <- phyloseq(
    otu_table(otu, taxa_are_rows = TRUE),
    tax,
    sample_data(samplesheet)
  )

  phylo
  # tax_fix(phylo, min_length = 0, unknowns = 'unidentified', anon_unique = F)
}