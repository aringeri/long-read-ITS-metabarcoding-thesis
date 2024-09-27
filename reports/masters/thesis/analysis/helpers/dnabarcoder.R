load_dna_barcoder_classification_nanoclust <- function(filepath) {
  read.csv(filepath, sep='\t', row.names = 1) %>%
    tibble::rownames_to_column('OTU') %>%
    separate(OTU, into=c('read', 'barcode', 'cluster', NA, 'size'), sep=';') %>%
    mutate_at(vars(read, barcode, cluster, size), \(str) gsub('.*=', '', str)) %>%
    filter(cluster != -1) %>%
    arrange(as.numeric(cluster)) %>%
    tibble::column_to_rownames("cluster")
}

load_dna_barcoder_classification_nanoclust_2 <- function(
  experiment="../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-09-12",
  sequence_type="nanoclust_abundant",
  reads_per_sample=2000,
  repetition=2,
  min_cluster_size=0
) {
  classFile <- list.files(
    glue('{experiment}/dnabarcoder/{sequence_type}/FULL_ITS/{reads_per_sample}/{repetition}/{min_cluster_size}/classify/'),
    pattern='*.unite2024ITS_BLAST.classification',
    full.names = T
  )
  load_dna_barcoder_classification_nanoclust(classFile[1])
}

read_dna_barcoder_classification_nanoclust <- function(filepath) {
  tax <- load_dna_barcoder_classification_nanoclust(filepath)

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
  # print(experiment)
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

load_nanoclust_phyloseq_2 <- function(
  samplesheet,
  experiment="../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14",
  sequence_type="nanoclust_abundant",
  reads_per_sample=2000,
  repetition=2)
{
  # print(experiment)
  otu <- read.csv(glue('{experiment}/hdbscan_clustering/FULL_ITS/{reads_per_sample}/{repetition}/otu_table/otu_table.tsv'), sep = '\t', row.names = 1)
  otu <- otu[order(as.numeric(rownames(otu))), ]
  otu <- otu[rownames(otu) != -1, ]

  classFile <- list.files(
    glue('{experiment}/dnabarcoder/{sequence_type}/FULL_ITS/{reads_per_sample}/{repetition}/classify/'),
    pattern='*.unite2024ITS_BLAST.classification',
    full.names = T
  )
  tax <- read_dna_barcoder_classification_nanoclust(classFile[1])

  phylo <- phyloseq(
    otu_table(otu, taxa_are_rows = TRUE),
    tax,
    sample_data(samplesheet)
  )

  phylo
  # tax_fix(phylo, min_length = 0, unknowns = 'unidentified', anon_unique = F)
}

load_nanoclust_phyloseq_3 <- function(
  samplesheet,
  experiment="../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-09-12",
  sequence_type="nanoclust_abundant",
  reads_per_sample=2000,
  repetition=2,
  min_cluster_size=0)
{
  # print(experiment)
  otu <- read.csv(glue('{experiment}/hdbscan_clustering/FULL_ITS/{reads_per_sample}/{repetition}/{min_cluster_size}/otu_table/otu_table.tsv'), sep = '\t', row.names = 1)
  otu <- otu[order(as.numeric(rownames(otu))), ]
  otu <- otu[rownames(otu) != -1, ]

  classFile <- list.files(
    glue('{experiment}/dnabarcoder/{sequence_type}/FULL_ITS/{reads_per_sample}/{repetition}/{min_cluster_size}/classify/'),
    pattern='*.unite2024ITS_BLAST.classification',
    full.names = T
  )
  tax <- read_dna_barcoder_classification_nanoclust(classFile[1])

  phylo <- phyloseq(
    otu_table(otu, taxa_are_rows = TRUE),
    tax,
    sample_data(samplesheet)
  )

  phylo
  # tax_fix(phylo, min_length = 0, unknowns = 'unidentified', anon_unique = F)
}

load_vsearch_phyloseq <- function(
  samplesheet,
  experiment="../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14",
  reads_per_sample=2000,
  repetition=2,
  consensus=F)
{
  vsearch_otus <- readRDS(glue('{experiment}/phyloseq/FULL_ITS/{reads_per_sample}/{repetition}/all_samples/all_samples.phyloseq.rds'))
  rep_seq <- if (consensus) 'vsearch_consensus' else 'vsearch'
  classFile <- list.files(
    glue('{experiment}/dnabarcoder/{rep_seq}/FULL_ITS/{reads_per_sample}/{repetition}/classify/'),
    pattern='*.unite2024ITS_BLAST.classification',
    full.names = T
  )
  tax_table <- read_dna_barcoder_classification_vsearch(classFile[1])

  vsearch_phylo <- phyloseq(vsearch_otus, tax_table, sample_data(samplesheet))
}