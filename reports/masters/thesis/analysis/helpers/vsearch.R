filter_taxa_by_thresh <- function(phylo, thresh) {
  t <- sum(sample_sums(phylo))*thresh
  filter_taxa(phylo, \(x) sum(x) > t, prune=TRUE)
}

filter_taxa_by_min_otu_size <- function(phylo, min_otu_size) {
  filter_taxa(phylo, \(x) sum(x) >= min_otu_size, prune = TRUE)
}

filter_otu_by_sample <- function(phylo, thresh) {
  transform_sample_counts(phylo, function(col) ifelse(col/sum(col) < thresh, 0, col) ) %>%
    filter_taxa(\(x) sum(x) > 0, prune = TRUE)
}