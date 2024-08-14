library(phyloseq)
library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(scales)

filter_taxa_by_thresh <- function(phylo, thresh) {
  t <- sum(sample_sums(phylo))*thresh
  filter_taxa(phylo, \(x) sum(x) > t, prune=TRUE)
}

output_dir <- '../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-08'

load_loss_stats <- function(output_dir = '../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-02') {
  thresh_range <- seq(0,0.02, length.out = 300)

  df <- data.frame(matrix(ncol = 6, nrow = 0))
  # colnames(df) <- c("sample_depth" ,"thresh" ,'nOTUs')
  sample_depths <- list.dirs(
    glue('{output_dir}/phyloseq/FULL_ITS'),
    recursive = FALSE, full.names = FALSE) %>%
    as.numeric() %>%
    sort()

  for (depth in sample_depths) {
    reps <- list.dirs(
      glue('{output_dir}/phyloseq/FULL_ITS/{depth}/'),
      recursive = FALSE, full.names = FALSE
    )
    stopifnot(length(reps) > 0)


    # df[depth, 'thresh'] <- thresh
    # print(glue('{depth} / {reps}'))
    for (rep in reps) {
      otus <- readRDS(glue('{output_dir}/phyloseq/FULL_ITS/{depth}/{rep}/all_samples/all_samples.phyloseq.rds'))
      vsearch_stats <- data.frame(thresh = thresh_range) %>%
        rowwise() %>%
        mutate(
          nOTUs = ntaxa(filter_taxa_by_thresh(otus, thresh)),
          loss = 1 - sum(sample_sums(filter_taxa_by_thresh(otus, thresh)))/sum(sample_sums(otus))
        )
      vsearch_stats$method <- 'vsearch'
      vsearch_stats$rep <- rep
      vsearch_stats$sample_depth <- depth

      nanoclust_stats <- read.csv(glue('{output_dir}/hdbscan_clustering_stats/FULL_ITS/{depth}/{rep}/min_cluster_size_stats.tsv'), sep='\t')
      nanoclust_stats <- nanoclust_stats[, c('thresh', 'loss', 'numOTUs')]
      nanoclust_stats <- nanoclust_stats[nanoclust_stats[, 'thresh'] <= 0.02,]
      colnames(nanoclust_stats) <- c('thresh', 'loss', 'nOTUs')
      nanoclust_stats$method <- 'nanoclust'
      nanoclust_stats$rep <- rep
      nanoclust_stats$sample_depth <- depth

      df <- rbind(df, rbind(vsearch_stats, nanoclust_stats))
    }
  }
  return(df)
}
combined_stats <- load_loss_stats(output_dir)

otu_counts <- combined_stats %>%
  ggplot(
    aes(x=thresh, y=nOTUs, colour=method)
  ) +
  geom_point(size=0.2) +
  stat_summary(fun=mean, geom="line") +
  geom_hline(aes(yintercept = 59, linetype='actual'))+#, show.legend =TRUE) +
  geom_vline(aes(colour='vsearch'), xintercept = 0.001) +
  geom_vline(aes(colour='nanoplot'), xintercept = 0.005) +
  facet_grid(. ~ sample_depth, labeller = as_labeller(\(x) paste0(60 * as.numeric(x), " (", paste0(x, " reads per sample"), ")"))) +
  scale_y_continuous(transform = 'log10', breaks = c(0,10,20, 40, 59, 100, 200, 400, 800, 1600, 2000)) +
  scale_x_continuous(labels = \(x) label_percent()(as.numeric(x))) +
  labs(
    x="Minimum cluster size threshold (proportion of library size)",
    y="Number of clusters"
  ) +
  scale_linetype_manual(name='', values=c("dashed"), drop=FALSE)

# K-mer size = 6
otu_loss <- combined_stats %>%
  ggplot(
    aes(x=thresh, y=loss, colour=method)
  ) +
  geom_point(size=0.2) +
  stat_summary(fun=mean, geom="line") +
  geom_vline(aes(colour='vsearch'), xintercept = 0.001) +
  geom_vline(aes(colour='nanoplot'), xintercept = 0.005) +
  facet_grid(. ~ sample_depth, labeller = as_labeller(\(x) paste0(60 * as.numeric(x), " (", paste0(x, " reads per sample"), ")"))) +
  scale_y_continuous(breaks=\(y) seq(0, max(y), 0.05), labels = \(y) label_percent()(as.numeric(y))) +
  scale_x_continuous(labels = \(x) label_percent()(as.numeric(x))) +
  labs(
    x="Minimum cluster size threshold (proportion of library size)",
    y="Read loss"
  )

combined_plots <- cowplot::plot_grid(otu_counts, otu_loss, ncol=1)

ggsave('./images/06-otu-count-and-loss.png', combined_plots)


# K-mer size = 5
if (FALSE) {
  k5_stats <- load_loss_stats('../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-k5-08-07')

  k5_plot <- cowplot::plot_grid(
    k5_stats %>%
      ggplot(
        aes(x=thresh, y=nOTUs, colour=method)
      ) +
      geom_point(size=0.2) +
      stat_summary(fun=mean, geom="line") +
      geom_hline(aes(yintercept = 59, linetype='actual'))+#, show.legend =TRUE) +
      geom_vline(aes(colour='vsearch'), xintercept = 0.001) +
      geom_vline(aes(colour='nanoplot'), xintercept = 0.005) +
      facet_grid(. ~ sample_depth, labeller = as_labeller(\(x) paste0(60 * as.numeric(x), " (", paste0(x, " reads per sample"), ")"))) +
      scale_y_continuous(transform = 'log10', breaks = c(0,10,20, 40, 59, 100, 200, 400, 800, 1600, 2000)) +
      scale_x_continuous(labels = \(x) label_percent()(as.numeric(x))) +
      labs(
        x="Minimum cluster size threshold (proportion of library size)",
        y="Number of clusters"
      ) +
      scale_linetype_manual(name='', values=c("dashed"), drop=FALSE),
    k5_stats %>%
      ggplot(
        aes(x=thresh, y=loss, colour=method)
      ) +
      geom_point(size=0.2) +
      stat_summary(fun=mean, geom="line") +
      geom_vline(aes(colour='vsearch'), xintercept = 0.001) +
      geom_vline(aes(colour='nanoplot'), xintercept = 0.005) +
      facet_grid(. ~ sample_depth, labeller = as_labeller(\(x) paste0(60 * as.numeric(x), " (", paste0(x, " reads per sample"), ")"))) +
      scale_y_continuous(breaks=\(y) seq(0, max(y), 0.05), labels = \(y) label_percent()(as.numeric(y))) +
      scale_x_continuous(labels = \(x) label_percent()(as.numeric(x))) +
      labs(
        x="Minimum cluster size threshold (proportion of library size)",
        y="Read loss"
      ),
    ncol=1
  )
  ggsave('./images/06-otu-count-and-loss-k5.png', k5_plot)
}


