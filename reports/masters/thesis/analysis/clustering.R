library(phyloseq)
library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)

problematic <- readRDS('../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-08/phyloseq/FULL_ITS/2500/1/all_samples/all_samples.phyloseq.rds')

# TODO update to latest (08-08)
# output_dir <- '../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-02'
output_dir <- '../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-08'

filter_taxa_by_thresh <- function(phylo, thresh) {
  t <- sum(sample_sums(phylo))*thresh
  filter_taxa(phylo, \(x) sum(x) > t, prune=TRUE)
}

load_otus <- function(filter_thresh = c(0, 0.006), output_dir) {
  stopifnot(length(filter_thresh) > 0)
  df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df) <- c("sample_depth"
                   ,"thresh"
                   ,'nOTUs')
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
      for (thresh in filter_thresh) {
        filt_otus <- filter_taxa_by_thresh(otus, thresh)
        # df[depth, paste0('rep', rep)] <- ntaxa(filt_otus)
        df <- df %>% add_row(sample_depth = depth, thresh = thresh, nOTUs = ntaxa(filt_otus))
      }
    }
  }
  return(df)
}

library(scales)

taxa_counts <- load_otus(c(0, 0.01, 0.005, 0.001, 0.0005, 0.0001), output_dir)
expected_otus <- data.frame(yintercept=59, expected=factor(59))

otu_counts_by_sample_depth <- taxa_counts %>%
  mutate_at(vars(sample_depth), \(d) d*60) %>%
  # mutate_at(vars(thresh), factor) %>%
  # mutate_at(vars(nOTUs), log) %>%
  ggplot(
    aes(x=sample_depth, y=nOTUs, colour=factor(thresh))
  ) +
  labs(
    title = "vsearch clustering",
    x="total library size (60 samples)",
    y="number of OTUs",
    colour = "Min OTU size\n(proportion of library size)",
    shape = "Reads per sample"
  ) +
  # stat_summary(aes(x=factor(sample_depth), y = nOTUs, group=thresh, linetype=thresh), fun=mean, geom="line") +
  stat_summary(aes(x=sample_depth, y = nOTUs, group=thresh), fun=mean, geom="line") +
  facet_grid(. ~ thresh,
             labeller = as_labeller(\(x) label_percent()(as.numeric(x)))) +
  # geom_boxplot(aes(x=factor(sample_depth), y=nOTUs, fill=thresh), position = 'identity') +
  geom_point(aes(x=sample_depth, y=nOTUs, shape=factor(sample_depth)), size=1) +
  # geom_hline(yintercept = 59, linetype='dashed') +
  geom_hline(aes(yintercept = 59, linetype=expected), data=expected_otus, show.legend =TRUE) +
  # scale_colour_discrete(labels=paste0(c('20', '50', '167', '1000', '2000', '2500'), ' per sample'))+
  scale_colour_discrete(labels=function(x) {label_percent()(as.numeric(x))}) +
  scale_shape_discrete(labels = function(x) { as.numeric(x)/60} ) +
  scale_x_continuous(transform = 'log10',
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_continuous(transform = 'log10', breaks = c(40, 59, 100, 200, 400, 800, 1600)) +
  scale_linetype_manual(name='Expected number of OTUs', values=c("dashed"), drop=FALSE) +
  guides(colour = guide_legend(override.aes = list(linetype = 0)),
         shape = guide_legend(override.aes = list(linetype = 0)))
  # annotate("text", x=log10(10^6), y=61, label="Expected number of OTUs (59)", size=3, color="black")
# otu_counts_by_sample_depth
ggsave('images/06-otu-count-vsearch.png',otu_counts_by_sample_depth)

# thresh_range <- seq(0,0.02, length.out = 300)
# vsearch_otus_1000 <- readRDS(glue('../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-02/phyloseq/FULL_ITS/2500/2/all_samples/all_samples.phyloseq.rds'))
#
# nanoclust_1000 <- read.csv('../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-02/hdbscan_clustering_stats/FULL_ITS/2500/2/min_cluster_size_stats.tsv', sep='\t')
# nanoclust_1000 <- nanoclust_1000[, c('thresh', 'loss', 'numOTUs')]
# nanoclust_1000 <- nanoclust_1000[nanoclust_1000[,'thresh'] <= 0.02,]
# colnames(nanoclust_1000) <- c('thresh', 'loss', 'nOTUs')
# nanoclust_1000$method <- 'nanoclust'
#
# vsearch_stats <- data.frame(thresh = thresh_range) %>%
#   rowwise() %>%
#   mutate(
#     nOTUs = ntaxa(filter_taxa_by_thresh(vsearch_otus_1000, thresh)),
#     loss = 1 - sum(sample_sums(filter_taxa_by_thresh(vsearch_otus_1000, thresh)))/sum(sample_sums(vsearch_otus_1000))
#   )
# vsearch_stats$method <- 'vsearch'
# combined <- rbind(vsearch_stats, nanoclust_1000)
# n_otus_by_thresh <- combined %>%
#   ggplot(
#     aes(x=thresh, colour = method)
#   ) +
#   geom_point(aes(y=nOTUs), size=0.5) +
#   stat_summary(aes(y = nOTUs), fun=mean, geom="line") +
#   # geom_point(aes(y=loss*500, colour=method), size=.5) +
#   scale_y_continuous(
#     "number of OTUs",
#     # transform = 'log10',
#     # sec.axis = sec_axis(~ . / 500, name = "loss", breaks = seq(0,1,.05)),
#     minor_breaks = \(x) seq(0, max(x), 10),
#     breaks = \(x) seq(0, max(x), 50)
#   ) +
#   scale_x_continuous(
#     minor_breaks = \(x) seq(0, max(x), 0.0025/2),
#     breaks = \(x) seq(0, max(x), 0.0025),
#     labels = \(x) label_percent()(as.numeric(x)),
#   ) +
#   labs(x="minimum OTU size threshold (proportion of library size)") +
#   geom_hline(yintercept = 59, linetype='dashed')
#
# loss_by_thresh <- combined %>%
#   ggplot(
#     aes(x=thresh, colour=method)
#   ) +
#   geom_point(aes(y=loss), size=0.5) +
#   stat_summary(aes(y = loss), fun=mean, geom="line") +
#   scale_y_continuous(
#     "read loss",
#     labels = \(y) label_percent()(as.numeric(y)),
#     breaks = \(y) seq(0, max(y), 0.05)
#   ) +
#   scale_x_continuous(
#     minor_breaks = \(x) seq(0, max(x), 0.0025/2),
#     breaks = \(x) seq(0, max(x), 0.0025),
#     labels = \(x) label_percent()(as.numeric(x)),
#   ) +
#   labs(
#     x="minimum OTU size threshold (proportion of library size)",
#   )
#
# cowplot::plot_grid(n_otus_by_thresh, loss_by_thresh)



load_nanoclust_stats <- function(filter_thresh = c(0, 0.0012), output_dir) {
  stopifnot(length(filter_thresh) > 0)
  df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df) <- c("sample_depth"
    ,"thresh"
    ,'nOTUs')
  sample_depths <- list.dirs(
    glue('{output_dir}/hdbscan_clustering_stats/FULL_ITS/'),
    recursive = FALSE, full.names = FALSE) %>%
    as.numeric() %>%
    sort()

  for (depth in sample_depths) {
    reps <- list.dirs(
      glue('{output_dir}/hdbscan_clustering_stats/FULL_ITS/{depth}/'),
      recursive = FALSE, full.names = FALSE
    )
    stopifnot(length(reps) > 0)

    for (rep in reps) {
      stats <- read.csv(glue('{output_dir}/hdbscan_clustering_stats/FULL_ITS/{depth}/{rep}/min_cluster_size_stats.tsv'), sep = '\t', stringsAsFactors = FALSE)

      stats_in_filt <- stats[stats[, 'thresh'] %in% filter_thresh & !(stats[, 'thresh'] != 0 & stats[,'min_cluster_size'] == 2), ]
      stats_in_filt$sample_depth <- depth
      stats_in_filt$rep <- rep
      # print(stats_in_filt)
      df <- rbind(df, stats_in_filt)
      # for (thresh in filter_thresh) {
      #   filt_otus <- filter_taxa_by_thresh(otus, thresh)
      #   # df[depth, paste0('rep', rep)] <- ntaxa(filt_otus)
      #   df <- df %>% add_row(sample_depth = depth, thresh = thresh, nOTUs = ntaxa(filt_otus))
      # }
    }
  }
  return(df)
}
#load_nanoclust_stats(c(0, 0.0006)) %>% View()

otu_count_nanoclust <- load_nanoclust_stats(c(0, 0.0006, 0.0012, 0.0050999999999999995, 0.0072, 0.01005), output_dir) %>%
  mutate_at(vars(sample_depth), \(d) d*60) %>%
  ggplot(
    aes(x=sample_depth, y=numOTUs, colour=factor(thresh))
  ) +
  labs(
    title = "NanoCLUST",
    x="total library size (60 samples)",
    y="number of clusters",
    colour = "Minimum cluster size\n(proportion of library size)",
    shape = "Reads per sample"
  ) +
  # stat_summary(aes(x=factor(sample_depth), y = numOTUs, group=thresh, linetype=thresh), fun=mean, geom="line") +
  stat_summary(aes(x=sample_depth, y = numOTUs, group=thresh), fun=mean, geom="line") +
  facet_grid(. ~ thresh,
             labeller = as_labeller(\(x) label_percent()(as.numeric(x)))) +
  # geom_boxplot(aes(x=factor(sample_depth), y=numOTUs, fill=thresh), position = 'identity') +
  geom_point(aes(x=sample_depth, y=numOTUs, shape=factor(sample_depth)), size=1) +
  # geom_hline(yintercept = 59, linetype='dashed') +
  geom_hline(aes(yintercept = 59, linetype=expected), data=expected_otus, show.legend =TRUE) +
  # scale_colour_discrete(labels=paste0(c('20', '50', '167', '1000', '2000', '2500'), ' per sample'))+
  scale_colour_discrete(labels=function(x) {label_percent()(as.numeric(x))}) +
  scale_shape_discrete(labels = function(x) { as.numeric(x)/60} ) +
  scale_x_continuous(transform = 'log10',
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_continuous(transform = 'log10', breaks = c(50, 59, 100, 200, 400, 800, 1600)) +
  scale_linetype_manual(name='Expected number of OTUs', values=c("dashed"), drop=FALSE) +
  guides(colour = guide_legend(override.aes = list(linetype = 0)),
         shape = guide_legend(override.aes = list(linetype = 0)))
# otu_count_nanoclust
ggsave('./images/06-otu-count-nanoclust.png', otu_count_nanoclust)

