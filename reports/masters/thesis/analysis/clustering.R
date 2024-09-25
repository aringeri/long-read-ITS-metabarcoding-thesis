library(phyloseq)
library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(scales)

source('./helpers/config.R')
source('./helpers/vsearch.R')


output_dir <- config$experiment_path
expected_species <- 55
n_samples <- 58

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

taxa_counts <- load_otus(c(0, 0.01, 0.005, 0.0015, 0.0005, 0.0001), output_dir)
expected_otus <- data.frame(yintercept=expected_species, expected=factor(expected_species))

otu_counts_by_sample_depth <- taxa_counts %>%
  mutate_at(vars(sample_depth), \(d) d*n_samples) %>%
  # mutate_at(vars(thresh), factor) %>%
  # mutate_at(vars(nOTUs), log) %>%
  ggplot(
    aes(x=sample_depth, y=nOTUs, colour=factor(thresh))
  ) +
  labs(
    title = "VSEARCH",
    x=glue("Total library size ({n_samples} isolates)"),
    y="Number of OTUs",
    colour = "Min OTU size\n(proportion of library size)",
    shape = "Reads per isolate"
  ) +
  # stat_summary(aes(x=factor(sample_depth), y = nOTUs, group=thresh, linetype=thresh), fun=mean, geom="line") +
  stat_summary(aes(x=sample_depth, y = nOTUs, group=thresh), fun=mean, geom="line") +
  facet_grid(. ~ thresh,
             labeller = as_labeller(\(x) label_percent()(as.numeric(x)))) +
  # geom_boxplot(aes(x=factor(sample_depth), y=nOTUs, fill=thresh), position = 'identity') +
  geom_point(aes(x=sample_depth, y=nOTUs, shape=factor(sample_depth)), size=1) +
  # geom_hline(yintercept = expected_species, linetype='dashed') +
  geom_hline(aes(yintercept = expected_species, linetype=expected), data=expected_otus, show.legend =TRUE) +
  # scale_colour_discrete(labels=paste0(c('20', '50', '167', '1000', '2000', '2500'), ' per sample'))+
  scale_colour_discrete(labels=function(x) {label_percent()(as.numeric(x))}) +
  scale_shape_discrete(labels = function(x) { as.numeric(x)/n_samples} ) +
  scale_x_continuous(transform = 'log10',
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_continuous(transform = 'log10', breaks = c(40, expected_species, 100, 200, 400, 800, 1600)) +
  coord_cartesian(ylim=c(40, 400)) +
  scale_linetype_manual(name='Expected number of OTUs', values=c("dashed"), drop=FALSE) +
  guides(colour = guide_legend(override.aes = list(linetype = 0)),
         shape = guide_legend(override.aes = list(linetype = 0))) +
  theme(aspect.ratio = 2)
  # annotate("text", x=log10(10^6), y=61, label="Expected number of OTUs (expected_species)", size=3, color="black")
otu_counts_by_sample_depth
ggsave('images/06-otu-count-vsearch.png', otu_counts_by_sample_depth)

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
output_dir <- '../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14'
nc <- load_nanoclust_stats(c(0, 0.0012, 0.006449999999999999, 0.0050999999999999995, 0.0072, 0.01005), output_dir)

otu_count_nanoclust <- load_nanoclust_stats(c(0, 0.0012, 0.006449999999999999, 0.0050999999999999995, 0.0072, 0.01005), output_dir) %>%
  mutate_at(vars(sample_depth), \(d) d*n_samples) %>%
  ggplot(
    aes(x=sample_depth, y=numOTUs, colour=factor(thresh))
  ) +
  labs(
    title = "UMAP + HDBSCAN",
    x=glue("Total library size ({n_samples} isolates)"),
    y="Number of OTUs",
    colour = "Minimum cluster size\n(proportion of library size)",
    shape = "Reads per isolate"
  ) +
  # stat_summary(aes(x=factor(sample_depth), y = numOTUs, group=thresh, linetype=thresh), fun=mean, geom="line") +
  stat_summary(aes(x=sample_depth, y = numOTUs, group=thresh), fun=mean, geom="line") +
  facet_grid(. ~ thresh,
    labeller = as_labeller(\(x)
                             ifelse(x == 0, "minimum cluster size = 2", label_percent()(as.numeric(x))))
  ) +
  # geom_boxplot(aes(x=factor(sample_depth), y=numOTUs, fill=thresh), position = 'identity') +
  geom_point(aes(x=sample_depth, y=numOTUs, shape=factor(sample_depth)), size=1) +
  # geom_hline(yintercept = 59, linetype='dashed') +
  geom_hline(aes(yintercept = expected_species, linetype=expected), data=expected_otus, show.legend =TRUE) +
  # scale_colour_discrete(labels=paste0(c('20', '50', '167', '1000', '2000', '2500'), ' per sample'))+
  scale_colour_discrete(labels=function(x) {label_percent()(as.numeric(x))}) +
  scale_shape_discrete(labels = function(x) { as.numeric(x)/n_samples} ) +
  scale_x_continuous(transform = 'log10',
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_continuous(transform = 'log10', breaks = c(50, expected_species, seq(60, 100, 10))) +
  scale_linetype_manual(name='Expected number of OTUs', values="dashed", drop=FALSE) +
  guides(colour = guide_legend(override.aes = list(linetype = 0)),
         shape = guide_legend(override.aes = list(linetype = 0))) +
  coord_cartesian(ylim=c(48, 100)) +
  theme(aspect.ratio = 2)
otu_count_nanoclust
ggsave('./images/06-otu-count-nanoclust.png', otu_count_nanoclust)

