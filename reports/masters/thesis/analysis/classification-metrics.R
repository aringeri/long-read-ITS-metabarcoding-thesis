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
source('helpers/config.R')

samplesheet <- read_samplesheet(config)

# nanoclust <- load_nanoclust_phyloseq(samplesheet, config$experiment_path, config$sample_depth, config$repetition)

calc_precision <- function(phylo, samplesheet) {
  tax_with_counts <- cbind.data.frame(
      tax_table(phylo),
      otu_table(phylo)
    ) %>%
    pivot_longer(barcode25:barcode89, names_to='barcode', values_to = 'count')

  classification_counts <- tax_with_counts %>%
    group_by(barcode) %>%
    filter(count > 0) %>%
    summarise(
      genus_classified = sum(count[genus != 'unidentified']),
      species_classified = sum(count[species != 'unidentified']),
      total = sum(count)
    )

  samplesheet_no_rownames <- samplesheet %>% tibble::rownames_to_column('barcode')

  precision_counts <- tax_with_counts %>%
    # group_by(barcode, genus) %>%
    filter(count > 0) %>%
    left_join(samplesheet_no_rownames, join_by(barcode), suffix = c('.actual', '.expected')) %>%
    summarise(
      # .by = c(barcode, genus.actual, genus.expected),
      # total = sum(count)
      .by = barcode,
      genus_correct = sum((genus_unite == genus)*count),
      species_correct = sum((species_unite == species.actual)*count),
      total = sum(count)
    ) %>%
    replace_na(list(species_correct = 0)) # NAs coming through for genus only samples

  left_join(classification_counts, precision_counts, join_by(barcode), suffix = c('', '.y'))
}

# precision_stats <- calc_precision(nanoclust, samplesheet)

# sample level
# cowplot::plot_grid(ncol = 1,
# precision_stats %>%
#   mutate(proportion = genus_classified / total) %>%
#   ggplot(aes(x=barcode, y=proportion)) +
#     geom_col() +
#     geom_hline(aes(yintercept = sum(genus_classified) / sum(total))) +
#     scale_x_discrete(
#       labels = \(x) sample_data(nanoclust)[x]$UpdatedName,
#       limits = rownames(sample_data(nanoclust))[order(sample_data(nanoclust)$UpdatedName)]
#     ) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
#   precision_stats %>%
#     mutate(precision = genus_correct / genus_classified) %>%
#     ggplot(aes(x=barcode, y=precision)) +
#     geom_col() +
#     geom_hline(aes(yintercept = sum(genus_correct) / sum(genus_classified))) +
#     scale_x_discrete(
#       labels = \(x) sample_data(nanoclust)[x]$UpdatedName,
#       limits = rownames(sample_data(nanoclust))[order(sample_data(nanoclust)$UpdatedName)]
#     ) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# )

load_all_precision_data <- function(config, samplesheet, load_phyloseq=load_nanoclust_phyloseq_3) {
  cluster_dir <- glue('{config$experiment_path}/hdbscan_clustering/FULL_ITS/')

  sample_depths <- list.dirs(
    cluster_dir,
    recursive = FALSE, full.names = FALSE) %>%
    as.numeric() %>%
    sort()

  df <- tibble()

  for (depth in sample_depths) {
    reps <- list.dirs(
      glue('{cluster_dir}/{depth}/'),
      recursive = FALSE, full.names = FALSE
    )
    for (rep in reps) {
      for (min_cluster_size in list.dirs(glue('{cluster_dir}/{depth}/{rep}'),recursive = FALSE, full.names = FALSE)) {
        print(glue("loading: {depth} {rep} {min_cluster_size}"))
        phylo <- load_phyloseq(samplesheet, config$experiment_path, depth, rep, min_cluster_size)
        precision_data <- calc_precision(phylo, samplesheet)
        precision_data$sample_depth <- depth
        precision_data$rep <- rep
        precision_data$min_cluster_size <- min_cluster_size
        df <- rbind(df, precision_data)
      }
    }
  }
  df
}


# all_df <- load_all_precision_data(config, samplesheet, load_nanoclust_phyloseq)
# all_df$method <- "nanoclust"
# all_df_vsearch <- load_all_precision_data(config, samplesheet, load_vsearch_phyloseq)
# all_df_vsearch$method <- "vsearch_abundant"

# cons_config <- Config$new(
#   experiment_path = "../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14-NC",
#   samplesheet_path = config$samplesheet_path,
#   sample_depth = config$sample_depth,
#   repetition = config$repetition
# )


all_df_cons <- load_all_precision_data(
  config, samplesheet,
  function(samplesheet, experiment, reads_per_sample, repetition, min_cluster_size) {
    load_nanoclust_phyloseq_3(samplesheet,
                            experiment = experiment,
                            sequence_type = 'nanoclust_consensus',
                            reads_per_sample = reads_per_sample,
                            repetition = repetition,
                            min_cluster_size = min_cluster_size
    )
  })
all_df_cons$method <- "nanoclust_consensus"

all_df <- load_all_precision_data(
  config, samplesheet,
  function(samplesheet, experiment, reads_per_sample, repetition, min_cluster_size) {
    load_nanoclust_phyloseq_3(samplesheet,
                              experiment = experiment,
                              sequence_type = 'nanoclust_abundant',
                              reads_per_sample = reads_per_sample,
                              repetition = repetition,
                              min_cluster_size = min_cluster_size
    )
  })
all_df$method <- "nanoclust_abundant"

precision_summary <- rbind(all_df, all_df_cons) %>%
  # filter(sample_depth %in% c(1000, 2000)) %>%
  filter(sample_depth > 160) %>%
  summarise(
    .by = c(method, sample_depth, rep, min_cluster_size),
    genus_classification_prop = sum(genus_classified) / sum(total),
    genus_classification_prop_by_depth = sum(genus_classified) / sum(sample_depth),
    species_classification_prop = sum(species_classified) / sum(total),
    species_classification_prop_by_depth = sum(species_classified) / sum(sample_depth),
    genus_precision = sum(genus_correct) / sum(genus_classified),
    species_precision = sum(species_correct) / sum(species_classified),
  ) #%>%
  # summarise(
  #   .by = c(method, sample_depth),
  #   genus_classification_prop = mean(genus_classification_prop),
  #   species_classification_prop = mean(species_classification_prop),
  #   genus_precision = mean(genus_precision),
  #   species_precision = mean(species_precision)
  # )

precision_plot <- cowplot::plot_grid(
  precision_summary %>%
    ggplot(aes(x=genus_precision, y=genus_classification_prop_by_depth, shape=factor(sample_depth), colour=factor(method))) +
      geom_point() +
      expand_limits(x=c(.75, .825), y=c(.8, 1)) +
      facet_wrap(~sample_depth + min_cluster_size) +
      scale_y_continuous(labels = scales::percent) +
      scale_x_continuous(labels = scales::percent) +
      scale_shape_discrete(name="reads per sample") +
      labs(x="Genera precision (%)", y="Genera classification proportion (%)") +
      theme(aspect.ratio=1)#,
  # precision_summary %>%
  #   ggplot(aes(x=species_precision, y=species_classification_prop, shape=factor(sample_depth), colour=factor(min_cluster_size))) +
  #   geom_point() +
  #   facet_wrap(~method) +
  #   expand_limits(x=c(.625, .825), y=c(.6, 1)) +
  #   scale_y_continuous(labels = scales::percent) +
  #   scale_x_continuous(labels = scales::percent) +
  #   scale_shape_discrete(name="reads per sample") +
  #   labs(x="Species precision (%)", y="Species classification proportion (%)") +
  #   theme(aspect.ratio=1)
)
precision_plot
ggsave('./images/06-precision-nanoclust-abundance.png', precision_plot)