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
source('helpers/vsearch.R')
source('helpers/manual-taxonomy.R')

samplesheet <- read_samplesheet(config)
nsamples <- 58

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
      genus_correct = sum((genus_unite == genus | (genus %in% accepted_synonyms$alt_genus & species_unite %in% accepted_synonyms$actual_species))*count),
      species_correct = sum((species_unite == species.actual | (species.actual %in% accepted_synonyms$alt_species & species_unite %in% accepted_synonyms$actual_species))*count),
      total = sum(count)
    ) %>%
    replace_na(list(species_correct = 0)) # NAs coming through for genus only samples

  left_join(classification_counts, precision_counts, join_by(barcode), suffix = c('', '.y'))
}

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

selected_min_cluster_sizes <- c(0, 0.001, 0.005)

precision_summary <- rbind(all_df, all_df_cons) %>%
  # filter(sample_depth %in% c(1000, 2000)) %>%
  filter(min_cluster_size %in% c(0, 0.001, 0.005)) %>%
  filter(sample_depth > 49) %>%
  mutate(
    library_size = sample_depth*nsamples,
    min_cluster_size = as.numeric(min_cluster_size),
    # method = ifelse(method == 'nanoclust_abundant', "most abundant", "consensus")
  ) %>%
  summarise(
    .by = c(method, library_size, rep, min_cluster_size),
    genus_classification_prop = sum(genus_classified) / sum(total),
    genus_classification_prop_by_depth = sum(genus_classified) / sum(sample_depth),
    species_classification_prop = sum(species_classified) / sum(total),
    species_classification_prop_by_depth = sum(species_classified) / sum(sample_depth),
    genus_precision = sum(genus_correct) / sum(genus_classified),
    species_precision = sum(species_correct) / sum(species_classified),
  )

precision_plot <-
  precision_summary %>%
    ggplot(aes(x=genus_precision, y=genus_classification_prop_by_depth,
               shape=factor(library_size), colour=factor(method, labels = c("most abundant", "consensus")))) +
      geom_point() +
      # expand_limits(x=c(.75, .825), y=c(.82, 1)) +
      expand_limits(y=c(.82, 1)) +
      facet_grid(
        rows=vars(min_cluster_size), cols=vars(method),
        labeller = as_labeller(\(x)
           ifelse(is.na(as.numeric(x)), ifelse(x == 'nanoclust_abundant', "most abundant", "consensus") , paste0("min OTU size (", scales::percent(as.numeric(x)), ")"))
        ),
      ) +
      scale_y_continuous(labels = scales::percent) +
      scale_x_continuous(labels = scales::percent, n.breaks = 3) +
      scale_shape_discrete(name="Library size") +
      scale_colour_discrete(name="Representative sequence") +
      labs( x="Genera precision (%)", y="Genera classification proportion (%)") +
      theme(aspect.ratio=1)

precision_plot
# ggsave('./images/06-precision-nanoclust-abundance.png', precision_plot)


vsearch_stats <- tibble()
for (sample_depth in c(20, 50, 167, 1000, 2000, 2500)) {
  for (t in c(0, 0.0005, 0.001, 0.0015)) {
    for (rep in 1:5) {
      # if (sample_depth==2000) {
        vsearch_cons <- load_vsearch_phyloseq(samplesheet, config$experiment_path, 2000, rep, consensus=T)
        vsearch_stats_cons <- calc_precision(vsearch_cons  %>% filter_taxa_by_thresh(t), samplesheet)
        vsearch_stats_cons$method <- 'consensus'
        vsearch_stats_cons$sample_depth <- 2000
        vsearch_stats_cons$rep <- rep
        vsearch_stats_cons$min_cluster_size <- t
      # }


      vsearch_abund <- load_vsearch_phyloseq(samplesheet, config$experiment_path, sample_depth, rep, consensus=F)
      vsearch_stats_abund <- calc_precision(vsearch_abund %>% filter_taxa_by_thresh(t), samplesheet)
      vsearch_stats_abund$method <- 'abundance'
      vsearch_stats_abund$sample_depth <- sample_depth
      vsearch_stats_abund$rep <- rep
      vsearch_stats_abund$min_cluster_size <- t

      vsearch_stats <- rbind(vsearch_stats, vsearch_stats_abund, vsearch_stats_cons)
    }
  }
}

vsearch_precision_summary <- vsearch_stats %>%
  # filter(sample_depth %in% c(1000, 2000)) %>%
  filter(sample_depth != 20) %>%
  filter(min_cluster_size %in% c(0, 0.0015)) %>%
  filter(method %in% c('abundance', 'consensus')) %>%
  # filter(sample_depth > 49) %>%
  mutate(library_size = sample_depth * nsamples) %>%
  summarise(
    .by = c(method, sample_depth, rep, min_cluster_size, library_size),
    genus_classification_prop = sum(genus_classified) / sum(total),
    genus_classification_prop_by_depth = sum(genus_classified) / sum(sample_depth),
    species_classification_prop = sum(species_classified) / sum(total),
    species_classification_prop_by_depth = sum(species_classified) / sum(sample_depth),
    genus_precision = sum(genus_correct) / sum(genus_classified),
    species_precision = sum(species_correct) / sum(species_classified),
  )

vsearch_abundance_precision <- vsearch_precision_summary %>%
  filter(method == 'abundance') %>%
  ggplot(aes(x=genus_precision, y=genus_classification_prop_by_depth, shape=factor(library_size))) +
  geom_point(aes(colour='most abundant')) + #colour='#35BBC0') +
  expand_limits(x=c(.76, .85), y=c(.82, 1)) +
  facet_grid(rows=vars(min_cluster_size), cols=vars(method),
    labeller = as_labeller(\(x)
      ifelse(is.na(as.numeric(x)), "most abundant", paste0("min OTU size (", scales::percent(as.numeric(x)), ")"))
    )
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  scale_shape_discrete(name="Library size") +
  labs(x="Genera precision (%)", y="Genera classification proportion (%)") +
  theme(aspect.ratio=1) +
  guides( colour='none')
vsearch_abundance_precision
# ggsave('./images/06-precision-vsearch-abundance.png', vsearch_abundance_precision)

precision_compare_plot <- cowplot::plot_grid(
  nrow = 1, scale = c(1, .9),
  align='h', axis='tb',
  precision_plot + labs(title="NanoCLUST"),
  vsearch_abundance_precision + labs(title="VSEARCH")
)
ggsave('./images/06-precision-even-abundance.png', precision_compare_plot)



vsearch_precision_summary %>%
  filter(method == 'consensus') %>%
  ggplot(aes(x=genus_precision, y=genus_classification_prop_by_depth, shape=factor(library_size), colour=factor(rep))) +
  geom_point() +
  # expand_limits(x=c(.75, .825), y=c(.8, 1)) +
  facet_grid(rows=vars(min_cluster_size), cols=vars(method),
             labeller = as_labeller(\(x)
                                      ifelse(is.na(as.numeric(x)), "most abundant", paste0("min OTU size (", scales::percent(as.numeric(x)), ")"))
             )
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  scale_shape_discrete(name="Library size") +
  labs(x="Genera precision (%)", y="Genera classification proportion (%)") +
  theme(aspect.ratio=1)