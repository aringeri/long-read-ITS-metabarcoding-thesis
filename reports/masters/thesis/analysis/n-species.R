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

samplesheet <- read_samplesheet(config)

# scenario <- 'even'
scenario <- 'uneven'

if (scenario == 'even') {
  experiment_path <- "../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-fixed-2-10-10"
  n_samples <- 58
  total_reads <- config$sample_depth * n_samples
  sample_depths <- c(167, 1000, 2000, 2500)
  min_cluster_sizes <- c(2, 5, 10, 20, 50, 100)
  min_otu_sizes <- c(2, 5, 10, 20, 50, 100)
  repetition <- 1
} else { #uneven
  experiment_path <- "../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-uneven-reps-lo-multi-t-09-25"
  n_samples <- 55
  sample_depths <- c(100, 200, 500, 1000, 2000)
  min_cluster_sizes <- c(0, 5, 10, 20, 50, 100)
  min_otu_sizes <- c(2, 5, 10, 20, 50, 100)
  repetition <- 1
}

loss_df <- data.frame(method=character(), min_cluster_size=double(), loss=double(), sample_depth=integer())
df <- NULL
for (sample_depth in sample_depths) {
  for (min_cluster_size in min_cluster_sizes) {
    if (min_cluster_size <= sample_depth) {
      # print(paste0(min_cluster_size, " ", sample_depth))
      phylo <- load_nanoclust_phyloseq_3(
        samplesheet,
        experiment = experiment_path,
        sequence_type = 'nanoclust_abundant',
        reads_per_sample = sample_depth,
        repetition = repetition,
        min_cluster_size = min_cluster_size)

      otu_tax_and_actual <- merge(tax_table(phylo), otu_table(phylo), by=0) %>%
        pivot_longer(barcode25:barcode89, names_to = 'barcode', values_to = 'count') %>%
        filter(count > 0) %>%
        rename(otu = 'Row.names') %>%
        left_join(
          (samplesheet %>%
            tibble::rownames_to_column("barcode") %>%
            select(barcode, order_unite, family_unite, genus_unite, species) %>%
            rename(
              actual_order = order_unite,
              actual_family = family_unite,
              actual_genus = genus_unite,
              actual_species = species
            )
          ),
          join_by(barcode)
        )
      c_size <- if (min_cluster_size == 0) { 2 } else { min_cluster_size }
      otu_tax_and_actual$min_cluster_size <- c_size
      otu_tax_and_actual$sample_depth <- sample_depth
      df <- rbind(df, otu_tax_and_actual)

      total_reads <- if (scenario == 'even') sample_depth*n_samples else 50*sample_depth + sum(5, 10, 20, 50, 100)

      loss_df <- loss_df %>%
        add_row(
          sample_depth = sample_depth,
          method = "nanoclust", min_cluster_size = c_size,  loss = 1 - sum(sample_sums(phylo))/total_reads
        )
    }
  }
}


df_vsearch <- NULL
for (sample_depth in sample_depths) {
  for (min_cluster_size in min_otu_sizes) {
    if (min_cluster_size <= sample_depth) {
      phylo <- load_vsearch_phyloseq(samplesheet,
                            experiment = experiment_path,
                            reads_per_sample = sample_depth,
                            repetition = repetition,
                            consensus = F
      )

      filtered_phylo <- filter_taxa_by_min_otu_size(phylo, min_cluster_size)
      otu_tax_and_actual <- merge(tax_table(filtered_phylo), otu_table(filtered_phylo), by=0) %>%
        pivot_longer(barcode25:barcode89, names_to = 'barcode', values_to = 'count') %>%
        filter(count > 0) %>%
        rename(otu = 'Row.names') %>%
        left_join(
          (samplesheet %>%
            tibble::rownames_to_column("barcode") %>%
            select(barcode, order_unite, family_unite, genus_unite, species) %>%
            rename(
              actual_order = order_unite,
              actual_family = family_unite,
              actual_genus = genus_unite,
              actual_species = species
            )
          ),
          join_by(barcode)
        )
      otu_tax_and_actual$min_cluster_size <- min_cluster_size
      otu_tax_and_actual$sample_depth <- sample_depth
      df_vsearch <- rbind(df_vsearch, otu_tax_and_actual)

      total_reads <- if (scenario == 'even') sample_depth*n_samples else 50*sample_depth + sum(5, 10, 20, 50, 100)

      loss_df <- loss_df %>%
        add_row(
          sample_depth = sample_depth,
          method = "vsearch", min_cluster_size = min_cluster_size, loss = 1 - sum(sample_sums(filtered_phylo))/total_reads
        )
    }
  }
}

expected_samples <- samplesheet[,c('family_unite', 'genus_unite', 'species')] %>%
  rename(family = 'family_unite', genus = 'genus_unite') %>%
  rows_append(tibble(family='Pucciniaceae', genus='Puccinia', species='Puccinia psidii')) %>%
  rows_append(tibble(family='Debaryomycetaceae', genus='Kurtzmaniella', species=NA))

lvls <- c("false positive", "unidentified (family and above)", "family", "genus", "species")
cols <- c("#F27971",  "#7F7F7F",  "#ECB176", "#E6E600", "#00A600")
names(cols) <- lvls

n_otus <- df_vsearch %>%
  group_by(sample_depth,min_cluster_size) %>%
  distinct(otu) %>%
  count()
n_otus

n_otus_NC <- df %>%
  group_by(sample_depth,min_cluster_size) %>%
  distinct(otu) %>%
  count()

label_facet <- function(sample_depth) {
  depth <- as.numeric(sample_depth)
  if (scenario == 'even') {
    paste(depth*n_samples, " (", depth, " reads per taxa)")
  } else {
    paste0(50*depth + sum(5, 10, 20, 50, 100), " reads in library")
  }
}

plot <- function(df, expected_samples, title, limits=c(0,255)) {
  # number of species
  n_species <- df %>%
    group_by(sample_depth, min_cluster_size, otu) %>%
    summarise_all(first) %>%
    group_by(
      sample_depth,
      min_cluster_size,
      species_identified = species != 'unidentified',
    ) %>%
    summarise(n=if(species_identified[1]) { n_distinct(species) } else {n()}) %>%
    group_by(sample_depth, min_cluster_size) %>%
    summarise(n = sum(n))

  n_otus <- df %>%
    group_by(sample_depth,min_cluster_size) %>%
    distinct(otu) %>%
    count()

  n_taxa_from_samples <- df %>%
    group_by(sample_depth, min_cluster_size, otu) %>%
    summarise_all(first) %>%
    mutate(
      id_level = ifelse(
        species != 'unidentified',
        'species',
        ifelse(
          genus != 'unidentified',
          'genus',
          ifelse(
            family != 'unidentified',
            'family',
            'unidentified (family and above)'
          )
        )
      )
    ) %>%
    group_by(
      sample_depth,
      min_cluster_size,
      id_level
    ) %>%
    summarise(
      n=if(id_level[1] == 'species') { n_distinct(species) } else { n() },
      n_matching = case_when(
        id_level[1] == 'species' ~ length(intersect(species, expected_samples$species)),
        id_level[1] == 'genus' ~ sum(genus %in% expected_samples$genus),
        id_level[1] == 'family' ~ sum(family %in% expected_samples$family),
        TRUE ~ n()
      ),
      false_positives = sum(n) - sum(n_matching)
    )

  fp_sums <- n_taxa_from_samples %>%
    group_by(sample_depth, min_cluster_size) %>%
    summarise(false_positives = sum(false_positives)) %>%
    rename(n_matching = false_positives) %>%
    mutate(id_level = 'false positive')

  stats <- rbind(n_taxa_from_samples, fp_sums)

  stats %>%
    # filter(min_cluster_size != 0)  %>%
    ggplot(aes(x=min_cluster_size)) +
    # geom_col(aes(y=n, fill=factor(id_level, levels=c("other", "family", "genus", "species")))) +
    # geom_bar(
    #   aes(y = n_matching, group = min_cluster_size),
    #   colour = "black",
    #   stat = "summary", fun = sum,
    #   # fill = "transparent",
    #   width=1
    # ) +
    geom_col(aes(y=n_matching, fill=factor(id_level, levels=lvls)), width=.7) +
    # geom_point(aes(y=n_matching, colour=factor(id_level, levels=c("other", "family", "genus", "species")))) +
    geom_hline(aes(yintercept=55, linetype='55'), show.legend = TRUE) +
    geom_line(aes(y=n, x=min_cluster_size, colour='blue'), data=n_otus, show.legend = TRUE) +
    labs(x="minimum OTU size", y="number of species", title = title) +
    scale_fill_manual(
      values = cols,
      name="OTU identification level",
      labels = c("false positive", "unidentified (family and above)", "family in sample set", "genus in sample set", "species in sample set")
    ) +
    # geom_text(aes(label=min_cluster_size*2000*58),y=0) +
    scale_linetype_manual(name='Expected number of species', values="dashed", drop=FALSE) +
    facet_grid(
      cols = vars(sample_depth),
      labeller = as_labeller(label_facet)
    ) +
    # scale_y_continuous(limits = limits, breaks = c(0,50,100,150, 200, 250)) +
    scale_colour_manual(name="Total OTUs", values='blue', label='') +
    scale_y_continuous(breaks = c(0,50,100,150, 200, 250)) +
    coord_cartesian(ylim=limits) +
    scale_x_continuous(transform = 'sqrt', breaks = c(2, 5, 10, 20, 50, 100)) +
    # scale_y_continuous(transform = 'log10', breaks=c(55, 100, 1000)) +
    guides( fill = guide_legend(override.aes = c(linetype = 0)) )
}

lims <- if (scenario == 'even') c(0,255) else c(0, 180)

n_species_plot <- cowplot::plot_grid(
  ncol = 1, axis='lr', align='hv',
  plot(df, expected_samples, title="UMAP + HDBSCAN (NanoCLUST)", limits=lims) + guides(fill='none', linetype='none'),
  plot(df_vsearch, expected_samples, title="VSEARCH", limits = lims) + guides( colour='none') ,
  loss_df %>%
    ggplot(aes(x=min_cluster_size, y=loss, colour=method)) +
      geom_line() +
      geom_point() +
      scale_y_continuous(labels = \(y) label_percent()(as.numeric(y))) +
      scale_x_continuous(transform = 'sqrt', name="minimum OTU size", breaks = c(2, 5, 10, 20, 50, 100)) +
      scale_colour_discrete(name="Clustering method") +
      coord_cartesian(ylim=c(0, 0.075)) +
      facet_grid(
        cols = vars(sample_depth),
        labeller = as_labeller(label_facet)
      )
)
n_species_plot

if (scenario =='even') {
  ggsave("images/06-n-species-even.png", n_species_plot)
} else {
  ggsave("images/06-n-species-uneven.png", n_species_plot)
}
