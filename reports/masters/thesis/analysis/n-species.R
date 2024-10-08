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
config$sample_depth <- 2000

total_reads <- config$sample_depth * 58
min_cluster_sizes <- c(0, 0.001, .0025, .005)
loss_df <- data.frame(method=character(), min_cluster_size=double(), loss=double())
df <- NULL
# for (sample_depths in c(2000, 2500)) {
  for (min_cluster_size in min_cluster_sizes) {
    phylo <- load_nanoclust_phyloseq_3(
      samplesheet,
      experiment = config$experiment_path,
      sequence_type = 'nanoclust_abundant',
      reads_per_sample = config$sample_depth,
      repetition = config$repetition,
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
    c_size <- if (min_cluster_size == 0) { 2/total_reads } else { min_cluster_size }
    otu_tax_and_actual$min_cluster_size <- c_size
    # otu_tax_and_actual$sample_depth <- sample_depth
    df <- rbind(df, otu_tax_and_actual)

    loss_df <- loss_df %>%
      add_row(method = "nanoclust", min_cluster_size = c_size,  loss=1 - (sum(sample_sums(phylo))/total_reads))
  }
# }


df_vsearch <- NULL
low_end <- 1/total_reads
for (min_cluster_size in c(low_end, 0.00017, 0.0005, 0.001, .0025, .005)) {
  phylo <- load_vsearch_phyloseq(samplesheet,
                        experiment = config$experiment_path,
                        reads_per_sample = config$sample_depth,
                        repetition = config$repetition,
                        consensus = F
  )

  filtered_phylo <- filter_taxa_by_thresh(phylo, min_cluster_size)
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
  df_vsearch <- rbind(df_vsearch, otu_tax_and_actual)

  loss_df <- loss_df %>%
    add_row(method = "vsearch", min_cluster_size = min_cluster_size,  loss=1 - (sum(sample_sums(filtered_phylo))/total_reads))
}

expected_samples <- samplesheet[,c('family_unite', 'genus_unite', 'species')] %>%
  rename(family = 'family_unite', genus = 'genus_unite') %>%
  rows_append(tibble(family='Pucciniaceae', genus='Puccinia', species='Puccinia psidii')) %>%
  rows_append(tibble(family='Debaryomycetaceae', genus='Kurtzmaniella', species=NA))

# number of OTUs
# n_otus <- df %>%
#   group_by(min_cluster_size) %>%
#   distinct(otu) %>%
#   count()
#
# # number of species
# n_species <- df %>%
#   group_by(min_cluster_size, otu) %>%
#   summarise_all(first) %>%
#   group_by(
#     min_cluster_size,
#     species_identified = species != 'unidentified',
#   ) %>%
#   summarise(n=if(species_identified[1]) { n_distinct(species) } else {n()}) %>%
#   group_by(min_cluster_size) %>%
#   summarise(n = sum(n))
#
#
# n_taxa_from_samples <- df %>%
#   group_by(min_cluster_size, otu) %>%
#   summarise_all(first) %>%
#   mutate(
#     id_level = ifelse(
#       species != 'unidentified',
#       'species',
#       ifelse(
#         genus != 'unidentified',
#         'genus',
#         ifelse(
#           family != 'unidentified',
#           'family',
#           'unidentified (family and above)'
#         )
#       )
#     )
#   ) %>%
#   group_by(
#     min_cluster_size,
#     id_level
#   ) %>%
#   summarise(
#        n=if(id_level[1] == 'species') { n_distinct(species) } else { n() },
#        n_matching = case_when(
#          id_level[1] == 'species' ~ length(intersect(species, expected_samples$species)),
#          id_level[1] == 'genus' ~ sum(genus %in% expected_samples$genus),
#          id_level[1] == 'family' ~ sum(family %in% expected_samples$family),
#          TRUE ~ n()
#        ),
#        false_positives = sum(n) - sum(n_matching)
#   )
#
# fp_sums <- n_taxa_from_samples %>%
#   group_by(min_cluster_size) %>%
#   summarise(false_positives = sum(false_positives)) %>%
#   rename(n_matching = false_positives) %>%
#   mutate(id_level = 'false positive')
#
# stats <- rbind(n_taxa_from_samples, fp_sums)

# df %>% filter(min_cluster_size==0 & species=='unidentified'& genus != 'unidentified' & !(genus %in% expected_samples$genus)) %>% group_by(otu) %>% summarise_all(first) %>% View()
#
# df %>% group_by(min_cluster_size) %>% distinct(species) %>% count()

lvls <- c("false positive", "unidentified (family and above)", "family", "genus", "species")
cols <- c("#F27971",  "#7F7F7F",  "#ECB176", "#E6E600", "#00A600")
names(cols) <- lvls
# stats %>%
#   ggplot(aes(x=min_cluster_size*2000*58)) +
#   # geom_col(aes(y=n, fill=factor(id_level, levels=c("other", "family", "genus", "species")))) +
#   geom_col(aes(y=n_matching, fill=factor(id_level, levels=lvls)), width=10) +
#   # geom_point(aes(y=n_matching, colour=factor(id_level, levels=c("other", "family", "genus", "species")))) +
#   geom_hline(aes(yintercept=55, linetype='55'), show.legend = TRUE) +
#   labs(x="minimum cluster size (reads)", y="number of species", title = "NanoCLUST (2000 reads per sample)") +
#   scale_fill_manual(
#     values = cols,
#     name="OTU identification level",
#     labels = c("species not from sample set", "unidentified (family and above)", "family in sample set", "genus in sample set", "species in sample set")
#   ) +
#   scale_y_continuous(limits = c(0, 205)) +
#   scale_linetype_manual(name='Expected number of species', values="dashed", drop=FALSE) +
#   guides( fill = guide_legend(override.aes = c(linetype = 0)) )
  # geom_line(aes(x=min_cluster_size*2000*58, y=n), data=n_otus) #+
  # geom_line(aes(x=min_cluster_size*2000*58, y=n), data=n_species, linetype='dotted', colour='red')


# nanoclust <- load_nanoclust_phyloseq_3(samplesheet,
#                           experiment = config$experiment_path,
#                           sequence_type = 'nanoclust_abundant',
#                           reads_per_sample = config$sample_depth,
#                           repetition = config$repetition,
#                           min_cluster_size = 0.005
# ) %>% tax_fix(min_length = 0, unknowns = 'unidentified')
#
# samples <- samplesheet[,c('family_unite', 'genus_unite', 'species')] %>%
#   rename(family = 'family_unite', genus = 'genus_unite')
#
# notus <- dim(tax_table(nanoclust))[1]
# unique_species <- tax_table(nanoclust)[,c('genus','species')] %>%
#   as.data.frame()%>%
#   distinct()
#
# sum(unique_species$species %in% sampled_species) / length(unique_species$species)
# unique_species[!(unique_species$species %in% sampled_species),]
# # 96/184 # in sample sheet at 0
# # 47/78 # in sample sheet at 0.001
# # 36/58 # in sample sheet at 0.005
#
# id_species_tax_table <- tax_table(nanoclust)[tax_table(nanoclust)[,'species'] != 'unidentified',]
# mismatch <- id_species_tax_table[!(id_species_tax_table[, 'species'] %in% sampled_species), ]
#
# merge(mismatch, otu_table(nanoclust), by=0) %>%
#   pivot_longer(barcode25:barcode89, names_to = 'barcode', values_to = 'count') %>%
#   filter(count > 0) %>%
#   left_join(
#     (samplesheet %>% tibble::rownames_to_column("barcode") %>% select(barcode, species)),
#     join_by(barcode)
#   ) %>%
#   View()

n_otus <- df_vsearch %>%
  group_by(min_cluster_size) %>%
  distinct(otu) %>%
  count()

plot <- function(df, expected_samples, title, limits=c(0,210)) {
  # number of species
  n_species <- df %>%
    group_by(min_cluster_size, otu) %>%
    summarise_all(first) %>%
    group_by(
      min_cluster_size,
      species_identified = species != 'unidentified',
    ) %>%
    summarise(n=if(species_identified[1]) { n_distinct(species) } else {n()}) %>%
    group_by(min_cluster_size) %>%
    summarise(n = sum(n))


  n_taxa_from_samples <- df %>%
    group_by(min_cluster_size, otu) %>%
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
    group_by(min_cluster_size) %>%
    summarise(false_positives = sum(false_positives)) %>%
    rename(n_matching = false_positives) %>%
    mutate(id_level = 'false positive')

  stats <- rbind(n_taxa_from_samples, fp_sums)

  stats %>%
    # filter(min_cluster_size != 0)  %>%
    ggplot(aes(x=min_cluster_size*total_reads)) +
    # geom_col(aes(y=n, fill=factor(id_level, levels=c("other", "family", "genus", "species")))) +
    geom_col(aes(y=n_matching, fill=factor(id_level, levels=lvls)), width=10) +
    # geom_point(aes(y=n_matching, colour=factor(id_level, levels=c("other", "family", "genus", "species")))) +
    geom_hline(aes(yintercept=55, linetype='55'), show.legend = TRUE) +
    labs(x="minimum cluster size (reads)", y="number of species", title = title) +
    scale_fill_manual(
      values = cols,
      name="OTU identification level",
      labels = c("species not from sample set", "unidentified (family and above)", "family in sample set", "genus in sample set", "species in sample set")
    ) +
    # geom_text(aes(label=min_cluster_size*2000*58),y=0) +
    scale_linetype_manual(name='Expected number of species', values="dashed", drop=FALSE) +
    # scale_y_sqrt() +
    scale_y_continuous(limits = limits) +
    scale_x_continuous(breaks = round(stats$min_cluster_size*total_reads)) +
    # scale_y_continuous(transform = 'log10', breaks=c(55, 100, 1000)) +
    guides( fill = guide_legend(override.aes = c(linetype = 0)) )
}

cowplot::plot_grid(
  ncol = 1, axis='lr', align='hv',
  plot(df, expected_samples, title="NanoCLUST") + guides(fill='none', linetype='none'),
  plot(df_vsearch, expected_samples, title="VSEARCH" ) ,
  loss_df %>%
    ggplot(aes(x=min_cluster_size*total_reads, y=loss, colour=method)) +
      geom_line() +
      geom_point() +
      scale_y_continuous(labels = \(y) label_percent()(as.numeric(y)))
)
