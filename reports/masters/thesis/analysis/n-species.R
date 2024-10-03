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

min_cluster_sizes <- c(0, 0.001, .0025, .005, .0065)
df <- NULL
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
  otu_tax_and_actual$min_cluster_size <- min_cluster_size
  df <- rbind(df, otu_tax_and_actual)
}

expected_samples <- samplesheet[,c('family_unite', 'genus_unite', 'species')] %>%
  rename(family = 'family_unite', genus = 'genus_unite') %>%
  rows_append(tibble(family='Pucciniaceae', genus='Puccinia', species='Puccinia psidii')) %>%
  rows_append(tibble(family='Debaryomycetaceae', genus='Kurtzmaniella', species=NA))

sampled_species <- c(samplesheet$species, 'Puccinia psidii')

# number of OTUs
n_otus <- df %>%
  group_by(min_cluster_size) %>%
  distinct(otu) %>%
  count()

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
          'other'
        )
      )
    )
  ) %>%
  group_by(
    min_cluster_size,
    id_level
    # species_identified = species != 'unidentified',
    # genus_identified = genus != 'unidentified',
    # family_identifed = family != 'unidentified',
    # order_identifed = order != 'unidentified',
    # class_identifed = class != 'unidentified',
    # phylum_identifed = phylum != 'unidentified',
    # kingdom_identifed = kingdom != 'unidentified'
  ) %>%
  # distinct(species) %>%
  # summarise(n=n_distinct(family, genus, species))
  # summarise(n=n_distinct(family, genus, species))
  # reframe(n=if_else(family_identifed, n_distinct(family, genus, species), n()))
  # summarise(n=if(family_identifed[1]) { n_distinct(family, genus, species)} else {n()})
  summarise(
       n=if(id_level[1] == 'species') { n_distinct(species) } else { n() },
       n_matching = case_when(
         id_level[1] == 'species' ~ length(intersect(species, expected_samples$species)),
         id_level[1] == 'genus' ~ sum(genus %in% expected_samples$genus),
         id_level[1] == 'family' ~ sum(family %in% expected_samples$family),
         TRUE ~ n()
       )
       # n_matching=if(species_identified[1]) {
       #   length(intersect(species, expected_samples$species))
       # } else if (genus_identified[1]) {
       #   sum(genus %in% expected_samples$genus)
       # } else if (family_identifed[1]){
       #   sum(family %in% expected_samples$family)
       # } else { 0 }
  )

# df %>% filter(min_cluster_size==0 & species=='unidentified'& genus != 'unidentified' & !(genus %in% expected_samples$genus)) %>% group_by(otu) %>% summarise_all(first) %>% View()
#
# df %>% group_by(min_cluster_size) %>% distinct(species) %>% count()

n_taxa_from_samples %>%
  ggplot(aes(x=min_cluster_size*2000*58)) +
  # geom_col(aes(y=n, fill=factor(id_level, levels=c("other", "family", "genus", "species")))) +
  geom_col(aes(y=n_matching, fill=factor(id_level, levels=c("other", "family", "genus", "species")))) +
  # geom_point(aes(y=n_matching, colour=factor(id_level, levels=c("other", "family", "genus", "species")))) +
  geom_hline(yintercept = 55, linetype='dotted') +
  # geom_line(aes(x=min_cluster_size*2000*58, y=n), data=n_otus) #+
  geom_line(aes(x=min_cluster_size*2000*58, y=n), data=n_species, linetype='dotted', colour='red')


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
