library(phyloseq)
library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(scales)
library(tidytext)
library(microViz)
library(ggpattern)

source('./helpers/dnabarcoder.R')
source('./helpers/config.R')
source('./helpers/vsearch.R')
source('./helpers/manual-taxonomy.R')

samplesheet <- read_samplesheet(config)

nanoclust <- load_nanoclust_phyloseq_3(
  samplesheet=samplesheet, experiment=config$experiment_path, reads_per_sample=config$sample_depth, repetition=config$repetition,
  sequence_type='nanoclust_abundant', min_cluster_size = '0.005') %>%
  tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)

nanoclust_tax <- load_dna_barcoder_classification_nanoclust_2(
  experiment=config$experiment_path,
  reads_per_sample=config$sample_depth,
  repetition=config$repetition,
  sequence_type='nanoclust_abundant',
  min_cluster_size = '0.005'
)

# nanoclust_abundant <- load_nanoclust_phyloseq_2(
#   samplesheet=samplesheet, experiment=config$experiment_path, reads_per_sample=config$sample_depth, repetition=config$repetition,
#   sequence_type='nanoclust_abundant') %>%
#   tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)

# vsearch_x <- load_nanoclust_phyloseq_2(
#   samplesheet=samplesheet, experiment=config$experiment_path, reads_per_sample=config$sample_depth, repetition=config$repetition,
#   sequence_type='vsearch') %>%
#   tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)

vsearch <- load_vsearch_phyloseq(samplesheet, config$experiment_path, config$sample_depth, config$repetition) %>%
  tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F) %>%
  filter_taxa_by_thresh(0.0015)
taxa_names(vsearch) <- 1:ntaxa(vsearch)

ggplot_splitting <- function(phylo, df, taxa_labels=FALSE, ncol=NULL, maxY=2200) {
  colours <- c("0" = "#00A600", "1" = "#E6E600", "2" = "#ECB176", "3" = "#F27971", "4" = "#7F7F7F")
  levs <- seq(0,4)

  df %>%
    ggplot( aes(x=OTU, y=count, fill=factor(match, levels = levs)) ) +
    geom_col() +
    # geom_col_pattern(aes(pattern_shape = factor(match)), pattern = 'pch', pattern_density = 0.5) +
    geom_text(aes(label=ifelse(manual_match, "*", '')), hjust=-2.7, vjust=-.2) +
    geom_label(aes(label = count, fill=factor(match, levels = levs)), vjust = -0.2, size = 3, show.legend = FALSE) +
    facet_wrap(~barcode,
               labeller = as_labeller(\(x)  paste0(
                 samplesheet[x, 'UpdatedName'], ' (BC ', gsub('barcode','', x), ')'
               )),
               scales='free_x',
               ncol = ncol) +
    tidytext::scale_x_reordered(
      labels = function(x) {
        if (taxa_labels) {
          paste0(tax_table(phylo)[reorder_func(x), 'species'], " (OTU ", reorder_func(x), ")")
        } else {
          reorder_func(x)
        }
      # labels = \(x) paste0("(", reorder_func(x), ")")
      },
    ) +
    scale_fill_manual(
      values = colours,
      labels = c("0"= "species","1"= "genus","2"= "family","3"= "unmatched", "4" ="unclassified kingdom"),
      name = "rank of match"
    ) +
    scale_colour_manual(values = colours, guide='none') +
    # scale_fill_gradientn(colours = terrain.colors(5)) +
    labs(x="Cluster with classification and ID", y="Abundance") +
    coord_cartesian(clip="off") +
    expand_limits(y=c(0, maxY)) +
    theme(
      # aspect.ratio=4/3,
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

organise_df <- function(phylo, subset) {
  species_level_correct <- function(species_unite, species) {
    (!is.na(species_unite) & species == species_unite) |
      # accepted_synonyms %>% filter(alt_species == species & actual_species == species_unite) %>% nrow() > 0
      species %in% accepted_synonyms$alt_species & species_unite %in% accepted_synonyms$actual_species
  }

  genus_level_correct <- function(genus_unite, species_unite, genus) {
    genus == genus_unite | genus %in% accepted_synonyms$alt_genus & species_unite %in% accepted_synonyms$actual_species
  }

  family_level_correct <- function(family_unite, species_unite, family) {
    family == family_unite | family %in% accepted_synonyms$alt_family & species_unite %in% accepted_synonyms$actual_species
  }

  cbind(
    otu_table(phylo),
    tax_table(phylo)[,c('family', 'genus', 'species')]
  ) %>% data.frame() %>%
    tibble::rownames_to_column('OTU') %>%
    pivot_longer(cols = sample_names(phylo), names_to='barcode', values_to='count' ) %>%
    mutate(count = as.integer(count)) %>%
    filter(count > 0) %>%
    left_join(
      samplesheet %>% tibble::rownames_to_column('barcode') %>% select(contains("_unite") | contains('barcode')),
      join_by(barcode)
    )%>%
    mutate(
      match = ifelse(
        species_level_correct(species_unite, species), 0,
        ifelse(genus_level_correct(genus_unite, species_unite, genus), 1,
          ifelse(family_level_correct(family_unite, species_unite, family), 2,
            ifelse(family == 'unclassified kingdom', 4,3))
        )
      ),
      manual_match =
        species %in% accepted_synonyms$alt_species & species_unite %in% accepted_synonyms$actual_species |
          genus %in% accepted_synonyms$alt_genus & species_unite %in% accepted_synonyms$actual_species |
          family %in% accepted_synonyms$alt_family & species_unite %in% accepted_synonyms$actual_species
    ) %>%
    filter(barcode %in% subset) %>%
    mutate(
      OTU = tidytext::reorder_within(OTU, -count, barcode)
    )
}

plot_splitting <- function(phylo, subset, taxa_labels=FALSE, ncol=NULL, maxY=2200) {
  df <- organise_df(phylo, subset)
  ggplot_splitting(phylo, df, taxa_labels, ncol, maxY)
}

# nanoclust_splitting <- plot_splitting(nanoclust, 1:5)
plot_splitting(nanoclust, sample_names(nanoclust)[1:30])
plot_splitting(nanoclust, sample_names(nanoclust)[31:70])

puccinias <- plot_splitting(nanoclust, paste0('barcode', c(25, 27, 28, 36)), taxa_labels = TRUE, ncol = 4, maxY=2100) +
  labs(title = "UMAP + HDBSCAN (NanoCLUST)", x='') +
  # geom_text(aes(label=ifelse(grepl('^49|^57', OTU), "*", '')), hjust=-2.7, vjust=-.2) +
  theme(aspect.ratio = 1)
puccinias
ggsave('images/06-cluster-splitting-nanoclust-puccinia.png', puccinias)

puccinias_vsearch <- plot_splitting(vsearch, paste0('barcode', c(25, 27, 28, 36)), taxa_labels = TRUE, ncol = 4, maxY=2100) +
  labs(title = "VSEARCH") +
  # geom_text(aes(label=ifelse(grepl('^5', OTU), "*", '')), hjust=-3, vjust=-.2) +
  theme(aspect.ratio = 1)
puccinias_vsearch
ggsave('images/06-cluster-splitting-vsearch-puccinia.png', puccinias_vsearch)


cryptococcus <- plot_splitting(nanoclust, paste0('barcode', 58:62), taxa_labels = TRUE, ncol = 4, maxY=2100) +
  labs(title = "UMAP + HDBSCAN (NanoCLUST)") +
  theme(aspect.ratio = 1) +
  labs(x="")
cryptococcus_vsearch <- plot_splitting(vsearch, paste0('barcode', 58:62), taxa_labels = TRUE, ncol = 4, maxY=2100) +
  labs(title = "VSEARCH") +
  theme(aspect.ratio = 1)

ggsave('images/06-cluster-splitting-nanoclust-cryptococcus.png', cryptococcus)
ggsave('images/06-cluster-splitting-vsearch-cryptococcus.png', cryptococcus_vsearch)


cand_nanoclust <- nanoclust
cand_tax <- as.matrix(tax_table(cand_nanoclust))
incertae <- cand_tax[,'genus'] == 'Debaryomycetaceae gen Incertae sedis'
cand_tax[incertae, c('family', 'genus', 'species')] <- 'Debaryomycetaceae family'
tax_table(cand_nanoclust) <- tax_table(cand_tax)
candidas <- plot_splitting(cand_nanoclust, paste0('barcode', c(47:57)), taxa_labels = F, ncol = 4, maxY=2100) +
  labs(title = "UMAP + HDBSCAN") +
  theme(aspect.ratio = 1) +
  labs(x="")
candidas

cand_vsearch <- vsearch
cand_tax <- as.matrix(tax_table(cand_vsearch))
incertae <- cand_tax[,'genus'] == 'Debaryomycetaceae gen Incertae sedis'
cand_tax[incertae, c('family', 'genus', 'species')] <- 'Debaryomycetaceae family'
tax_table(cand_vsearch) <- tax_table(cand_tax)
plot_splitting(cand_vsearch, paste0('barcode', c(47:57)), taxa_labels = F, ncol = 4, maxY=2100) +
  labs(title = "VSEARCH")

candidas_plot_nanoclust <-
  plot_splitting(cand_nanoclust, paste0('barcode', c(47, 49, 53, 54)), taxa_labels = T, ncol = 5, maxY=2100) +
    theme(aspect.ratio = 1) +
    labs(title="UMAP + HDBSCAN (NanoCLUST)", x="")
ggsave('images/06-cluster-splitting-nanoclust-candida.png', candidas_plot_nanoclust)


candidas_plot_vsearch <-
    plot_splitting(cand_vsearch, paste0('barcode',  c(47, 49, 53, 54)), taxa_labels = T, ncol = 5, maxY=2100) +
      theme(aspect.ratio = 1) +
      labs(title="VSEARCH")
ggsave('images/06-cluster-splitting-vsearch-candida.png', candidas_plot_vsearch)


botrytis_plot_nanoclust <- plot_splitting(nanoclust, paste0("barcode", c(30, 31, 32)), taxa_labels = T) +
  labs(title="UMAP + HDBSCAN (NanoCLUST)") +
  theme(aspect.ratio = 1) +
  labs(x="")
botrytis_plot_vsearch <- plot_splitting(vsearch, paste0("barcode", c(30, 31, 32)), taxa_labels = T) +
  labs(title="VSEARCH") +
  theme(aspect.ratio = 1)
ggsave('images/06-cluster-splitting-nanoclust-botrytis.png', botrytis_plot_nanoclust)
ggsave('images/06-cluster-splitting-vsearch-botrytis.png', botrytis_plot_vsearch)

plot_splitting(vsearch, sample_names(vsearch)[1:30])
plot_splitting(vsearch, sample_names(vsearch)[31:70])

ggsave('images/06-cluster-splitting-nanoclust-with-tax-1-30.png', plot_splitting(nanoclust, sample_names(nanoclust)[1:30]))
ggsave('images/06-cluster-splitting-nanoclust-with-tax-31.png', plot_splitting(nanoclust, sample_names(nanoclust)[31:70]))

ggsave('images/06-cluster-splitting-vsearch-with-tax-1-30.png', plot_splitting(vsearch, sample_names(vsearch)[1:30]))
ggsave('images/06-cluster-splitting-vsearch-with-tax-31.png', plot_splitting(vsearch, sample_names(vsearch)[31:70]))

write_tax_to_table <- function(phylo) {
  tax_to_show <- tax_table(phylo) %>% as.matrix()
  tax_to_show[tax_to_show[,'genus'] == 'Debaryomycetaceae gen Incertae sedis', c('family', 'genus', 'species') ] <- 'Debaryomycetaceae family'
  tax_to_show[tax_to_show[,'genus'] == 'Pucciniaceae gen Incertae sedis', c('family', 'genus', 'species') ] <- 'Pucciniaceae family'
  tax_to_show <- merge(tax_to_show, nanoclust_tax[, c('score', 'cutoff', 'confidence')], by = 0) %>% arrange(Row.names) %>%
    rename(OTU.ID = 'Row.names')
  tax_to_show
}

write_tax_to_table(nanoclust) %>%
  write.csv('tables/all-otus-tax.csv', row.names = F)

write_tax_to_table(vsearch) %>%
  write.csv('tables/all-otus-tax-vsearch.csv', row.names = F)

left_join(tax_to_show, nanoclust_tax, join_by())



nanoclust_uneven <- load_nanoclust_phyloseq_3(
  samplesheet=samplesheet, experiment="../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-uneven-reps-09-16", reads_per_sample=1000, repetition=config$repetition,
  sequence_type='nanoclust_consensus', min_cluster_size = '0.005') %>%
  tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)

plot_splitting(nanoclust_uneven, sample_names(nanoclust_uneven)[1:30])
plot_splitting(nanoclust_uneven, sample_names(nanoclust_uneven)[31:70])

ntaxa(nanoclust_uneven)
plot_splitting(nanoclust_uneven, paste0('barcode', c(77, 35, 39, 71, 57, 43)))

# seq_type_comparison <- "../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14-NC"
# nanoclust_ma <- load_nanoclust_phyloseq_2(samplesheet, seq_type_comparison, sequence_type="nanoclust_abundant",reads_per_sample=1000, repetition=1) %>%
#   tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)
# nanoclust_cons <- load_nanoclust_phyloseq_2(samplesheet, seq_type_comparison, sequence_type="nanoclust_consensus",reads_per_sample=1000, repetition=1) %>%
#   tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)
#
# plot_splitting(nanoclust_cons, sample_names(nanoclust_cons)[31:70])
# plot_splitting(nanoclust_ma, sample_names(nanoclust_ma)[31:70])
#
