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

samplesheet <- read_samplesheet(config)

nanoclust <- load_nanoclust_phyloseq_3(
  samplesheet=samplesheet, experiment=config$experiment_path, reads_per_sample=config$sample_depth, repetition=config$repetition,
  sequence_type='nanoclust_abundant', min_cluster_size = '0.005') %>%
  tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)

# nanoclust_abundant <- load_nanoclust_phyloseq_2(
#   samplesheet=samplesheet, experiment=config$experiment_path, reads_per_sample=config$sample_depth, repetition=config$repetition,
#   sequence_type='nanoclust_abundant') %>%
#   tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)

# vsearch_x <- load_nanoclust_phyloseq_2(
#   samplesheet=samplesheet, experiment=config$experiment_path, reads_per_sample=config$sample_depth, repetition=config$repetition,
#   sequence_type='vsearch') %>%
#   tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)

vsearch <- load_vsearch_phyloseq(samplesheet, config$experiment_path, config$sample_depth, config$repetition) %>%
  filter_taxa_by_thresh(0.0015)

ggplot_splitting <- function(phylo, df, taxa_labels=FALSE, ncol=NULL, maxY=2200) {
  colours <- c("0" = "#00A600", "1" = "#E6E600", "2" = "#ECB176", "3" = "#F27971")

  df %>%
    ggplot( aes(x=OTU, y=count, fill=factor(match, levels = seq(0,3))) ) +
    geom_col() +
    # geom_col_pattern(aes(pattern_shape = factor(match)), pattern = 'pch', pattern_density = 0.5) +
    geom_label(aes(label = count, fill=factor(match, levels = seq(0,3))), vjust = -0.2, size = 3, show.legend = FALSE) +
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
      }
    ) +
    scale_fill_manual(
      values = colours,
      labels = c("0"= "species","1"= "genus","2"= "family","3"= "unmatched"),
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
      # match = ifelse(
      #   !is.na(species_unite) & species == species_unite, "species",
      #   ifelse(genus == genus_unite, "genus",
      #     ifelse(family == family_unite, "family",  "mismatch")
      #   )
      # )
      match = ifelse(
        !is.na(species_unite) & species == species_unite, 0,
        ifelse(genus == genus_unite, 1,
          ifelse(family == family_unite, 2,  3)
        )
      ),
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
  labs(title = "UMAP + HDBSCAN")
puccinias
ggsave('images/06-cluster-splitting-nanoclust-puccinia.png', puccinias)

taxa_names(vsearch) <- 1:ntaxa(vsearch)
puccinias_vsearch <- plot_splitting(vsearch, paste0('barcode', c(25, 27, 28, 36)), taxa_labels = TRUE, ncol = 4, maxY=2100) +
  labs(title = "VSEARCH")
puccinias_vsearch
ggsave('images/06-cluster-splitting-vsearch-puccinia.png', puccinias_vsearch)


cryptococcus <- plot_splitting(nanoclust, paste0('barcode', 58:62), taxa_labels = TRUE, ncol = 4, maxY=2100) +
  labs(title = "UMAP + HDBSCAN")
cryptococcus_vsearch <- plot_splitting(vsearch, paste0('barcode', 58:62), taxa_labels = TRUE, ncol = 4, maxY=2100) +
  labs(title = "VSEARCH")
# ggsave('images/06-cluster-splitting-nanoclust-cryptococcus.png', cryptococcus)
# ggsave('images/06-cluster-splitting-vsearch-cryptococcus.png', cryptococcus_vsearch)

plot_splitting(vsearch, sample_names(vsearch)[1:30])
plot_splitting(vsearch, sample_names(vsearch)[31:70])

# # ggsave('images/06-cluster-splitting-nanoclust-with-tax-6.png', plot_splitting(nanoclust, 1:6))
# # ggsave('images/06-cluster-splitting-nanoclust-with-tax-ex.png', plot_splitting(nanoclust, 2))
ggsave('images/06-cluster-splitting-nanoclust-with-tax-1-30.png', plot_splitting(nanoclust, sample_names(nanoclust)[1:30]))
ggsave('images/06-cluster-splitting-nanoclust-with-tax-31.png', plot_splitting(nanoclust, sample_names(nanoclust)[31:70]))

# seq_type_comparison <- "../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14-NC"
# nanoclust_ma <- load_nanoclust_phyloseq_2(samplesheet, seq_type_comparison, sequence_type="nanoclust_abundant",reads_per_sample=1000, repetition=1) %>%
#   tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)
# nanoclust_cons <- load_nanoclust_phyloseq_2(samplesheet, seq_type_comparison, sequence_type="nanoclust_consensus",reads_per_sample=1000, repetition=1) %>%
#   tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)
#
# plot_splitting(nanoclust_cons, sample_names(nanoclust_cons)[31:70])
# plot_splitting(nanoclust_ma, sample_names(nanoclust_ma)[31:70])
#
