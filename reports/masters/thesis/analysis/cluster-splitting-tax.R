library(phyloseq)
library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(scales)
library(tidytext)
library(microViz)

source('./helpers/dnabarcoder.R')
source('./helpers/config.R')

samplesheet <- read_samplesheet(config)

nanoclust <- load_nanoclust_phyloseq(samplesheet, config$experiment_path, config$sample_depth, config$repetition) %>%
  tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)

ggplot_splitting <- function(phylo, df, taxa_labels=FALSE) {
  df %>%
    ggplot( aes(x=OTU, y=count, fill=match) ) +
    geom_col() +
    geom_text(aes(label = count, colour=match), vjust = -0.5) +
    facet_wrap(~barcode,
               labeller = as_labeller(\(x)  paste0(
                 samplesheet[x, 'UpdatedName'], ' (', gsub('barcode','', x), ')'
               )),
               scales='free_x') +
    tidytext::scale_x_reordered(
      labels = function(x) {
        if (taxa_labels) {
          paste0(tax_table(phylo)[reorder_func(x), 'species'], " (OTU ", reorder_func(x), ")")
        } else {
          paste0("(", reorder_func(x), ")")
        }
      # labels = \(x) paste0("(", reorder_func(x), ")")
      }
    ) +
    scale_fill_manual(values = c("green4", "green3", "red4", "green2")) +
    scale_colour_manual(values = c("green4", "green3", "red4", "green2"), guide='none') +
    labs(x="Cluster", y="Abundance") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
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
      match = ifelse(
        !is.na(species_unite) & species == species_unite, "species",
        ifelse(genus == genus_unite, "genus",
          ifelse(family == family_unite, "family",  "mismatch")
        )
      )
    ) %>%
    filter(barcode %in% subset) %>%
    mutate(
      OTU = tidytext::reorder_within(OTU, -count, barcode)
    )
}

plot_splitting <- function(phylo, subset, taxa_labels=FALSE) {
  df <- organise_df(phylo, subset)
  ggplot_splitting(phylo, df, taxa_labels)
}

# nanoclust_splitting <- plot_splitting(nanoclust, 1:5)
plot_splitting(nanoclust, sample_names(nanoclust)[1:30])
plot_splitting(nanoclust, sample_names(nanoclust)[31:70])

# plot_splitting(nanoclust, paste0("barcode", c(57, 66, 68)), taxa_labels=TRUE)

# ggsave('images/06-cluster-splitting-nanoclust-with-tax-6.png', plot_splitting(nanoclust, 1:6))
# ggsave('images/06-cluster-splitting-nanoclust-with-tax-ex.png', plot_splitting(nanoclust, 2))

