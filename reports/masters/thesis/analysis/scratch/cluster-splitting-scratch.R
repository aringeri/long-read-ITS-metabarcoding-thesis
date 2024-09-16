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

load_abundance_vs_consensus <- function(config, subset, thresh='0', ncol=NULL, taxa_labels=TRUE) {
  nanoclust <- load_nanoclust_phyloseq_3(
    samplesheet=samplesheet, experiment=config$experiment_path, reads_per_sample=config$sample_depth, repetition=config$repetition,
    sequence_type='nanoclust_consensus', thresh=thresh) %>%
    # ) %>%
    tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)

  nanoclust_abundant <- load_nanoclust_phyloseq_3(
    samplesheet=samplesheet, experiment=config$experiment_path, reads_per_sample=config$sample_depth, repetition=config$repetition,
    sequence_type='nanoclust_abundant', thresh=thresh) %>%
    # )%>%
    tax_fix(min_length = 0, unknowns = 'unidentified', anon_unique = F)

  experiment=config$experiment_path
  reads_per_sample=config$sample_depth
  repetition=config$repetition
  sequence_type='nanoclust_consensus'

  classFile <- list.files(
    glue('{experiment}/dnabarcoder/{sequence_type}/FULL_ITS/{reads_per_sample}/{repetition}/{thresh}/classify/'),
    # glue('{experiment}/dnabarcoder/{sequence_type}/FULL_ITS/{reads_per_sample}/{repetition}/classify/'),
    pattern='*.unite2024ITS_BLAST.classification',
    full.names = T
  )
  tax <- load_dna_barcoder_classification_nanoclust(classFile[1])

  tax_abundant <- load_dna_barcoder_classification_nanoclust(list.files(
    glue('{experiment}/dnabarcoder/nanoclust_abundant/FULL_ITS/{reads_per_sample}/{repetition}/{thresh}/classify/'),
    # glue('{experiment}/dnabarcoder/nanoclust_abundant/FULL_ITS/{reads_per_sample}/{repetition}/classify/'),
    pattern='*.unite2024ITS_BLAST.classification',
    full.names = T
  )[1])

  list(
   consensus = plot_splitting(nanoclust, tax, subset,  taxa_labels=taxa_labels, ncol=ncol) + labs(title = "Consensus"),
   abundance= plot_splitting(nanoclust_abundant, tax_abundant, subset,  taxa_labels=taxa_labels, ncol=ncol) + labs(title = "Abundance")
  )
}


ggplot_splitting <- function(phylo, df, taxa_labels=FALSE, ncol=NULL, maxY=2200) {
  colours <- c("0" = "#00A600", "1" = "#E6E600", "2" = "#ECB176", "3" = "#F27971")

  df %>%
    ggplot( aes(x=OTU, y=count, fill=factor(match, levels = seq(0,3))) ) +
    geom_col() +
    # geom_col_pattern(aes(pattern_shape = factor(match)), pattern = 'pch', pattern_density = 0.5) +
    geom_label(aes(label = count, fill=factor(match, levels = seq(0,3))), vjust = -0.2, size = 3, show.legend = FALSE) +
    geom_text(
      aes(label = ifelse(is.na(score), NA, paste0("(", score, ", ", confidence, ")"))),
      vjust = 1.5, size = 2.5, show.legend = FALSE) +
    # geom_text(aes(label=count), x=0, y=0) +
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

organise_df <- function(phylo, tax, subset) {
  cbind(
    otu_table(phylo),
    tax_table(phylo)[,c('family', 'genus', 'species')],
    tax[, c('score', 'confidence')]
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
      score = ifelse(count > 0.25 * 2000, score, NA),
      confidence = ifelse(count > 0.25 * 2000, confidence, NA)
    ) %>%
    filter(barcode %in% subset) %>%
    mutate(
      OTU = tidytext::reorder_within(OTU, -count, barcode)
    )
}

plot_splitting <- function(phylo, tax, subset, taxa_labels=FALSE, ncol=NULL, maxY=2200) {
  df <- organise_df(phylo, tax, subset)
  ggplot_splitting(phylo, df, taxa_labels, ncol, maxY)
}

plot_splitting(nanoclust, tax, paste0('barcode', c(25, 27, 28, 36)),  taxa_labels=TRUE, ncol=4)
plot_splitting(nanoclust_abundant, tax_abundant, paste0('barcode', c(25, 27, 28, 36)),  taxa_labels=TRUE, ncol=4)

# puccinias
load_abundance_vs_consensus(config, paste0('barcode', c(25, 27, 28, 36)), ncol=4)
# cryptococcus
load_abundance_vs_consensus(config, paste0('barcode', 58:62), thresh='0.0065', ncol=4)

load_abundance_vs_consensus(config, paste0('barcode', c(31,32,46)), ncol=4)


load_abundance_vs_consensus(config, paste0('barcode', 25:55), taxa_labels = FALSE)
load_abundance_vs_consensus(config, paste0('barcode', 56:86), taxa_labels = FALSE)

