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

# samplesheet <- read.csv('../../../../experiments/66-fungal-isolate-ONT/Rep1_Run2_ITS_Fungal_database_samplesheet.csv', header = TRUE, row.names = 1)
samplesheet <- read.csv(config$samplesheet_path, header = TRUE, row.names = 1)
rownames(samplesheet) <- paste0("barcode", rownames(samplesheet))

experiment <- "../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14"

nanoclust <- load_nanoclust_phyloseq(samplesheet, config$experiment_path, config$sample_depth, config$repetition)

tax_with_counts <- cbind.data.frame(
    tax_table(nanoclust),
    otu_table(nanoclust)
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

genus_classification_proportion <- sum(classification_counts$genus_classified) / sum(classification_counts$total)
genus_classification_proportion
sum(classification_counts$species_classified) / sum(classification_counts$total)

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

genus_precision <- sum(precision_counts$genus_correct) / sum(precision_counts$total)
genus_precision

cowplot::plot_grid(ncol = 1,
classification_counts %>%
  mutate(proportion = genus_classified / total) %>%
  ggplot(aes(x=barcode, y=proportion)) +
    geom_col() +
    geom_hline(yintercept =  genus_classification_proportion) +
    scale_x_discrete(
      labels = \(x) sample_data(nanoclust)[x]$UpdatedName,
      limits = rownames(sample_data(nanoclust))[order(sample_data(nanoclust)$UpdatedName)]
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  precision_counts %>%
    mutate(precision = genus_correct / total) %>%
    ggplot(aes(x=barcode, y= precision)) +
    geom_col() +
    geom_hline(yintercept = genus_precision) +
    scale_x_discrete(
      labels = \(x) sample_data(nanoclust)[x]$UpdatedName,
      limits = rownames(sample_data(nanoclust))[order(sample_data(nanoclust)$UpdatedName)]
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
)
