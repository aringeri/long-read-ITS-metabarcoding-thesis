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

samplesheet <- read.csv('../../../../experiments/66-fungal-isolate-ONT/Rep1_Run2_ITS_Fungal_database_samplesheet.csv', header = TRUE, row.names = 1)
rownames(samplesheet) <- paste0("barcode", rownames(samplesheet))

experiment <- "../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-14"

nanoclust <- load_nanoclust_phyloseq(samplesheet, experiment, 2000, 2)


tax_with_counts <- cbind.data.frame(
  tax_table(nanoclust),
  otu_table(nanoclust)
) %>% pivot_longer(barcode25:barcode89, names_to='barcode', values_to = 'count')

genus_classification_counts <- tax_with_counts %>%
  group_by(barcode) %>%
  filter(count > 0) %>%
  summarise(
    classified=sum(count[genus!='unidentified']),
    total=sum(count)
  )

genus_classification_proportion <- sum(genus_classification_counts$classified) / sum(genus_classification_counts$total)
genus_classification_proportion

samplesheet_genus <- samplesheet %>%
  tibble::rownames_to_column('barcode') %>%
  separate(Sample, into='genus', sep='_', remove = FALSE)

genus_classification_precision_counts <- tax_with_counts %>%
  # group_by(barcode, genus) %>%
  filter(count > 0) %>%
  left_join(samplesheet_genus, join_by(barcode), suffix = c('.actual', '.expected')) %>%
  summarise(
    # .by = c(barcode, genus.actual, genus.expected),
    # total = sum(count)
    .by = barcode,
    correct = sum((genus.actual == genus.expected)*count),
    total = sum(count)
  )

genus_precision <- sum(genus_classification_precision_counts$correct) / sum(genus_classification_precision_counts$total)
genus_precision


cowplot::plot_grid(ncol = 1,
genus_classification_counts %>%
  mutate(proportion = classified / total) %>%
  ggplot(aes(x=barcode, y=proportion)) +
    geom_col() +
    geom_hline(yintercept =  genus_classification_proportion) +
    scale_x_discrete(
      labels = \(x) sample_data(nanoclust)[x]$Sample,
      limits = rownames(sample_data(nanoclust))[order(sample_data(nanoclust)$Sample)]
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  genus_classification_precision_counts %>%
    mutate(precision = correct / total) %>%
    ggplot(aes(x=barcode, y= precision)) +
    geom_col() +
    geom_hline(yintercept = genus_precision) +
    scale_x_discrete(
      labels = \(x) sample_data(nanoclust)[x]$Sample,
      limits = rownames(sample_data(nanoclust))[order(sample_data(nanoclust)$Sample)]
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
)

# tax_with_counts %>%
#   # group_by(barcode, genus) %>%
#   filter(count > 0) %>%
#   left_join(samplesheet_genus, join_by(barcode), suffix = c('.actual', '.expected')) %>%
#   summarise(
#     .by = c(barcode, species, genus.actual, genus.expected, Sample),
#     total = sum(count)
#   ) %>%
#   arrange(barcode) %>% View()

