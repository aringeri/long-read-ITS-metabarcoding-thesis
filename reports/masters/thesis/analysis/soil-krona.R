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
library(decontam)
library("MiscMetabar")


source('helpers/dnabarcoder.R')
source('helpers/vsearch.R')
source('helpers/config.R')

soil_full_dataset <- '../../../../experiments/camille-soil-samples/outputs/c-soil-singleton-10-17/'
soil_subset_dataset <- '../../../../experiments/camille-soil-samples/outputs/c-soil-NC-10-16/'
experiment <- soil_full_dataset

otus <- readRDS(glue('{experiment}/phyloseq/FULL_ITS/ALL_READS/1/all_samples/all_samples.phyloseq.rds'))
classFile <- list.files(
  glue('{experiment}/dnabarcoder/vsearch/FULL_ITS/ALL_READS/1/classify/'),
  pattern='*.unite2024ITS_BLAST.classification',
  full.names = T
)
soil_full_classifications <- read.csv(classFile, sep='\t', row.names = 1) %>%
  tibble::rownames_to_column('OTU') %>%
  separate(OTU, into=c('read', 'barcode', 'size'), sep=';') %>%
  mutate_at(vars(read, barcode, size), \(str) gsub('.*=', '', str)) %>%
  tibble::column_to_rownames('read')

sample_ids <- data.frame(
  row.names = colnames(otus),
  name = gsub('sample_', '', colnames(otus)),
  is_control = grepl('Con', colnames(otus))
) %>%
  mutate(type = substr(name, 1, 2))
tax_table <- read_dna_barcoder_classification_vsearch(classFile[1])
soil_full_phylo <- phyloseq(otu_table(otus), tax_table, sample_data(sample_ids))

contams <- isContaminant(soil_full_phylo, method='prevalence', neg='is_control')#, threshold=0.5)

soil_full_phylo_noncontam <- prune_taxa(!contams$contaminant, soil_full_phylo)

p <- soil_full_phylo_noncontam %>%
  filter_otu_by_sample(0.0015) %>%
  prune_samples(!sample_ids$is_control, .)

# allocasuarina
AL <- p %>% prune_samples(sample_data(p)$type == 'AL', .)
AP <- p %>% prune_samples(sample_data(p)$type == 'AP', .)

# eucalyptus
EC <- p %>% prune_samples(sample_data(p)$type == 'EC', .)
EV <- p %>% prune_samples(sample_data(p)$type == 'EV', .)

krona(AL, file="krona/AL.html")
krona(AP, file="krona/AP.html")
krona(EC, file="krona/EC.html")
krona(EV, file="krona/EV.html")
