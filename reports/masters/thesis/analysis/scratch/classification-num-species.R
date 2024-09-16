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

source('../helpers/dnabarcoder.R')
source('../helpers/config.R')

samplesheet <- read_samplesheet(config)

phylo <- load_nanoclust_phyloseq_3(samplesheet = samplesheet, config$experiment_path, min_cluster_size='0.001')
ntaxa(phylo)

tax <- as.data.frame(tax_table(phylo))
tax[tax[, 'species'] == 'unidentified', 'species'] <- paste0('unidentified species ', 1:(sum(tax[, 'species'] == 'unidentified')))
tax[tax[, 'genus'] == 'unidentified', 'genus'] <- paste0('unidentified genus ', 1:(sum(tax[, 'genus'] == 'unidentified')))

tax_table(phylo) <- tax_table(as.matrix(tax))

species <- tax_glom(phylo, 'species')
ntaxa(species)
plot_bar(species, fill = 'species')


genus <- tax_glom(phylo, 'genus')
ntaxa(genus)
plot_bar(genus, fill='genus')
