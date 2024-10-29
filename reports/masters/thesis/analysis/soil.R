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

n_identified <- sum(taxa_sums(prune_taxa(as.logical(tax_table[,'kingdom'] != 'unidentified'), soil_full_phylo)))
n_total <- 1848487
n_unidentified / n_total

contams <- isContaminant(soil_full_phylo, method='prevalence', neg='is_control')#, threshold=0.5)

soil_full_phylo_noncontam <- prune_taxa(!contams$contaminant, soil_full_phylo)

soil_full_phylo_noncontam %>%
  # filter_otu_by_sample(0.0015) %>%
  plot_richness()
  estimate_richness(measures = 'Observed')

soil_full_phylo_noncontam %>%
  filter_otu_by_sample(0.0015) %>%
  plot_heatmap()



?ordinate(soil_full_phylo_noncontam)

p <- soil_full_phylo_noncontam %>%
  filter_otu_by_sample(0.0015) %>%
  prune_samples(!sample_ids$is_control, .)


plot_bar(p %>%
  transform_sample_counts(\(x) x/sum(x)) %>%
  # prune_taxa(as.logical(tax_table(p)[, 'phylum'] != 'unidentified'), .) %>%
  tax_glom('class')
  ,
  # x='type',
 fill='class')# + facet_wrap(~ name)

o <- ordinate(p, distance='bray')
plot_ordination(p, o, color='type')

# make unique taxa if unidentified at any taxonomic level
df <- data.frame()
tax_df <- tax_table %>% as.data.frame()
for (i in 1:nrow(tax_table)) {
  new_row <- tax_df[i,]
  new_row[, tax_df[i,] == 'unidentified' ] <- rownames(new_row)
  df <- rbind(df,
    new_row
  )
}
new_tax_table <- tax_table(as.matrix(df))
new_soil_full_phylo <- phyloseq(new_tax_table, otu_table(soil_full_phylo))
tax_glom(new_soil_full_phylo, taxrank='species')
tax_table(new_soil_full_phylo)[, 'species'] %>% head()

stats <- data.frame()
min_otu_threshes <- seq(from = 0, to=0.1, by=.001)
for (thresh in min_otu_threshes) {
  filtered <- filter_otu_by_sample(new_soil_full_phylo, thresh)
  stats <- rbind(stats,
    data.frame(
      min_out = thresh,
      n_otus = filtered %>% ntaxa(),
      n_species = filtered %>% tax_table() %>% as.data.frame() %>% select(species) %>% distinct() %>% summarise(n=n())
    )
  )
}

stats %>%
  ggplot(aes(x=min_out, y=n_otus)) +
  geom_point() +
  geom_line(aes(y=n), colour='red')

filt <- soil_full_phylo #%>% filter_taxa_by_min_otu_size(700)
filt



ntaxa(soil_full_phylo)
estimate_richness(soil_full_phylo, measures = 'Observed')
soil_full_phylo


top_30_otus <- taxa_sums(soil_full_phylo) %>% data.frame(count = .) %>%
  tibble::rownames_to_column("OTU") %>%
  arrange(desc(count)) %>%
  slice_head(n=30)

pruned_top_30 <- prune_taxa(top_30_otus$OTU, soil_full_phylo) %>%
  tax_fix(unknowns = 'unidentified', anon_unique = F)

unite2024 <- read.csv('/Users/alex/Documents/uni/fungi-metabarcoding/databases/unite2024/unite2024ITS.classification', sep = '\t')


x <- tax_table(pruned_top_30)[, 'species'] %>%
  merge(top_30_otus, by.x = 0, by.y = 'OTU') %>%
  rename(OTU = Row.names) %>%
  arrange(desc(count)) %>%
  left_join(
    soil_full_classifications %>% tibble::rownames_to_column("OTU") %>%
      select(OTU, score, cutoff, ReferenceID),
    join_by(OTU)
  )

common <- x %>%
  left_join(unite2024, join_by(ReferenceID == id)) %>%
  select(ReferenceID:species.y) %>%
  filter(ReferenceID != '') %>%
  distinct(ReferenceID, .keep_all = TRUE) %>%
  tibble::column_to_rownames("ReferenceID") %>%
  rename(species=species.y)

common['KP889944', c('order', 'family', 'genus')] <- 'unidentified'
common['DQ309136', c('class', 'order', 'family', 'genus')] <- 'unidentified'
common['UDB0343814', c('family', 'genus')] <- 'unidentified'
common['UDB0757809', 'genus'] <- 'unidentified'

uniteRefs <- tax_table(as.matrix(common)) %>%
  # View()
  tax_fix(unknowns = 'unidentified') %>%
  as.data.frame() %>%
  tibble::rownames_to_column('ReferenceID')

top_30_full_soil_with_closest_match <- x %>%
  rename(classification = species) %>%
  left_join(
    uniteRefs %>%
      select(ReferenceID, species) %>%
      rename(closest_match = species),
    join_by(ReferenceID)
  ) %>%
  select(classification:closest_match) %>%
  arrange(desc(count))

write.csv(top_30_full_soil_with_closest_match, file='./tables/soil-top30-otus-vsearch.csv', row.names = F)

load_tax_nc <- function(experiment, min_cluster_size) {
  classFile <- list.files(
    glue('{experiment}/dnabarcoder/nanoclust_abundant/FULL_ITS/4000/1/{min_cluster_size}/classify/'),
    pattern='*.unite2024ITS_BLAST.classification',
    full.names = T
  )
  tax <- read.csv(classFile[1], sep='\t', row.names = 1) %>%
    tibble::rownames_to_column('OTU') %>%
    separate(OTU, into=c('read', 'cluster', 'barcode', NA, 'size'), sep=';') %>%
    mutate_at(vars(read, barcode, cluster, size), \(str) gsub('.*=', '', str)) %>%
    filter(cluster != -1) %>%
    arrange(as.numeric(cluster)) %>%
    tibble::column_to_rownames("cluster")
  tax
}

load_nc <- function(experiment, min_cluster_size) {
  otu <- read.csv(glue('{experiment}/hdbscan_clustering/FULL_ITS/4000/1/{min_cluster_size}/otu_table/otu_table.tsv'), sep = '\t', row.names = 1)
  otu <- otu[order(as.numeric(rownames(otu))), ]
  otu <- otu[rownames(otu) != -1, ]

  tax <- load_tax_nc(experiment, min_cluster_size)
  tax <- tax_table(as.matrix(tax[, !(names(tax) %in% c('read', 'barcode', 'size', 'ReferenceID', 'rank', 'score', 'cutoff', 'confidence'))]))

  phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax)
}

subset_nc <- load_nc(soil_subset_dataset, 2)
sample_data(subset_nc) <- sample_ids
tax_nc <- load_tax_nc(soil_subset_dataset, 2)

nc_decontam <- isContaminant(subset_nc, method="prevalence", neg='is_control')
nc_p_decontam <- prune_taxa(!nc_decontam$contaminant, subset_nc)

plot_richness(nc_p_decontam, measures='Observed')

nc_o <- nc_p_decontam %>% prune_samples(!sample_ids$is_control, .) #%>%
  # filter_otu_by_sample(0.0015)

o_nc <- ordinate(nc_o, method='PCoA')
plot_scree(o_nc)
plot_ordination(nc_o, o_nc, color = 'type', label='name')

plot_net(
  nc_p_decontam %>% prune_samples(!sample_ids$is_control, .),
  point_label = "name", color = 'type', laymeth="graphopt")

top_30_otus_nc <- taxa_sums(subset_nc) %>% data.frame(count = .) %>%
  tibble::rownames_to_column("OTU") %>%
  arrange(desc(count)) %>%
  slice_head(n=30)

pruned_top_30_nc <- prune_taxa(top_30_otus_nc$OTU, subset_nc) %>%
  tax_fix(unknowns = 'unidentified', anon_unique = F)



top_30_nc_w_cutoff <- top_30_otus_nc %>%
  merge(
    tax_table(pruned_top_30_nc)[, 'species'] %>% as.data.frame() %>% rename(classification = species),
    by.x = 'OTU', by.y = 0) %>%
  merge(tax_nc %>% select(ReferenceID, score, cutoff),
    by.x = 'OTU', by.y = 0)


closest_match_nc <- left_join(top_30_nc_w_cutoff, unite2024, join_by(ReferenceID == id)) %>%
  select(ReferenceID:species) %>%
  filter(ReferenceID != '') %>%
  distinct(ReferenceID, .keep_all = TRUE) %>%
  tibble::column_to_rownames("ReferenceID") %>%
  as.matrix() %>% tax_table() %>% tax_fix(unknowns = 'unidentified') %>%
  as.data.frame() %>% tibble::rownames_to_column("ReferenceID") %>%
  select(ReferenceID, species)

top_30_nc_w_cutoff %>%
  left_join(closest_match_nc, join_by("ReferenceID")) %>%
  arrange(desc(count)) %>%
  rename(closest_match = species) %>%
  select(classification, count, score, cutoff, ReferenceID, closest_match) %>%
  write.csv(file='./tables/soil-top30-otus-nanoclust.csv', row.names = F)



# plot proportion identified
prop_id_vsearch <- soil_full_phylo %>%
  # filter_otu_by_sample(0.0015) %>%
  # filter_taxa_by_min_otu_size(50) %>%
  # transform_sample_counts(\(x) x/sum(x)) %>%
  # prune_samples('sample_AL4', .) %>%
  tax_glom(taxrank="kingdom") %>%
  plot_bar(fill="kingdom != \"unidentified\"") +
  # scale_y_continuous(name='Relative abundance') +
  labs(title="VSEARCH (full dataset)") +
  scale_fill_discrete(name='kingdom', labels = c('unidentified', 'identified')) +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
prop_id_vsearch

prop_id_nc <- subset_nc %>%
  # transform_sample_counts(\(x) x/sum(x)) %>%
  # prune_samples('sample_AL4', .) %>%
  tax_glom(taxrank="phylum") %>%
  plot_bar(fill="phylum != \"unidentified\"") +
  scale_fill_discrete(name='phylum', labels = c('unidentified', 'identified')) +
  # scale_y_continuous(name='Relative abundance') +
  labs(title="NanoCLUST (100K reads)") +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

unident_plot <- cowplot::plot_grid(ncol = 2, align='v', axis='tb',
  prop_id_vsearch + guides(fill='none'),
  prop_id_nc
)
unident_plot

ggsave('images/06-soil-unidentified.png', prop_id_vsearch)

subset_nc %>%
  filter_otu_by_sample(0.005) %>%
  tax_fix(unknowns = 'unidentified', anon_unique = F)  %>%
  prune_taxa(rownames(tax_table(subset_nc)[tax_table(subset_nc)[, 'kingdom'] != 'unidentified', ]), x=.) %>%
  transform_sample_counts(\(x) x/sum(x)) %>%
  # plot_bar( facet_grid = 'kingdom', fill='phylum')
  plot_bar( facet_grid = 'phylum == "Ascomycota"', fill='genus') +
  guides(fill='none')

soil_full_phylo %>%
  filter_otu_by_sample(0.005) %>%
  tax_fix(unknowns = 'unidentified', anon_unique = F)  %>%
  prune_taxa(rownames(tax_table(soil_full_phylo)[tax_table(soil_full_phylo)[, 'kingdom'] != 'unidentified', ]), x=.) %>%
  transform_sample_counts(\(x) x/sum(x)) %>%
  # tax_glom(taxrank = 'class') %>%
  plot_bar( facet_grid = 'phylum == "Ascomycota"', fill='genus') +
  guides(fill='none')

# tax_count <- cbind(tax_table(fix), taxa_sums(fix))
# un_id <- tax_count[tax_count[,'species'] %in% 1:700, ]
# un_id %>% as.data.frame() %>% mutate(count = as.numeric(V8)) %>% arrange(count) %>% summarise(x= sum(count))
#
# load_nanoclust_phyloseq_3(NULL, experiment = soil_subset_dataset, sequence_type="nanoclust_abundant", 4000, 1, 2)

