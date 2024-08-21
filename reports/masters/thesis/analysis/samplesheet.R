library(tidyr)
library(dplyr)



samplesheet <- read.csv('../../../../experiments/66-fungal-isolate-ONT/ITS_Fungal_database_samplesheet_updated_2024-08.csv')

fmt_names <- samplesheet %>%
  tibble() %>%
  filter(UpdatedName != 'Negative_control') %>%
  separate_wider_delim(UpdatedName, names=c('genus', 'species_2', 'strain'), delim='_', too_many = 'merge', too_few = 'align_start', cols_remove = F) %>%
  mutate(species = paste0(genus, " ", species_2))

uniteDB <- read.csv('/Users/alex/Documents/uni/fungi-metabarcoding/databases/unite2024/unite2024ITS.classification', sep = '\t')
distinct_db <- distinct(uniteDB, species, .keep_all = T) %>%
  rename_with(~ paste0(.x, "_unite"))

species_level <- fmt_names %>% rename_with(~ paste0(.x, "_sample")) %>%
  left_join(
    distinct_db,
    join_by(species_sample == species_unite),
    keep = T
  )

genus_level <- species_level[is.na(species_level[,'id_unite']), ] %>%
  select(!contains("_unite")) %>%
  left_join(
    distinct(distinct_db, genus_unite, .keep_all = T),
    join_by(genus_sample == genus_unite),
    keep = T
  ) %>%
  mutate(
    species_unite = NA,
    SH_unite = NA
  )

species_level[is.na(species_level[,'id_unite']), ] <- genus_level

write.csv(
  species_level %>%
    select(!c(genus_sample, species_2_sample, strain_sample, id_unite)) %>%
    rename_with(~ gsub("_sample", "", .x))
  ,
  '../../../../experiments/66-fungal-isolate-ONT/samplesheet-unite-tax.csv',
  row.names = F
)
