library(tidyr)
library(dplyr)

source('./helpers/config.R')

samplesheet <- read_samplesheet(config) %>%
  tibble::rownames_to_column("barcode") %>%
  mutate(barcode = gsub('barcode', '', barcode))

fmted <- read.csv('../../../../experiments/66-fungal-isolate-ONT/sample-stats.tsv', sep='\t') %>%
  separate_wider_delim(file, '.', names_sep = "_") %>%
  mutate(
    BC = gsub('BC', '', file_4),
    sample = gsub('_', ' ', file_3)
  ) %>%
  arrange(BC) %>%
  left_join(samplesheet, join_by(BC == barcode)) %>%
  mutate(
    sample = ifelse(OriginalName != UpdatedName, paste0(sample, " (", gsub('_', ' ', UpdatedName), ")"), sample)
  ) %>%
  select(BC, sample, num_seqs) %>%
  rename(Barcode = 'BC', Sample = 'sample',  'Number of raw reads'=num_seqs)

write.csv(fmted, '../../../../experiments/66-fungal-isolate-ONT/sample-stats-fmt.csv', row.names = F)
read.csv('../../../../experiments/66-fungal-isolate-ONT/sample-stats-fmt.csv') %>%
  rename_with(\(col) gsub('\\.', ' ', col)) %>% View()