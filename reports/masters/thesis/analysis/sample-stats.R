library(tidyr)


read.csv('../../../../experiments/66-fungal-isolate-ONT/sample-stats.tsv', sep='\t') %>%
  separate_wider_delim(file, '.', names_sep = "_") %>%
  mutate(
    BC = gsub('BC', '', file_4),
    sample = gsub('_', ' ', file_3)
  ) %>%
  select(BC, sample, num_seqs) %>%
  arrange(BC) %>%
  rename(Barcode = 'BC', Sample = 'sample',  'Number of reads'=num_seqs)