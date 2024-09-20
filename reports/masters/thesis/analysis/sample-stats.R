library(tidyr)
library(dplyr)
library(glue)

source('./helpers/config.R')
source('./helpers/qc-stats.R')

samplesheet <- read_samplesheet(config) %>%
  tibble::rownames_to_column("barcode") %>%
  mutate(barcode = gsub('barcode', '', barcode))

fmted <- read_counts_by_sample_before_and_after_qc(samplesheet = samplesheet) %>%
  mutate(
    qcd_reads = ifelse(is.na(qcd_reads), NA, paste0(qcd_reads, " (", qcd_prop, ")"))
  ) %>%
  select(BC, sample, num_seqs, qcd_reads) %>%
  rename(
    Barcode = 'BC',
    Sample = 'sample',
    'Number of raw reads'=num_seqs,
    'Number of reads after QC'=qcd_reads,
  )

write.csv(fmted, '../../../../experiments/66-fungal-isolate-ONT/sample-stats-fmt.csv', row.names = F)
read.csv('../../../../experiments/66-fungal-isolate-ONT/sample-stats-fmt.csv') %>%
  rename_with(\(col) gsub('\\.', ' ', col)) %>% View()



