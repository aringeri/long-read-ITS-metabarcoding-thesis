library(tidyr)
library(dplyr)
library(glue)

source('./helpers/config.R')
source('./helpers/qc-stats.R')

samplesheet <- read_samplesheet(config) %>%
  tibble::rownames_to_column("barcode") %>%
  mutate(barcode = gsub('barcode', '', barcode))

excluded_even <- c(
  'Aspergillus_niger' = '44',
  'Cryptococcus_albidus' = '58',
  'Galactomyces_geotrichum' = '64',
  'Meyerozyma_guillermondii' = '69',
  'Yarrowia_lipolytica' = '82',
  'Fusarium_proliferatum' = '63',
  'Puccinia_triticina' = '26'
)

excluded_uneven <- c(
  'Aspergillus_niger' = '44',
  'Cryptococcus_albidus' = '58',
  'Galactomyces_geotrichum' = '64',
  'Meyerozyma_guillermondii' = '69',
  'Yarrowia_lipolytica' = '82',
  'Fusarium_proliferatum' = '63',
  'Puccinia_triticina' = '26',
  'Austropuccinia_psidii' = '36',
  'Cryptococcus_gattii_VG_III' = '60',
  'Cryptococcus_neoformans_VNI' = '61'
)

fmted <- read_counts_by_sample_before_and_after_qc(samplesheet = samplesheet) %>%
  mutate(
    qcd_reads = ifelse(is.na(qcd_reads), NA, paste0(qcd_reads, " (", qcd_prop, ")")),
  ) %>%
  select(BC, sample, num_seqs, qcd_reads) %>%
  mutate(
    excluded_even = BC %in% excluded_even,
    excluded_uneven = BC %in% excluded_uneven
  ) %>%
  rename(
    Barcode = 'BC',
    Sample = 'sample',
    'Number of raw reads'=num_seqs,
    'Number of reads after QC'=qcd_reads,
  )

write.csv(fmted, '../../../../experiments/66-fungal-isolate-ONT/sample-stats-fmt.csv', row.names = F)
read.csv('../../../../experiments/66-fungal-isolate-ONT/sample-stats-fmt.csv') %>%
  rename_with(\(col) gsub('\\.', ' ', col)) %>% View()



