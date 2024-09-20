library(tidyr)
library(dplyr)
library(glue)

source('./helpers/config.R')

samplesheet <- read_samplesheet(config) %>%
  tibble::rownames_to_column("barcode") %>%
  mutate(barcode = gsub('barcode', '', barcode))

gather_qc_stats <- function(experiment_path) {
  post_qc_stats <- data.frame(bc = NULL, qcd_reads = NULL)

  qc_dir <- glue('{experiment_path}/QC/05-post_chimera_filtering/FULL_ITS/nanoplot')

  sample_dirs <- list.dirs(
    qc_dir,
    recursive = FALSE, full.names = FALSE) %>%
    grep('BC.*', ., value = TRUE)

  for (sample_dir in sample_dirs) {
    sample_stats <- read.csv(glue('{qc_dir}/{sample_dir}/NanoStats.txt'), sep = '\t', row.names = 1)
    bc <- sub('BC', '', substring(sample_dir, regexpr('BC.*', sample_dir)))
    post_qc_stats[bc, 'qcd_reads'] <- sample_stats['number_of_reads', 'dataset']
  }
  post_qc_stats %>% tibble::rownames_to_column('BC')
}

qc_stats <- gather_qc_stats('../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-all-sample-qc')

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
  left_join(qc_stats, join_by(BC)) %>%
  mutate(
    qcd_reads = ifelse(is.na(qcd_reads), NA, paste0(qcd_reads, " (", scales::percent(as.numeric(qcd_reads) / as.numeric(num_seqs), accuracy = 0.1), ")"))
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



