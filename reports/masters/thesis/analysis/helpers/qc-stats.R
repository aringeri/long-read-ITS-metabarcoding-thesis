library(tidyr)
library(dplyr)
library(glue)

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

read_counts_by_sample_before_and_after_qc <- function(
  qc_experiment_path='../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-all-sample-qc',
  raw_sample_stats_tsv='../../../../experiments/66-fungal-isolate-ONT/sample-stats.tsv',
  samplesheet
) {
  qc_stats <- gather_qc_stats(qc_experiment_path)

  stats <- read.csv(raw_sample_stats_tsv, sep='\t') %>%
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
      qcd_prop = ifelse(is.na(qcd_reads), NA, scales::percent(as.numeric(qcd_reads) / as.numeric(num_seqs), accuracy = 0.1))
    )
  return(stats)
}