library(tidyr)
library(dplyr)
library(glue)
library(ggplot2)

source('./helpers/config.R')
source('./helpers/qc-stats.R')

samplesheet <- read_samplesheet(config) %>%
  tibble::rownames_to_column("barcode") %>%
  mutate(barcode = gsub('barcode', '', barcode))

fmted <- read_counts_by_sample_before_and_after_qc(samplesheet=samplesheet)


cutoff <- 2500
sample_read_distribution <- cowplot::plot_grid(
  fmted %>%
    ggplot(aes(x=num_seqs+1)) +
    geom_histogram(bins=66) +
    # scale_x_continuous(n.breaks = 16, name = "Number of reads") +
    scale_x_log10(breaks = c(0, 10, 100, 1000, cutoff, 10000, 40000, 110000), limits = c(1, max(fmted$num_seqs)), name="Number of reads") +
    scale_y_continuous(limits = c(0,16),name = "Number of samples") +
    labs(title='A) Raw Reads') +
    theme(aspect.ratio = 1),
  fmted %>%
    ggplot(aes(x=as.numeric(qcd_reads)+1)) +
    geom_histogram(bins=66) +
    geom_vline(xintercept = cutoff, linetype='dashed') +
    scale_x_log10(breaks = c(0, 10, 100, 1000, cutoff, 10000, 40000, 110000), limits = c(1, max(fmted$num_seqs)), name="Number of reads") +
    scale_y_continuous(limits = c(0,16), name = "Number of samples") +
    labs(title='B) Reads after QC') +
    theme(aspect.ratio = 1)
)
sample_read_distribution

ggsave('images/06-read-count-distribution-by-sample.png', sample_read_distribution)



