library(glue)
library(ggplot2)
library(scales)
library(magrittr)
library(dplyr)

source('./helpers/config.R')

output_dir <- config$experiment_path
output_dir <- '../../../../experiments/camille-soil-samples/outputs/c-soil-10-03'
# output_dir <- '../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-all-sample-qc'

stages <- list(
  A_raw=list(dir=glue('{output_dir}/QC/01-raw_reads_all_samples/nanoplot/all-samples/NanoStats.txt'), name="Raw Reads"),
  B_adapter_trim=list(dir=glue('{output_dir}/QC/XX-adapter-trimming/nanoplot/all-samples/NanoStats.txt'), name="Adapter Trimming\n(Dorado)"),
  C_primer_trim=list(dir=glue('{output_dir}/QC/02-post_primer_trimming/nanoplot/all-samples/NanoStats.txt'), name="Primer Trimming\n(cutadapt)"),
  D_its_extraction=list(dir=glue('{output_dir}/QC/03-post_its_extraction/FULL_ITS/nanoplot/all-samples/NanoStats.txt'), name="Full ITS Extraction\n(ITSxpress)"),
  E_quality_filtering=list(dir=glue('{output_dir}/QC/04-post_quality_filtering/FULL_ITS/nanoplot/all-samples/NanoStats.txt'), name="Quality/Length Filtering\n(chopper)"),
  F_chimera_filtering=list(dir=glue('{output_dir}/QC/05-post_chimera_filtering/FULL_ITS/nanoplot/all-samples/NanoStats.txt'), name="Chimera Filtering\n(VSEARCH)")
)

stages <- list(
  A_raw=list(dir=glue('01-raw_reads_all_samples'), name="Raw Reads"),
  B_adapter_trim=list(dir=glue('XX-adapter-trimming'), name="Adapter Trimming\n(Dorado)"),
  C_primer_trim=list(dir=glue('02-post_primer_trimming'), name="Primer Trimming\n(cutadapt)"),
  D_its_extraction=list(dir=glue('03-post_its_extraction/FULL_ITS'), name="Full ITS Extraction\n(ITSxpress)"),
  E_quality_filtering=list(dir=glue('04-post_quality_filtering/FULL_ITS'), name="Quality/Length Filtering\n(chopper)"),
  F_chimera_filtering=list(dir=glue('05-post_chimera_filtering/FULL_ITS'), name="Chimera Filtering\n(VSEARCH)")
)

get_read_counts_for_all_stages <- function(output_dir, stages) {
  new_df <- data.frame()
  for (i in seq_along(stages)) {
    stage <- names(stages)[[i]]
    stage_subdir <- stages[[stage]][['dir']]
    stage_dir <- glue('{output_dir}/QC/{stage_subdir}/nanoplot')

    for (sample in list.dirs(stage_dir ,recursive = FALSE, full.names = FALSE)) {
      if (sample != 'all-samples') {
        # s_bc <- stringr::str_split_1(sample, ".BC")
        s_bc <- stringr::str_split_1(sample, "sample_")

        stage_read_data <- read.csv(glue('{stage_dir}/{sample}/NanoStats.txt'), sep='\t', row.names = 1, stringsAsFactors = FALSE)
        stage_read_data <- data.frame(t(stage_read_data), stringsAsFactors = FALSE) %>%
          mutate_at(vars(number_of_reads), as.numeric)
        stage_read_data$sample_name <- s_bc[1]
        stage_read_data$barcode <- s_bc[2]
        stage_read_data$stage <- stage
        stage_read_data$stage_name <- stages[[stage]][['name']]
        if (i > 1) {
          prev_stage_nreads <- new_df[new_df$stage == names(stages)[[i-1]] & new_df$barcode == s_bc[2], 'number_of_reads']
          stage_read_data$pct_loss <- - ((prev_stage_nreads - stage_read_data$number_of_reads) / prev_stage_nreads)
        } else {
          stage_read_data$pct_loss <- 0
        }

        new_df <- rbind(new_df, stage_read_data)
      }
    }
  }
  new_df
}

new_df <- get_read_counts_for_all_stages(output_dir, stages)

# df <- data.frame()
# for (i in seq_along(stages)) {
#   stage <- names(stages)[[i]]
#   stage_read_data <- read.csv(stages[[stage]][['dir']], sep='\t', row.names = 1, stringsAsFactors = FALSE)
#   stage_read_data <- data.frame(t(stage_read_data), stringsAsFactors = FALSE) %>%
#     mutate_at(vars(number_of_reads), as.numeric)
#   stage_read_data$stage <- stage
#   stage_read_data$stage_name <- stages[[stage]][['name']]
#   if (i > 1) {
#     stage_read_data$pct_loss <- - ((df[df$stage == names(stages)[[i-1]], 'number_of_reads'] - stage_read_data$number_of_reads) / df[df$stage == names(stages)[[i-1]], 'number_of_reads'])
#   } else {
#     stage_read_data$pct_loss <- 0
#   }
#   df <- rbind(df, stage_read_data)
# }

# read_loss <- df %>%
#   ggplot(
#     aes(x=stage,
#         y=number_of_reads,
#         colour = pct_loss
#     )
#   ) +
#   labs(x="Pipeline Stage", y="Number of Reads") +
#   geom_col() +
#   geom_text(aes(label = label_percent(0.01)(pct_loss)), vjust = -0.5) +
#   scale_x_discrete(labels = \(x) lapply(stages[match(x, names(stages))], \(e) e[['name']])) +
#   scale_color_continuous(name = 'Relative Read Loss', labels = scales::percent) +
#   scale_y_continuous(breaks = pretty_breaks(n=10))

  # scale_colour_viridis_c(option='plasma')
  # scale_colour_distiller(palette = "Spectral", direction = 1)

# omitted <- c(44, 58, 64, 69, 82, 63, 26)
omitted <- c()

sumrz <- new_df %>%
  mutate( included = !(barcode %in% omitted) ) %>%
  group_by(stage, included) %>%
  filter(included == T) %>%
  summarise(
    number_of_reads = mean(number_of_reads),
    mean_pct_loss = mean(pct_loss),
    barcode = NA
  )

low_samples <- new_df %>%
  filter(stage =='F_chimera_filtering' & number_of_reads < 45000) %>%
  select(barcode) %>%
  distinct()

plot <- new_df %>%
  # mutate( included = !(barcode %in% omitted) ) %>%
  # filter(barcode %in% c(25:30)) %>%
  ggplot(aes(x=stage, y=number_of_reads, colour=ifelse(barcode %in% low_samples$barcode, barcode, 'other'))) +
    geom_jitter(size=2, alpha=0.9, width = 0.1) +
    geom_crossbar(aes(ymin=number_of_reads, ymax=number_of_reads), data=sumrz, width=0.3, show.legend = F) +
    # geom_label(aes(label=round(number_of_reads)), data = sumrz, hjust = -.7) +
    geom_text(aes(x=stage, y=number_of_reads, label=scales::label_percent(0.01)(mean_pct_loss)), data = sumrz, hjust = -.7, show.legend = F) +
    # geom_boxplot() +
    # geom_point() +
    # facet_wrap(~barcode, ncol =1)
    scale_x_discrete(labels = \(x) lapply(stages[match(x, names(stages))], \(e) e[['name']])) +
    # scale_y_continuous(n.breaks = 20) +
    # scale_y_sqrt(breaks=c(0, 1000, 10000, 50000, 100000)) +
    scale_y_sqrt(breaks = c(0, 1000, 50000, 100000, 150000, 200000, 250000, 300000, 350000)) +
    expand_limits(y=0) +
    scale_colour_discrete(name='Sample with low read count', breaks=low_samples$barcode) +
    labs(x="Pipeline Stage", y="Number of Reads") +
    theme(aspect.ratio=.7)
plot

# ggsave('images/06-read-loss-by-stage.png', read_loss)
ggsave('images/06-read-loss-by-stage-by-sample.png', plot)
ggsave('images/06-read-loss-by-stage-by-sample-soil.png', plot)


rc <- get_read_counts_for_all_stages('../../../../experiments/camille-soil-samples/outputs/c-soil-10-03', stages)
rc %>%
  select(barcode, number_of_reads, stage) %>%
  filter(stage == 'A_raw' | stage == 'F_chimera_filtering') %>%
  tidyr::pivot_wider(names_from = stage, values_from = number_of_reads) %>%
  mutate(
    qcd_prop = scales::percent(as.numeric(F_chimera_filtering) / as.numeric(A_raw), accuracy = 0.1)
  ) %>%
  mutate(F_chimera_filtering = paste0(F_chimera_filtering, " (", qcd_prop, ")")) %>%
  select(barcode, A_raw, F_chimera_filtering) %>%
  rename(Sample=barcode, 'Number of raw reads'=A_raw, "Number of reads after QC"=F_chimera_filtering) %>%
  write.csv('../../../../experiments/camille-soil-samples/soil-sample-stats-fmt.csv', row.names = F)