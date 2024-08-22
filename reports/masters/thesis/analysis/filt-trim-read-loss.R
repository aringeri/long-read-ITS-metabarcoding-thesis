library(glue)
library(ggplot2)
library(scales)
library(magrittr)
library(dplyr)

source('./helpers/config.R')

# output_dir <- '../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-02'
# output_dir <- '../../../../experiments/66-fungal-isolate-ONT/outputs/isolate-even-reps-08-08'
output_dir <- config$experiment_path
# output_dir <- '/Users/alex/repos/long-read-ITS-metabarcoding/output/isolate-sub5000'

stages <- list(
  A_raw=list(dir=glue('{output_dir}/QC/01-raw_reads_all_samples/nanoplot/all-samples/NanoStats.txt'), name="Raw Reads"),
  B_adapter_trim=list(dir=glue('{output_dir}/QC/XX-adapter-trimming/nanoplot/all-samples/NanoStats.txt'), name="Adapter Trimming\n(Dorado)"),
  C_primer_trim=list(dir=glue('{output_dir}/QC/02-post_primer_trimming/nanoplot/all-samples/NanoStats.txt'), name="Primer Trimming\n(cutadapt)"),
  D_its_extraction=list(dir=glue('{output_dir}/QC/03-post_its_extraction/FULL_ITS/nanoplot/all-samples/NanoStats.txt'), name="Full ITS Extraction\n(ITSxpress)"),
  E_quality_filtering=list(dir=glue('{output_dir}/QC/04-post_quality_filtering/FULL_ITS/nanoplot/all-samples/NanoStats.txt'), name="Quality/Length Filtering\n(chopper)"),
  F_chimera_filtering=list(dir=glue('{output_dir}/QC/05-post_chimera_filtering/FULL_ITS/nanoplot/all-samples/NanoStats.txt'), name="Chimera Filtering\n(VSEARCH)")
)

df <- data.frame()
for (i in seq_along(stages)) {
  stage <- names(stages)[[i]]
  stage_read_data <- read.csv(stages[[stage]][['dir']], sep='\t', row.names = 1, stringsAsFactors = FALSE)
  stage_read_data <- data.frame(t(stage_read_data), stringsAsFactors = FALSE) %>%
    mutate_at(vars(number_of_reads), as.numeric)
  stage_read_data$stage <- stage
  stage_read_data$stage_name <- stages[[stage]][['name']]
  if (i > 1) {
    stage_read_data$pct_loss <- - ((df[df$stage == names(stages)[[i-1]], 'number_of_reads'] - stage_read_data$number_of_reads) / df[df$stage == names(stages)[[i-1]], 'number_of_reads'])
  } else {
    stage_read_data$pct_loss <- 0
  }
  df <- rbind(df, stage_read_data)
}

read_loss <- df %>%
  ggplot(
    aes(x=stage,
        y=number_of_reads,
        colour = pct_loss
    )
  ) +
  labs(x="pipeline stage", y="number of reads") +
  geom_col() +
  geom_text(aes(label = label_percent(0.01)(pct_loss)), vjust = -0.5) +
  scale_x_discrete(labels = \(x) lapply(stages[match(x, names(stages))], \(e) e[['name']])) +
  scale_color_continuous(name = 'relative read loss', labels = scales::percent) +
  scale_y_continuous(breaks = pretty_breaks(n=10))

  # scale_colour_viridis_c(option='plasma')
  # scale_colour_distiller(palette = "Spectral", direction = 1)

ggsave('images/06-read-loss-by-stage.png', read_loss)