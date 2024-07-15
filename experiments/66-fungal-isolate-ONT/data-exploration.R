library(phyloseq)
library(tidyr)

samplesheet <- read.csv('Rep1_Run2_ITS_Fungal_database_samplesheet.csv', header = TRUE, row.names = 1)
rownames(samplesheet) <- paste0("barcode", rownames(samplesheet))

isolate.98.FULL_ITS <- readRDS('outputs/isolate-98pct-full300-06-28/phyloseq/FULL_ITS/clustered-by-region/blast/clustered-by-region.phyloseq.rds')
sample_data(isolate.98.FULL_ITS) <- samplesheet
isolate.98.FULL_ITS
isolate.98.FULL_ITS.total <- sum(sample_sums(isolate.98.FULL_ITS))
isolate.98.FULL_ITS.total*0.0004
filter_taxa(isolate.98.FULL_ITS, \(x) sum(x) > isolate.98.FULL_ITS.total*0.0004, prune=TRUE)
  1 - (filter_taxa(isolate.98.FULL_ITS, \(x) sum(x) > isolate.98.FULL_ITS.total*0.0004, prune=TRUE) |>
  sample_sums() |> sum() / isolate.98.FULL_ITS.total)
# isolate.98.FULL_ITS.totals <- sum(taxa_sums(isolate.98.FULL_ITS))
# filter_taxa(isolate.98.FULL_ITS, \(x) sum(x) > isolate.98.FULL_ITS.totals*.0003, prune=TRUE) |>
#   view_tax() |> View()
# tax_glom(isolate.98.FULL_ITS, taxrank = "Species")

#isolate.98.filt.FULL_ITS <- readRDS('outputs/isolate-98pct-full300-06-28/phyloseq/')
#sample_data(isolate.98.filt.FULL_ITS) <- samplesheet


# isolate.97.FULL_ITS <- readRDS('outputs/isolate-97pct-06-24/phyloseq/FULL_ITS/clustered-by-region/blast/clustered-by-region.phyloseq.rds')
# sample_data(isolate.97.FULL_ITS) <- samplesheet
# isolate.97.FULL_ITS
# tax_glom(isolate.97.FULL_ITS, taxrank = "Species")
# plot_bar(tax_glom(isolate.97.FULL_ITS, taxrank = "Species"), x='sample_Sample', fill='Species')


# x <- read.csv('outputs/isolate-97pct-full300-07-01/vsearch-otu-map/FULL_ITS/all_samples/all_samples.otus.tsv', sep = '\t', row.names = 1) |>
#   otu_table(taxa_are_rows = TRUE)


# isolate.97.filt.FULL_ITS <- readRDS('outputs/isolate-97pct-full300-06-25/phyloseq/FULL_ITS/clustered-by-region/blast/clustered-by-region.phyloseq.rds')
# isolate.97.filt.FULL_ITS <- readRDS('outputs/isolate-97pct-full300-07-01/phyloseq/FULL_ITS/clustered-by-region/blast/clustered-by-region.phyloseq.rds')
isolate.97.filt.FULL_ITS <- readRDS('outputs/isolate-even-run2-07-12/phyloseq/FULL_ITS/clustered-by-region/blast/clustered-by-region.phyloseq.rds')
sample_data(isolate.97.filt.FULL_ITS) <- samplesheet

isolate.97.filt.FULL_ITS
plot_bar(isolate.97.filt.FULL_ITS)
plot_bar(tax_glom(isolate.97.filt.FULL_ITS, taxrank = "Species"), x='sample_Sample', fill='Species')

total_reads <- sum(sample_sums(isolate.97.filt.FULL_ITS))

view_tax <- function(phylo) {
  m <- merge(
    tax_table(phylo),
    data.frame(otu_table(phylo)),
    by=0
  )
  m <- m[, -1]
  colnames(m) <- c(colnames(tax_table(phylo)), sample_data(isolate.97.filt.FULL_ITS)$Sample)
  m
}

min_pct <- 0.0006

filter_taxa(isolate.97.filt.FULL_ITS, \(x) sum(x) > total_reads*min_pct, prune=TRUE) |>
  view_tax() |> View()
filter_taxa(isolate.97.filt.FULL_ITS, \(x) sum(x) > total_reads*min_pct, prune=TRUE)
total_reads*min_pct
1 - (sum(sample_sums(filter_taxa(isolate.97.filt.FULL_ITS, \(x) sum(x) > total_reads*min_pct, prune=TRUE))) / total_reads)

f <- \(t) sum(filter_taxa(isolate.97.filt.FULL_ITS, \(x) sum(x) > total_reads*t))
g <- Vectorize(f)
f.2 <- Vectorize(\(t) 100 * sum(sample_sums(filter_taxa(isolate.97.filt.FULL_ITS, \(x) sum(x) > total_reads*t, prune = TRUE))) / total_reads)
curve(g, from = 0.00005, to=0.025, n=200)
curve(f.2, from = 0.00005, to=0.025, n=200, add=TRUE,col='red')
legend(x="topright", legend = c("# of OTUs", "% reads"), lty=c(1,1), col=c("black", "red"))

plot_OTU_threshold(isolate.97.filt.FULL_ITS, 0.00005, 0.03)
plot_OTU_threshold(isolate.98.FULL_ITS, 0.00005, 0.03)

plot_OTU_threshold <- function(phylo, from, to) {
  total_reads <- sum(sample_sums(phylo))
  f <- \(t) sum(filter_taxa(phylo, \(x) sum(x) > total_reads*t))
  g <- Vectorize(f)
  f.2 <- Vectorize(\(t) 100 * sum(sample_sums(filter_taxa(phylo, \(x) sum(x) > total_reads*t, prune = TRUE))) / total_reads)
  curve(g, from = from, to=to, n=200)
  curve(f.2, from = from, to=to, n=200, add=TRUE,col='red')
  legend(x="topright", legend = c("# of OTUs", "% reads"), lty=c(1,1), col=c("black", "red"))
}

# nanoclust <- read.csv('outputs/nanoclust-kmer6/non-chimeras-clustering.tsv', sep = '\t', header = 1)
nanoclust_lsu <- read.csv('outputs/nanoclust-kmer6/isolate-06-25-kmer6/non-chimeras-clustering.tsv', sep = '\t', header = 1)
# nanoclust_fullits <- read.csv('outputs/nanoclust-kmer6/isolate-97pct-07-01/fullits.all_samples.nonchimeras.rerep.tsv', sep = '\t', header = 1)
nanoclust_fullits <- read.csv('outputs/nanoclust-kmer6/isolate-even-run2-07-12/clusters.tsv', sep = '\t', header = 1)
nanoclust_min72 <- read.csv('outputs/nanoclust-kmer6/isolate-even-run2-07-12/clusters-min72.tsv', sep = '\t', header = 1)
sums_min72 <- nanoclust_min72 %>% dplyr::count(bin_id)
sums_min72[1, 'n'] / sum(sums_min72$n)

nanoclust_min150 <- read.csv('outputs/nanoclust-kmer6/isolate-even-run2-07-12/clusters-min150.tsv', sep = '\t', header = 1)
sums_min150 <- nanoclust_min150 %>% dplyr::count(bin_id)
sums_min150[1, 'n'] / sum(sums_min150$n)

nanoclust_min200 <- read.csv('outputs/nanoclust-kmer6/isolate-even-run2-07-12/clusters-min200.tsv', sep = '\t', header = 1)
sums_min200 <- nanoclust_min200 %>% dplyr::count(bin_id)
sums_min200[1, 'n'] / sum(sums_min200$n)

nanoclust_min250 <- read.csv('outputs/nanoclust-kmer6/isolate-even-run2-07-12/clusters-min250.tsv', sep = '\t', header = 1)
sums_min250 <- nanoclust_min250 %>% dplyr::count(bin_id)
sums_min250[1, 'n'] / sum(sums_min250$n)

nanoclust_min350 <- read.csv('outputs/nanoclust-kmer6/isolate-even-run2-07-12/clusters-min350.tsv', sep = '\t', header = 1)
sums_min350 <- nanoclust_min350 %>% dplyr::count(bin_id)
sums_min350[1, 'n'] / sum(sums_min350$n)


View(nanoclust_fullits[nanoclust['bin_id'] != -1, ])
factor(nanoclust_fullits$bin_id)

# sum(nanoclust['bin_id'] == 666)
nanoclust_sums_fullits <- nanoclust_fullits |> dplyr::count(bin_id)
dim(nanoclust_sums_fullits)
hist(x = nanoclust_sums_fullits$n)

nanoclust_total_reads <- sum(nanoclust_sums_fullits$n)
nanoclust_total_reads*.001
nanoclust_max_clusters <- nanoclust_sums_fullits[nanoclust_sums_fullits$n > nanoclust_total_reads*.001, ]
dim(nanoclust_max_clusters)[1]
sum(nanoclust_max_clusters$n) / nanoclust_total_reads

nanoclust_w_barcode <- tidyr::separate(nanoclust_fullits, read, sep=';barcodelabel=', into=c("read id", "barcode"))

tabulated <- table(nanoclust_w_barcode[, c('barcode', 'bin_id')]) |> data.frame()
class(tabulated)

selected_clusters <- nanoclust_sums_fullits[nanoclust_sums_fullits$n > 464, ]$bin_id
tabulated[tabulated$Freq != 0 & (tabulated$bin_id %in% selected_clusters) & tabulated$bin_id != -1, ] |> View()


# derep.OTUS <- read.table('outputs/isolate-even-derep-07-15/vsearch-cluster/FULL_ITS/all_samples/all_samples.otu.tsv.gz', sep="\t", header=TRUE, comment.char = "", row.names=1)
# p <- phyloseq(otu_table(derep.OTUS, taxa_are_rows = TRUE))
p <- readRDS('outputs/isolate-even-derep-07-15/phyloseq/FULL_ITS/clustered-by-region/blast/clustered-by-region.phyloseq.rds')
sample_data(p) <- samplesheet

sum(sample_sums(p))*.0006
ntaxa(filter_taxa(p, \(x) sum(x) > 72, prune=TRUE))
1 - sum(sample_sums(filter_taxa(p, \(x) sum(x) > 72, prune=TRUE))) / sum(sample_sums(p))
plot_bar(p)
plot_bar(filter_taxa(p, \(x) sum(x) > 72, prune=TRUE), fill='Species', x = 'sample_Sample')
p1 <- filter_taxa(p, \(x) sum(x) > 72, prune=TRUE)
p2 <- prune_samples(sample_names(p1)[1:10], p1)
p3 <- prune_taxa(taxa_sums(p2) > 3, p2 )
plot_bar( p3, fill='Species', x = 'sample_Sample')
sample_data(p)
plot_OTU_threshold(p, 0.00005, 0.03)

