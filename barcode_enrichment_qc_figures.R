library(dplyr)
library(monocle)
library(stringr)
library(grid)
source('helper_functions.R')
source('paths.R')
source('constants.R')

## Load in barcodes
mock_version1_ko_barcodes = read.delim(file.path(TEMP_DATA_DIR, "ko_barcodes/singles/version1/mock_ko_barcodes.txt"))
dox_treated_100nm_version1_ko_barcodes = read.delim(file.path(TEMP_DATA_DIR, "ko_barcodes/singles/version1/dox_treated_100nm_ko_barcodes.txt"))
dox_treated_100nm_rep2_version1_ko_barcodes = read.delim(file.path(TEMP_DATA_DIR, "ko_barcodes/singles/version1/dox_treated_100nm_rep2_ko_barcodes.txt"))
dox_treated_500nm_version1_ko_barcodes = read.delim(file.path(TEMP_DATA_DIR, "ko_barcodes/singles/version1/dox_treated_500nm_ko_barcodes.txt"))



# Load in data about mutants
guide_metadata = read.delim(file.path(DATA_DIR, "guide_gene_associations_17bp.txt"), header=F, col.names=c('gene', 'guide', 'category'))
target_metadata = read.delim(file.path(DATA_DIR, "gene_functional_categories.txt"), header=F, col.names=c('gene', 'functional_category'))
guide_metadata = merge(guide_metadata, target_metadata, by='gene')

## Get rid of unprocessed barcodes
mock_version1_ko_barcodes = process_barcodes(mock_version1_ko_barcodes, guide_metadata)
dox_treated_100nm_version1_ko_barcodes = process_barcodes(dox_treated_100nm_version1_ko_barcodes, guide_metadata)
dox_treated_100nm_rep2_version1_ko_barcodes = process_barcodes(dox_treated_100nm_rep2_version1_ko_barcodes, guide_metadata)
dox_treated_500nm_version1_ko_barcodes = process_barcodes(dox_treated_500nm_version1_ko_barcodes, guide_metadata)

# Make a plot of where the thresholds lie for each sample
mock_version1_ko_barcodes$sample = "mock V1"
dox_treated_100nm_version1_ko_barcodes$sample = "dox treated 100nm V1"
dox_treated_100nm_rep2_version1_ko_barcodes$sample = "dox treated 100nm rep2 V1"
dox_treated_500nm_version1_ko_barcodes$sample = "dox treated 500nm V1"

combined_data_singles = rbind(mock_version1_ko_barcodes, dox_treated_100nm_version1_ko_barcodes, dox_treated_100nm_rep2_version1_ko_barcodes, dox_treated_500nm_version1_ko_barcodes)

# Reorder the levels for plotting
combined_data_singles$sample = factor(combined_data_singles$sample, levels = c("mock V1", "dox treated 100nm V1", "dox treated 100nm rep2 V1", "dox treated 500nm V1"))

# Singles plot
ggplot(combined_data_singles, aes(x=proportion, y=log10(read_count))) +
    geom_point(size=0.15, alpha=0.5) +
    geom_vline(xintercept=KO_ASSIGNMENT_PROPORTION_THRESHOLD, color="red", size=1) +
    geom_hline(yintercept=log10(KO_ASSIGNMENT_READS_THRESHOLD), color="red", size=1, alpha=0.5) +
    facet_wrap(~sample, ncol=4) +
    xlab('barcode proportion') +
    ylab('log10(reads from barcode)') +
    theme_cfg() +
    theme(panel.margin = unit(1, "lines")) +
    ggsave(file.path(SUPPLEMENT_FIGURE_DIR, "barcode_enrichment_qc.singles.png"), height=5, width=8)
