library(monocle)
library(argparse)
library(dplyr)

# Define script
parser = ArgumentParser(description="Script to shuffle a set fraction of cells in a CDS and perform DEG test.")
parser$add_argument('cds', help='CDS to use.')
parser$add_argument('output_file', help='Text file to output DEG results to.')
parser$add_argument('--swap_rate', type="double", help='Swap rate applied prior to DEG.')
parser$add_argument('--treatment', help='Condition to subset from the CDS (treatment column must be present in pData).')
parser$add_argument('--seed', type='integer', help='Random seed to use in sampling.')
parser$add_argument('--use_as_is', action='store_true', help='Set to turn off dispersion and size factor estimates and instead use the CDS object as passed to script. Useful for very large datasets where dispersion estimates can cause memory spikes.')
args = parser$parse_args()


# Load in and process data
aggregated_cds = readRDS(args$cds)

args$swap_rate = as.numeric(args$swap_rate)
args$seed = as.numeric(args$seed)

aggregated_cds = aggregated_cds[, row.names(subset(pData(aggregated_cds), !is.na(gene) & treatment == args$treatment & guide_count == 1))]

if (! args$use_as_is) {
	aggregated_cds = estimateSizeFactors(aggregated_cds)
	aggregated_cds = estimateDispersions(aggregated_cds)
}

# Simulate swaps
## Keep half the data for each mutant

# If the dataset is very large, downsample to 25K cells just so
# it doesn't take forever to run
set.seed(0)
if (nrow(pData(aggregated_cds)) > 25000) { 
    cells_to_keep = pData(aggregated_cds) %>%
        group_by(gene) %>%
        sample_frac(6000 / nrow(pData(aggregated_cds))) %>%
        ungroup()

    aggregated_cds = aggregated_cds[,cells_to_keep$cell]
}

set.seed(args$seed)
pData(aggregated_cds)$id = row.names(pData(aggregated_cds))

cells_to_keep = pData(aggregated_cds) %>%
	group_by(gene) %>%
	sample_frac(1 - args$swap_rate) %>%
	ungroup()

cells_to_shuffle = colnames(aggregated_cds)[! colnames(aggregated_cds) %in% cells_to_keep$id]
pData(aggregated_cds)[cells_to_shuffle, "gene"] = sample(pData(aggregated_cds)[cells_to_shuffle, "gene"])

# Do the DE test and write out the results
expressed_genes <- row.names(fData(aggregated_cds)[rowSums(exprs(aggregated_cds) > 0) > 50 ,])
de_results = differentialGeneTest(aggregated_cds[expressed_genes, ], fullModelFormulaStr="~gene", cores=8)
de_results$swap_rate = args$swap_rate
de_results$seed = args$seed
write.table(de_results, file=args$output_file, quote=FALSE, row.names=FALSE, sep="\t")
