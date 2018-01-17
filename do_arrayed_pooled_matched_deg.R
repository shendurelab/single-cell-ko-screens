library(monocle)

library(dplyr)
library(argparse)

parser = argparse::ArgumentParser(description='Script to perform matched DE tests for the the arrayed and pooled experiments.')
parser$add_argument('cds', help='RDS file with CDS to use')
parser$add_argument('output_file', help='Output RDS file with DEG test results.')
parser$add_argument('--seed', default=0, type="integer", help='Seed to use for random sampling.')
args = parser$parse_args()

do_matched_deg = function(cds1, cds2, seed, cores=detectCores()) {
	sampled_cds1_cells = c()
	sampled_cds2_cells = c()

	for (gene in unique(cds1$gene)) {
		cds1_count = sum(cds1$gene == gene)
		cds2_count = sum(cds2$gene == gene)
		sampled_count = min(cds1_count, cds2_count)

		set.seed(seed)
		cds1_cells = sample(colnames(cds1[, cds1$gene == gene]), sampled_count)
		cds2_cells = sample(colnames(cds2[, cds2$gene == gene]), sampled_count)
		sampled_cds1_cells = c(sampled_cds1_cells, cds1_cells)
		sampled_cds2_cells = c(sampled_cds2_cells, cds2_cells)
	}
	
	closeAllConnections()
	cds1_deg_results = differentialGeneTest(cds1[, sampled_cds1_cells], fullModelFormulaStr="~gene", cores=cores)
	cds1_deg_results$sample = '1'
	
	closeAllConnections()
	cds2_deg_results = differentialGeneTest(cds2[, sampled_cds2_cells], fullModelFormulaStr="~gene", cores=cores)
	cds2_deg_results$sample = '2'

	all_deg_results = rbind(cds1_deg_results, cds2_deg_results)
	all_deg_results$seed = seed
	return(all_deg_results)
}


# Load in data from initial screen
initial_screen_data = readRDS('temp_data/aggregated_cds.initial_screens.rds')

# Subset to just the dox treated experiments
initial_screen_data = initial_screen_data[, pData(initial_screen_data)$treatment != "mock"]

# Get list of expressed genes
expressed_genes = row.names(initial_screen_data)[rowSums(exprs(initial_screen_data) > 0) > 50]

# Estimate size factors
initial_screen_data = estimateSizeFactors(initial_screen_data)
initial_screen_data = estimateDispersions(initial_screen_data)

# Subset only to the set of genotypes reliably detected in both screens and estimate dispersions
genotypes_to_use = c("CBFB", "NCOR1", "PTEN", "TP53")
initial_screen_data = initial_screen_data[, initial_screen_data$gene %in% genotypes_to_use]

# Subset out each screen
arrayed = initial_screen_data[, pData(initial_screen_data)$condition == "arrayed_dox_500nm"]
pooled = initial_screen_data[, pData(initial_screen_data)$condition == "pooled_dox_500nm"]

deg_results = do_matched_deg(arrayed[expressed_genes, ], pooled[expressed_genes, ], seed=args$seed, cores=8)
saveRDS(deg_results, file=args$output_file)
