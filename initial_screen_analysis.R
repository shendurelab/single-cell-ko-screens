library(monocle)

library(dplyr)
source('helper_functions.R')

# Do TSNE on each set of cells separately
do_tsne = function(cds) {
	cds = estimateSizeFactors(cds)
	cds = estimateDispersions(cds)

	cds = cds[expressed_genes, !is.na(pData(cds)$gene)]

	cds = reduceDimension(cds, reduction_method='tSNE', max_components = 2, num_dim=10, norm_method='log', pseudo_expr=1)
	return(cds)
}

plot_tsne = function(cds, color_by='gene') {
	tSNE_dim_coords = reducedDimA(cds)
	data_df = data.frame(t(tSNE_dim_coords[c(1, 2), ]))
	data_df[, color_by] = pData(cds)[, color_by]

	ggplot(data_df, aes(X1, X2)) +
		geom_point(aes_string(color=color_by)) +
		theme_classic()
}

# Load in data from initial screen
initial_screen_data = readRDS('temp_data/aggregated_cds.initial_screens.rds')

# Subset to just the dox treated experiments
initial_screen_data = initial_screen_data[, pData(initial_screen_data)$treatment != "mock"]

# Get list of expressed genes
expressed_genes = row.names(initial_screen_data)[rowSums(exprs(initial_screen_data) > 0) > 50]

# Estimate size factors
initial_screen_data = estimateSizeFactors(initial_screen_data)

# Subset only to the set of genotypes reliably detected in both screens
genotypes_to_use = c("CBFB", "NCOR1", "PTEN", "TP53")
initial_screen_data = initial_screen_data[, initial_screen_data$gene %in% genotypes_to_use]

# Subset out each screen and estimate dispersions
arrayed = initial_screen_data[, pData(initial_screen_data)$condition == "arrayed_dox_500nm"]
pooled = initial_screen_data[, pData(initial_screen_data)$condition == "pooled_dox_500nm"]

arrayed = estimateDispersions(arrayed)
pooled = estimateDispersions(pooled)

# Markers for TP53 KOs show clear patterns
p53_markers = row.names(subset(fData(pooled), gene_short_name %in% c("CDKN1A", "TP53I3")))

plot_genes_violin(arrayed[p53_markers, ], color_by='gene', grouping='gene', log_scale=T, plot_trend=T, min_expr = 0.1, relative_expr = TRUE) +
	theme_cfg() +
	scale_fill_manual(values=c("CBFB"='#d3d3d3', "NCOR1"='#d3d3d3', "PTEN"='#d3d3d3', "TP53"='#d3d3d3')) +
	scale_y_log10(limits=c(0.1, 15)) +
	ylab('expression') +
	xlab('target') +
	ggsave('supplemental_figures/initial_screens.arrayed.markers.pdf', width=4, height=4)

plot_genes_violin(pooled[p53_markers, ], color_by='gene', grouping='gene', log_scale=T, plot_trend=T, min_expr = 0.1, relative_expr = TRUE) +
	theme_cfg() +
	scale_fill_manual(values=c("CBFB"='red', "NCOR1"='red', "PTEN"='red', "TP53"='red')) +
	ylab('expression') +
	xlab('target') +
	ggsave('supplemental_figures/initial_screens.pooled.markers.pdf', width=4, height=4)

# Load in results from matched DEG tests
deg_results = do.call(rbind, lapply(Sys.glob('temp_data/initial_screen_analysis_deg/*.rds'), readRDS))

deg_counts = deg_results %>%
	dplyr::filter(qval < 0.05) %>%
	dplyr::group_by(sample, seed) %>%
	dplyr::summarize(count = n())

deg_counts$sample = ifelse(deg_counts$sample == 1, 'arrayed', 'pooled')

ggplot(deg_counts, aes(sample, count)) +
	geom_boxplot(fill='#d3d3d3', color="#CCCCCC", alpha=0.2) +
	geom_point(aes(fill=sample), color="black", pch=21, size=2.5, position = position_jitter(w = 0.25, h = 0)) +
	scale_fill_manual(values=c("arrayed"="#d3d3d3", "pooled"="red")) +
	theme_cfg() +
	ylab('DEG count') +
	xlab('experiment') +
	ggsave('supplemental_figures/initial_screens.matched_deg.pdf', width=4, height=4)

