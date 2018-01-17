library(monocle)
library(ggplot2)
library(stringr)
source('helper_functions.R')

# Get dox treated cells
aggregated_cds = readRDS('temp_data/aggregated_cds.rds')
aggregated_cds.dox = aggregated_cds[, aggregated_cds$treatment == 'dox_100nm']

# Estimate size factors and dispersions
aggregated_cds.dox = estimateSizeFactors(aggregated_cds.dox)
aggregated_cds.dox = estimateDispersions(aggregated_cds.dox)

# Get expressed genes and filter
expressed_genes = row.names(exprs(aggregated_cds.dox)[rowSums(exprs(aggregated_cds.dox) > 0) > 50, ])

aggregated_cds.dox = aggregated_cds.dox[expressed_genes, ]

# Now do TSNE
aggregated_cds.dox = reduceDimension(aggregated_cds.dox, reduction_method='tSNE', num_dim=20, norm_method='log', pseudo_expr=1, max_components=2)

# Make a TSNE plot highlighting relevant groups of cells
tsne_coordinates = as.data.frame(t(reducedDimA(aggregated_cds.dox)))
colnames(tsne_coordinates) = c('tsne_1', 'tsne_2')

tsne_coordinates$assignment = aggregated_cds.dox$gene

tsne_coordinates$assignment[str_detect(tsne_coordinates$assignment, 'TP53')] = 'TP53'

tsne_coordinates$assignment[tsne_coordinates$assignment == 'NONTARGETING'] = 'NTC'

tsne_coordinates$assignment[is.na(tsne_coordinates$assignment)] = 'unassigned'

tsne_coordinates$assignment[! tsne_coordinates$assignment %in% c("TP53", 'unassigned', 'NTC')] = 'other'

custom_colors = c('TP53'='red', other='#d3d3d3', 'NTC'='#53585F', 'unassigned'='#4daf4a')

ggplot(tsne_coordinates, aes(tsne_1, tsne_2)) +
	geom_point(aes(color=assignment), size=0.3, alpha=0.6) +
	theme_cfg() +
	scale_color_manual(values=custom_colors, guide=F) +
	xlab('TSNE 1') +
	ylab('TSNE 2') +
	ggsave('figures/tsne_by_assignment.png', height=3.5, width=3.5)

# Cluster cells and make marker plots
aggregated_cds.dox = clusterCells(aggregated_cds.dox)

# Cluster 8 is TP53 
plot_cell_clusters(aggregated_cds.dox) + 
	theme_cfg() +
	ggsave('diagnostic_plots/tsne_by_cluster.png')

# Now make marker plots
TP53_markers = c("CDKN1A","TP53I3")

TP53_main_text_ids = row.names(subset(fData(aggregated_cds.dox), gene_short_name %in% TP53_markers))

TP53_target_cds_subset <- aggregated_cds.dox[c(TP53_main_text_ids),pData(aggregated_cds.dox)$gene %in% c("TP53","NONTARGETING")]

TP53_target_cds_subset$gene[TP53_target_cds_subset$gene == 'NONTARGETING'] = 'NTC'

# Reassign NTC labels for cases where not in informative clusters
enriched_clusters = c(8)
pData(TP53_target_cds_subset)$label <- ifelse(pData(TP53_target_cds_subset)$Cluster %in% enriched_clusters, "TP53 Cluster", "Other")

TP53_target_cds_subset$interaction = interaction(TP53_target_cds_subset$label, TP53_target_cds_subset$gene, sep=' ')

# Relevel for plot
TP53_target_cds_subset$interaction = factor(TP53_target_cds_subset$interaction, levels=c("Other NTC", "All TP53", 'Other TP53', 'TP53 Cluster TP53'))

# Labels for plot ordering
#tp53_labels = c("NTC", expression(paste("All", italic(TP53))), expression(paste("Other", italic(TP53))), expression(paste(italic(TP53), "Cluster")))

tp53_labels = c("NTC", expression(paste("All", ' TP53')), expression(paste("Other", ' TP53')), expression(paste('TP53', " Cluster")))

tp53_colors = c("Other NTC"="#53585F", "All TP53"="#feb24c", "TP53 Cluster TP53"="red", " NTC"="#53585F", "Other TP53"="#ffeda0")


# Note show_combined shows an All TP53 category in addition to the other groups
plot_genes_violin(TP53_target_cds_subset[TP53_main_text_ids, ], grouping = "interaction", color_by = "interaction", ncol = 3, plot_trend = T, min_expr = 0.1, relative_expr = TRUE, log_scale = TRUE, show_combined = c("TP53")) + 
    scale_fill_manual(values=tp53_colors) +
    guides(fill=FALSE) +
    scale_x_discrete(labels=tp53_labels, name="Gene KO") +
    theme_cfg() +
    #theme(strip.text = element_text(face = "italic"), axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    ggsave('figures/tsne_by_assignment_tp53_markers.png', width = 5, height = 3)


