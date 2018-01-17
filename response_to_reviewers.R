#############################################################################
# To help provide more confidence in the pooled vs. arrayed experiment, make 
# a TSNE plot of the two experiments colored by TP53 vs. other as well as by some
# basic markers of mitosis
#############################################################################
library(monocle)
library(dplyr)
source('helper_functions.R') 

initial_screens = readRDS('temp_data/aggregated_cds.initial_screens.rds')  

cds.arrayed = initial_screens[, initial_screens$condition == 'arrayed_dox_500nm']
cds.pooled = initial_screens[, initial_screens$condition == 'pooled_dox_500nm']

do_tsne = function(cds, num_dim=20, num_cells_exprs=50) {
	cds = estimateSizeFactors(cds)
	cds = estimateDispersions(cds)
	expressed_genes = row.names(exprs(cds)[rowSums(exprs(cds) > 0) > num_cells_exprs, ])
	cds = cds[expressed_genes, ]
	cds = reduceDimension(cds, reduction_method='tSNE', norm_method='log', pseudo_expr=1, max_components=2, num_dim=num_dim)
	cds = clusterCells(cds)
	return(cds)
}

cds.arrayed = do_tsne(cds.arrayed, num_dim=12, num_cells_exprs=50)
cds.pooled = do_tsne(cds.pooled, num_dim=12, num_cells_exprs=50)

# Highlight a few genes within each tsne
ccnb2 = 'ENSG00000157456'
tp53i3 = 'ENSG00000115129'

plot_genes = function(cds, gene) {
	reduced_dim = as.data.frame(t(reducedDimA(cds)))
	colnames(reduced_dim) = c('tsne1', 'tsne2')
	reduced_dim$gene = exprs(cds)[gene, ] / cds$Size_Factor
	reduced_dim$Cluster = pData(cds)$Cluster

	reduced_dim = reduced_dim %>%
		dplyr::group_by(Cluster) %>%
		dplyr::mutate(score = mean(gene)) %>%
		dplyr::ungroup()
	ggplot(reduced_dim, aes(tsne1, tsne2)) +
		geom_point(aes(color=score), size=0.5) +
		theme_cfg() + 
		scale_color_viridis()
}
# Highlight the TP53 labeled cells from each experiment
plot_genes(cds.arrayed, ccnb2) + ggsave('reviewer_figures/arrayed_screen.tsne.ccnb2.png', height=4, width=4)
plot_genes(cds.pooled, ccnb2) + ggsave('reviewer_figures/pooled_screen.tsne.ccnb2.png', height=4, width=4)

plot_genes(cds.arrayed, tp53i3) + ggsave('reviewer_figures/arrayed_screen.tsne.tp53i3.png', height=4, width=4)
plot_genes(cds.pooled, tp53i3) + ggsave('reviewer_figures/pooled_screen.tsne.tp53i3.png', height=4, width=4)


plot_genotype = function(cds, selected_gene) {
	reduced_dim = as.data.frame(t(reducedDimA(cds)))
	colnames(reduced_dim) = c('tsne1', 'tsne2')
	reduced_dim$gene = pData(cds)$gene

	reduced_dim$gene[is.na(reduced_dim$gene)] = 'unassigned'
	reduced_dim$gene[grepl(selected_gene, reduced_dim$gene)] = selected_gene
	reduced_dim$gene[!grepl(selected_gene, reduced_dim$gene) & !reduced_dim$gene == 'unassigned'] = 'other'

	custom_colors = c('TP53'='red', other='#d3d3d3', 'NTC'='#53585F', 'unassigned'='#9758A3')

	ggplot(reduced_dim, aes(tsne1, tsne2)) +
		geom_point(aes(color=gene), size=0.5) +
		theme_cfg() + 
		scale_color_manual(values=custom_colors, guide=FALSE)
}

plot_genotype(cds.arrayed, 'TP53') + ggsave('reviewer_figures/arrayed_screen.tsne.genotypes.png', height=4, width=4)
plot_genotype(cds.pooled, 'TP53') + ggsave('reviewer_figures/pooled_screen.tsne.genotypes.png', height=4, width=4)

# quantify the proportion of 'other' genotypes in the TP53 cluster for the pooled screen to show close to 50:50
plot_cell_clusters(cds.arrayed) + ggsave('reviewer_figures/arrayed_screen.tsne.clusters.png') 
plot_cell_clusters(cds.pooled) + ggsave('reviewer_figures/pooled_screen.tsne.clusters.png')

# cluster 4 is the tp53 cluster so quantify within that
pData(cds.pooled) %>% 
	dplyr::filter(Cluster == 4) %>%
	dplyr::group_by(grepl('TP53', gene)) %>%
	dplyr::summarize(count = n())

pData(cds.pooled) %>% 
	dplyr::filter(Cluster == 4) %>%
	dplyr::group_by(is.na(gene)) %>%
	dplyr::summarize(count = n())

# cluster 3 is the TP53 cluster so quantify within that
pData(cds.arrayed) %>% 
	dplyr::filter(Cluster == 3) %>%
	dplyr::group_by(grepl('TP53', gene)) %>%
	dplyr::summarize(count = n())

pData(cds.arrayed) %>% 
	dplyr::filter(Cluster == 3) %>%
	dplyr::group_by(is.na(gene)) %>%
	dplyr::summarize(count = n())

# we observe 41% of assigned cells are TP53 in pooled but 99.4% in arrayed

# Confirm that there are no genes that show significant deviation from the distribution of cell counts over the different clusters (confirms that there is not some non-random reason for this 42% that could be due to other targets)
pval_list = c()
gene_list = unique(subset(pData(cds.pooled), guide_count == 1)$gene)
for (ko in gene_list) {
	# Distribution of all cells over all clusters
	null.dist = table(pData(cds.pooled) %>% select(Cluster))

	# Distribution of the gene discounting anything with TP53 included unless testing TP53
	gene.dist = table(pData(cds.pooled)[pData(cds.pooled)[, ko] & (ko == 'TP53' | !pData(cds.pooled)[, 'TP53']),] %>% select(Cluster))

	pval_list = c(pval_list, chisq.test(p=null.dist, rescale.p = TRUE, x=gene.dist)$p.val)
}

gene_list[p.adjust(pval_list, method = 'BH') < 0.05] # only one gene (TP53) is signficant here

####################################################################################
# One reviewer asked about cells that contain multiple guides (how many vs. singles)
# they specifically asked about the pooled experiment from the pLGB-scKO vector
####################################################################################
# counts of cells with various numbers of guides in this experiment
dim(subset(pData(cds.pooled), guide_count == 1 & !is.na(guide_count)))
dim(subset(pData(cds.pooled), guide_count == 2 & !is.na(guide_count)))
dim(subset(pData(cds.pooled), guide_count == 3 & !is.na(guide_count)))
dim(subset(pData(cds.pooled), guide_count == 4 & !is.na(guide_count)))
dim(subset(pData(cds.pooled), guide_count > 4 & !is.na(guide_count)))
dim(subset(pData(cds.pooled), is.na(guide_count)))

# proportion of cells with single guides that we observe in this experiment
dim(subset(pData(cds.pooled), guide_count == 1 | is.na(guide_count))) / dim(subset(pData(cds.pooled), guide_count > 1 & !is.na(guide_count)))

# This is the expected (80%) of cells with a single guide in a selected populaton at our approximate MOI
dpois(1, lambda=0.45)/(dpois(2, lambda=0.45) + dpois(3, lambda=0.45) + dpois(4, lambda=0.45) + dpois(1, lambda=0.45))

#############################################################################
# For the response, get the number of guides per target for each screen
# (this confirmed that one guide was used per target in arrayed and that those same
# guides were also included in the pooled screen)
# All the same guides from the pooled screen were present in the CROP-seq screen
#############################################################################
# arrayed screen guides
pData(initial_screens) %>%
	filter(! is.na(gene) & guide_count == 1 & screen == 'arrayed') %>% 
	select(gene, barcode) %>%
	unique()

# pooled screen guides
pData(initial_screens) %>%
	filter(! is.na(gene) & guide_count == 1 & screen == 'pooled') %>% 
	select(gene, barcode) %>%
	unique() %>%
	arrange(gene)

#############################################################################
# Confirm that we see more cells assigned in the true TP53 null population
# in the assignments with no targeted enrichment than in the enriched library
# to support theory for unexpected result in doxorubicin treated sample DEG without
# targeted enrichment
#############################################################################
all_assignments_cds = readRDS('temp_data/aggregated_cds.rds')  
no_enrichment_cds = readRDS('temp_data/aggregated_cds.no_enrichment.rds')
all_assignments_cds = all_assignments_cds[, all_assignments_cds$treatment == 'dox_100nm']
no_enrichment_cds = no_enrichment_cds[, no_enrichment_cds$treatment == 'dox_100nm'] 


all_assignments_cds = do_tsne(all_assignments_cds, num_dim=20, num_cells_exprs=50)
no_enrichment_cds = do_tsne(no_enrichment_cds, num_dim=20, num_cells_exprs=50)

plot_cell_clusters(all_assignments_cds) + ggsave('reviewer_figures/cropseq_screen.enrichment.tsne.clusters.png') # cluster 8 is the cluster

no_enrichment_cds = clusterCells(no_enrichment_cds, delta_threshold=5, rho_threshold=40)
plot_cell_clusters(no_enrichment_cds) + ggsave('reviewer_figures/cropseq_screen.no_enrichment.tsne.clusters.png') # cluster 9 is the cluster

plot_genotype(all_assignments_cds, 'TP53') + ggsave('reviewer_figures/cropseq_screen.enrichment.tsne.genotypes.png', height=4, width=4)
plot_genotype(no_enrichment_cds, 'TP53') + ggsave('reviewer_figures/cropseq_screen.no_enrichment.tsne.genotypes.png', height=4, width=4)

# cluster 8 is the tp53 cluster so quantify within that
pData(all_assignments_cds) %>% 
dplyr::group_by(Cluster == 8, is.na(gene)) %>%
	dplyr::summarize(count = n())

# cluster 9 is the TP53 cluster so quantify within that
pData(no_enrichment_cds) %>% 
	dplyr::group_by(Cluster == 9, is.na(gene)) %>%
	dplyr::summarize(count = n())

