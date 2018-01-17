
###### Load packages ######
# Load necessary packages for single cell RNA-Seq analysis including packages for downstream Gene Ontology Analysis
library(monocle)
library(plyr)
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)

library(piano)
library(tidyr)
library(viridis)
library(ggrepel)
source('helper_functions.R')

# Colors for plots
cluster_colors=c("1"="#e41a1c", "2"="#377eb8", "3"="#4daf4a", "4"="#984ea3", "5"="#ff7f00", "6"="#a65628", "NONE"="black")

# Load in data
informative_filelist <- Sys.glob("temp_data/pairwise_deg/informative/singles/cropseq/*.txt")
all_pairwise_informative_deg_results_gene = do.call(rbind, lapply(informative_filelist, read.delim))
all_pairwise_informative_deg_results_gene <- all_pairwise_informative_deg_results_gene[,-1]

dox_100nm_cds <- readRDS("temp_data/barcode_enrichment/dox_100nm_cds.rds")

# Preprocess RNA-seq data
dox_100nm_expressed_genes <- row.names(fData(dox_100nm_cds)[rowSums(exprs(dox_100nm_cds) > 0) > 50 ,])
dox_100nm_informative_cds <- dox_100nm_cds[,pData(dox_100nm_cds)$informative == TRUE]

preprocess_cds <- function(cds){
    cds <- detectGenes(cds, min_expr = 0.5)
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
}

dox_100nm_informative_cds <- preprocess_cds(dox_100nm_informative_cds)

# Get the cluster-target associations where enrichment was found
dox.cluster_target_pairs = pData(dox_100nm_informative_cds) %>% filter(informative == TRUE & gene != "NONTARGETING") %>% select(gene, Cluster) %>% unique()

# Overlap between DEGs and TP53 set of DEGs for cells where their target-cluster pair was found to be enriched
dox_singles_deg_results <- all_pairwise_informative_deg_results_gene[all_pairwise_informative_deg_results_gene$column_value == "dox_100nm",]

dox_singles_significant_degs <- dox_singles_deg_results[dox_singles_deg_results$qval < 0.05,]

dox_singles_significant_genes = unique(dox_singles_significant_degs$id)

dox.tp53_genes = unique(subset(dox_singles_significant_degs, target == 'TP53')$id)


dox.tp53_deg_overlap = dox_singles_significant_degs %>%
                        group_by(target) %>%
                        summarize(overlap=sum(id %in% dox.tp53_genes), not_overlap=sum(! id %in% dox.tp53_genes), proportion=overlap / (overlap + not_overlap)) %>%
                        arrange(proportion) %>%
                        mutate(target = factor(target, levels=target))

dox.tp53_deg_overlap = merge(dox.tp53_deg_overlap, dox.cluster_target_pairs, by.x='target', by.y='gene')

## Dox
ggplot(dox.tp53_deg_overlap, aes(target, proportion)) +
        geom_bar(stat='identity', color='black', aes(fill=Cluster)) +
        theme_cfg() +
        xlab('target') +
        ylab('proportion\nTP53 DEG Overlap') +
        ylim(0, 1) +
        coord_flip() +
        scale_fill_manual(values=cluster_colors, guide=FALSE) +
        ggsave('supplemental_figures/dox.tp53_deg_overlaps.pdf', height=3.5, width=2.5)

# Now PC analysis
getMeanTargetExpression <- function(cds, DEGs){

	gene_list <- sort(unique(pData(cds)$gene))

	Mean_expression_byGene <- list()
	Mean_expression <- list()

	i <-  1

	for (gene in gene_list){

	    cds_subset <- cds[DEGs,pData(cds)$gene == gene]
	    cds_subset_expr <- exprs(cds_subset)
	    cds_subset_expr <- Matrix::t(Matrix::t(cds_subset_expr)/sizeFactors(cds_subset))
	    cds_subset_Mean_expr <- as.matrix(rowMeans(cds_subset_expr), ncol = 1)

	    colnames(cds_subset_Mean_expr) <- as.factor(gene)

    	Mean_expression_byGene[[i]] <- cds_subset_Mean_expr
   	 names(Mean_expression_byGene[[i]]) <- as.factor(gene)
   	 i <- i + 1

	}

	Mean_expression <- do.call("cbind",Mean_expression_byGene)

    	return(Mean_expression)
  }

dox_100nm_normalized_matrix <- getMeanTargetExpression(dox_100nm_informative_cds, dox_singles_significant_genes)

getScaledMatrix <- function(matrix){
    scaled_matrix <- scale(t(scale(t(matrix))))
    scaled_matrix[scaled_matrix > 3] <- 3
    scaled_matrix[scaled_matrix < -3] <- -3
    scaled_matrix[is.na(scaled_matrix)] <- 0

    return(scaled_matrix)
}

dox_100nm_normalized_matrix_scaled = getScaledMatrix(dox_100nm_normalized_matrix)
dox_100nm_pca_matrix <- t(dox_100nm_normalized_matrix_scaled)
dox_pca <- prcomp(dox_100nm_pca_matrix, center = T, scale. = T)

dox_target_loadings <- data.frame(dox_pca$x,
.names = row.names(dox_pca$x))

dox_gene_loadings <- data.frame(dox_pca$rotation,
.names = row.names(dox_pca$rotation))

dox_target_loadings = merge(dox_target_loadings, dox.cluster_target_pairs, by.x="row.names", by.y="gene", all.x=T)
dox_target_loadings$Cluster = as.character(dox_target_loadings$Cluster)
dox_target_loadings$Cluster[is.na(dox_target_loadings$Cluster)] = "NONE"
dox_target_loadings$.names[dox_target_loadings$.names == 'NONTARGETING'] = 'NTC'

ggplot(dox_target_loadings, aes(x = PC1)) +
    geom_histogram(color='black', aes(fill = Cluster)) +
    theme_cfg() +
    scale_fill_manual(values=cluster_colors, guide=FALSE) +
    ggsave("supplemental_figures/dox.informative.pc1_histogram.png", width = 3.5, height = 1.75)


# Load necessary functions
source("jose_helper_functions/loadGSCSafe.R")
source('jose_helper_functions/plot_gsea_go.R')
source("jose_helper_functions/GSA_hyper_helper_functions.R")

## Load Gene Set Collections
hallmarksGSC <- loadGSCSafe(file="data/gmt_files/h.all.v6.0.symbols.gmt")

getGeneLoadingGSAhyperResults <- function(cds, expressed_genes, pca_res, gsc, PC_list, cutoff){

    Ensembl_GSAlist <- as.matrix(fData(cds[expressed_genes])$gene_short_name)
    rownames(Ensembl_GSAlist)<-row.names(fData(cds[expressed_genes]))
    colnames(Ensembl_GSAlist) <- c("gene_short_name")
    Ensembl_GSAlist<-Ensembl_GSAlist[,1]
    Ensembl_GSAlist<-unique(toupper(Ensembl_GSAlist))
    length(Ensembl_GSAlist)

    GSAhyper.list <- list()

    positive_gene_loadings_GSAhyper <- list()
    negative_gene_loadings_GSAhyper <- list()

    aload <- abs(pca_res$rotation)
    sweep(aload, 2, colSums(aload), "/")

    gene_loadings <- data.frame(pca_res$rotation,
    .names = row.names(pca_res$rotation))

    print(PC_list)

    for(PC in PC_list){

    top_positive_gene_loadings <- row.names(gene_loadings[gene_loadings[,PC] > cutoff[1],])
    top_negative_gene_loadings <- row.names(gene_loadings[gene_loadings[,PC] < cutoff[2],])

        print(length(top_positive_gene_loadings))
        print(length(top_negative_gene_loadings))

    positive_genes <-  as.matrix(fData(cds[top_positive_gene_loadings])$gene_short_name)
    rownames(positive_genes)<-row.names(fData(cds[top_positive_gene_loadings]))
    colnames(positive_genes) <- c("gene_short_name")
    positive_genes<-positive_genes[,1]
    positive_genes<-unique(toupper(positive_genes))
    length(positive_genes)

    negative_genes <-  as.matrix(fData(cds[top_negative_gene_loadings])$gene_short_name)
    rownames(negative_genes)<-row.names(fData(cds[top_negative_gene_loadings]))
    colnames(negative_genes) <- c("gene_short_name")
    negative_genes<-negative_genes[,1]
    negative_genes<-unique(toupper(negative_genes))
    length(negative_genes)

    positive_gene_loadings_GSAhyper[[PC]] <- collect_gsa_hyper_results_All(positive_genes,
                                                                           Ensembl_GSAlist,
                                                                           gsc)

    negative_gene_loadings_GSAhyper[[PC]] <- collect_gsa_hyper_results_All(negative_genes,
                                                                           Ensembl_GSAlist,
                                                                           gsc)
        }

    GSAhyper.list[["positive_loadings"]] <- positive_gene_loadings_GSAhyper
    GSAhyper.list[["negative_loadings"]] <- negative_gene_loadings_GSAhyper

    return(GSAhyper.list)
}

Dox_geneLoadings_GSAhyper <- getGeneLoadingGSAhyperResults(dox_100nm_cds, dox_100nm_expressed_genes, dox_pca,
                                                            hallmarksGSC, c("PC1"), c(0.02,-0.02))

gsea_plot <- function(GSAhyper_list, qval_cutoff, pattern, sample, gsc, PC){

    GSAhyper_df.positive <- as.data.frame(GSAhyper_list[['positive_loadings']][[PC]]$diffExp$p.adj)
    GSAhyper_df.positive$loading= 'positive'

    GSAhyper_df.positive$gene_set = row.names(GSAhyper_df.positive)
    colnames(GSAhyper_df.positive) <- c("qval", 'loading', 'gene_set')

    GSAhyper_df.negative <- as.data.frame(GSAhyper_list[['negative_loadings']][[PC]]$diffExp$p.adj)
    GSAhyper_df.negative$loading= 'negative'
    GSAhyper_df.negative$gene_set = row.names(GSAhyper_df.negative)
    colnames(GSAhyper_df.negative) <- c("qval", 'loading', 'gene_set')

    GSAhyper_df = rbind(GSAhyper_df.positive, GSAhyper_df.negative)




    if(is.null(pattern) == FALSE){
        GSAhyper_df$gene_set <- str_replace(string = GSAhyper_df$gene_set, pattern = pattern, replace = "")
    }

    GSAhyper_df$gene_set = str_replace_all(GSAhyper_df$gene_set, '_', ' ')

    GSAhyper_df_cutoff <- GSAhyper_df %>% filter(qval < qval_cutoff) %>% arrange(desc(qval)) %>%
    mutate(gene_set = factor(gene_set, levels = gene_set))

    ggplot(GSAhyper_df_cutoff, aes(x = gene_set, y = -log10(qval))) +
        geom_bar(stat = "identity", color='black', aes(fill = loading)) +
        coord_flip() +
        facet_wrap(~loading, ncol=1, scales='free') +
        xlab('gene set')
}

gsea_plot(Dox_geneLoadings_GSAhyper, 0.01, pattern = "HALLMARK_", sample = "Dox_100nm", gsc = "hallmarksGSC", PC='PC1') +
    theme_cfg() +
    scale_fill_manual(values=c('positive'="#53585F", 'negative'='#53585F'), guide=FALSE) +
    ggsave('supplemental_figures/dox.pc1.loadings.gsea.pdf', height=3.5, width=4.5)
