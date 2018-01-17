
###### Load packages ######
# Load necessary packages for single cell RNA-Seq analysis including packages for downstream Gene Ontology Analysis
#library(monocle)
suppressPackageStartupMessages({library(devtools)
library(monocle)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(colorRamps)
library(tidyr)
library(stringr)})
source('helper_functions.R')

gene_barcode_association <- read.table("data/jose_gene_barcode_association.txt", sep = "\t")

initial.target.level.chisq.qval <- readRDS("temp_data/barcode_enrichment/initial.target.level.chisq.qval.rds")
initial.guide.level.chisq.qval <- readRDS("temp_data/barcode_enrichment/initial.guide.level.chisq.qval.rds")

# Load results that contain TSNE for each
cds.v1.mock = readRDS("temp_data/barcode_enrichment/mock_cds.rds")
cds.v1.dox.100nm = readRDS("temp_data/barcode_enrichment/dox_100nm_cds.rds")

# TSNE figures
## mock
cluster_colors=c("1"="#e41a1c", "2"="#377eb8", "3"="#4daf4a", "4"="#984ea3", "5"="#ff7f00", "6"="#a65628")
tsne_height=3
tsne_width=3
tsne_cell_size = 0.25

plot_cell_clusters(cds.v1.mock, x = 1, y = 2, cell_size = tsne_cell_size, color_by = "Cluster") +
        theme_cfg() +
	scale_color_manual(values=cluster_colors) +
	xlab('TSNE 1') +
	ylab('TSNE 2') +
	theme(legend.position = "none") +
        ggsave('supplemental_figures/mock.jose_enrichment.tsne.1_2.png', height=tsne_height, width=tsne_width)

plot_cell_clusters(cds.v1.mock, x = 1, y = 3, cell_size = tsne_cell_size, color_by = "Cluster") +
        theme_cfg() +
        scale_color_manual(values=cluster_colors) +
        xlab('TSNE 1') +
        ylab('TSNE 3') +
        theme(legend.position = "none") +
        ggsave('supplemental_figures/mock.jose_enrichment.tsne.1_3.png', height=tsne_height, width=tsne_width)

plot_cell_clusters(cds.v1.mock, x = 2, y = 3, cell_size = tsne_cell_size, color_by = "Cluster") +
        theme_cfg() +
        scale_color_manual(values=cluster_colors) +
        xlab('TSNE 2') +
        ylab('TSNE 3') +
        theme(legend.position = "none") +
        ggsave('supplemental_figures/mock.jose_enrichment.tsne.2_3.png', height=tsne_height, width=tsne_width)

## Dox
plot_cell_clusters(cds.v1.dox.100nm, x = 1, y = 2, cell_size = tsne_cell_size, color_by = "Cluster") +
        theme_cfg() +
        scale_color_manual(values=cluster_colors) +
        xlab('TSNE 1') +
        ylab('TSNE 2') +
        theme(legend.position = "none") +
        ggsave('supplemental_figures/dox.jose_enrichment.tsne.1_2.png', height=tsne_height, width=tsne_width)

plot_cell_clusters(cds.v1.dox.100nm, x = 1, y = 3, cell_size = tsne_cell_size, color_by = "Cluster") +
        theme_cfg() +
        scale_color_manual(values=cluster_colors) +
        xlab('TSNE 1') +
        ylab('TSNE 3') +
        theme(legend.position = "none") +
        ggsave('supplemental_figures/dox.jose_enrichment.tsne.1_3.png', height=tsne_height, width=tsne_width)

plot_cell_clusters(cds.v1.dox.100nm, x = 2, y = 3, cell_size = tsne_cell_size, color_by = "Cluster") +
        theme_cfg() +
        scale_color_manual(values=cluster_colors) +
        xlab('TSNE 2') +
        ylab('TSNE 3') +
        theme(legend.position = "none") +
        ggsave('supplemental_figures/dox.jose_enrichment.tsne.2_3.png', height=tsne_height, width=tsne_width)


# Guide enrichment summary figures
mock_target_chisq.qval <- as.data.frame(initial.target.level.chisq.qval[["V1 mock"]])
mock_target_chisq.qval$barcode <- row.names(mock_target_chisq.qval)
colnames(mock_target_chisq.qval) <- c("qval","target")

mock_guide_chisq.qval <- as.data.frame(initial.guide.level.chisq.qval[["V1 mock"]])
mock_guide_chisq.qval$barcode <- row.names(mock_guide_chisq.qval)
colnames(mock_guide_chisq.qval) <- c("qval","barcode")

dox_target_chisq.qval <- as.data.frame(initial.target.level.chisq.qval[["V1 dox 100nm"]])
dox_target_chisq.qval$barcode <- row.names(dox_target_chisq.qval)
colnames(dox_target_chisq.qval) <- c("qval","target")

dox_guide_chisq.qval <- as.data.frame(initial.guide.level.chisq.qval[["V1 dox 100nm"]])
dox_guide_chisq.qval$barcode <- row.names(dox_guide_chisq.qval)
colnames(dox_guide_chisq.qval) <- c("qval","barcode")

mock_guide_plot_qval_df <- merge(gene_barcode_association, mock_guide_chisq.qval, by.x = "barcode", by.y = "barcode")

dox_guide_plot_qval_df <- merge(gene_barcode_association, dox_guide_chisq.qval, by.x = "barcode", by.y = "barcode")

mock_guide_plot_qval_df <- mock_guide_plot_qval_df %>% group_by(gene) %>% mutate(min=min(qval), average=mean(qval)) %>% arrange(desc(min), desc(average))

dox_guide_plot_qval_df= dox_guide_plot_qval_df %>% group_by(gene) %>% mutate(min=min(qval), average=mean(qval)) %>% arrange(desc(min), desc(average))

mock_guide_plot_qval_df$gene <- factor(mock_guide_plot_qval_df$gene, levels = unique(mock_guide_plot_qval_df$gene))

dox_guide_plot_qval_df$gene <- factor(dox_guide_plot_qval_df$gene, levels = unique(dox_guide_plot_qval_df$gene))

ggplot(mock_guide_plot_qval_df, aes(x = gene, y = -log10(qval), fill = -log10(qval) > 1.3)) + 
geom_point(pch=21, alpha=0.5) + 
theme_cfg() +
geom_hline(yintercept = 1.3, color = "grey80", linetype="dashed") +
xlab('target') +
ylab('-log10(qval) guide') +
coord_flip() +
scale_fill_manual(name = 'q < 0.05', values = setNames(c('red','black'),c(T, F))) +
theme(legend.position = "none") + 
ggsave("supplemental_figures/mock_guide_level_chisq.png", width = 2.25, height = 4.5)

ggplot(dox_guide_plot_qval_df, aes(x = gene, y = -log10(qval), fill = -log10(qval) > 1.3)) + 
geom_point(pch=21, alpha=0.5, position=position_jitter(width=0.25)) +
theme_cfg() +
geom_hline(yintercept = 1.30, color = "grey80", linetype="dashed") +
xlab('target') +
ylab('-log10(qval) guide') +
coord_flip() +
scale_fill_manual(name = 'q < 0.05', values = setNames(c('red','black'),c(T, F))) +
theme(legend.position = "none") + 
ggsave("supplemental_figures/dox_guide_level_chisq.png", width = 2.25, height = 4.5)

mock_target_chisq.qval <- mock_target_chisq.qval[order(mock_target_chisq.qval$qval, decreasing = TRUE),]

mock_target_chisq.qval$target <- factor(mock_target_chisq.qval$target, levels = mock_target_chisq.qval$target)

ggplot(mock_target_chisq.qval[order(mock_target_chisq.qval$qval, decreasing = TRUE),], aes(x = target, y = -log10(qval), fill = -log10(qval) > 1.3)) + 
geom_point(pch=21) + 
theme_cfg() +
xlab('target') +
ylab('-log10(qval) target') +
coord_flip() +
geom_hline(yintercept = 1.30, color = "grey80", linetype="dashed") +
scale_fill_manual(name = 'q < 0.05', values = setNames(c('red','black'),c(T, F))) +
theme(legend.position = "none") + 
ggsave("supplemental_figures/mock_target_level_chisq.png", width = 2.25, height = 4.5)

dox_target_chisq.qval <- dox_target_chisq.qval[order(dox_target_chisq.qval$qval, decreasing = TRUE),]

dox_target_chisq.qval$target <- factor(dox_target_chisq.qval$target, levels = dox_target_chisq.qval$target)

ggplot(dox_target_chisq.qval, aes(x = target, y = -log10(qval), fill = -log10(qval) > 1.3)) + 
geom_point(pch=21) +
xlab('target') +
ylab('-log10(qval) target') +
coord_flip() + 
geom_hline(yintercept = 1.30, color = "grey80", linetype="dashed") +
scale_fill_manual(name = 'q < 0.05', values = setNames(c('red','black'),c(T, F))) +
theme(legend.position = "none") + 
theme_cfg() +
ggsave("supplemental_figures/dox_target_level_chisq.png", width = 2.25, height = 4.5)
