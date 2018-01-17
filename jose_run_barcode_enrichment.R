suppressPackageStartupMessages({
    library(RColorBrewer)
    library(plyr)
    library(dplyr)
    library(tidyr)
    library(piano)
    library(ggplot2)
    library(glmnet)
    library(monocle)
    library(parallel)
    
    library(reshape2)
    library(scales)
    library(viridis)
    library(argparse)
})
source('helper_functions.R')

simulate_swap = function(aggregated_cds, swap_rate, seed=0) {
    set.seed(seed)
    pData(aggregated_cds)$id = row.names(pData(aggregated_cds))

    cells_to_keep = pData(aggregated_cds) %>%
        group_by(gene) %>%
        sample_frac(1 - swap_rate) %>%
        ungroup()

    # Swap the barcode and gene sequences (retaining their order)
    cells_to_shuffle = colnames(aggregated_cds)[! colnames(aggregated_cds) %in% cells_to_keep$id]
    shuffled_cells = sample(cells_to_shuffle)
    pData(aggregated_cds)[cells_to_shuffle, "gene"] = pData(aggregated_cds)[shuffled_cells, "gene"]
    pData(aggregated_cds)[cells_to_shuffle, "barcode"] = pData(aggregated_cds)[shuffled_cells, "barcode"]
    return(aggregated_cds)
}

parser = argparse::ArgumentParser(description='Script to perform barcode enrichment.')
parser$add_argument('--cds', default='temp_data/aggregated_cds.rds', help='RDS file with CDS object to use.')
parser$add_argument('--target_level_chisq', default='temp_data/barcode_enrichment/initial.target.level.chisq.qval.rds', help='RDS file with target level CHISQ results.')
parser$add_argument('--guide_level_chisq', default='temp_data/barcode_enrichment/initial.guide.level.chisq.qval.rds', help='RDS file with guide level CHISQ results.')
parser$add_argument('--processed_cds', default="temp_data/barcode_enrichment/aggregated_cds.barcode_enrichment.rds", help='RDS file with processed cds after barcode enrichment.')
parser$add_argument('--cds_mock', default="temp_data/barcode_enrichment/mock_cds.rds", help='Processed CDS from barcode enrichment for mock condition.')
parser$add_argument('--cds_dox', default= "temp_data/barcode_enrichment/dox_100nm_cds.rds", help='Processed CDS from barcode enrichment for dox condition.')
parser$add_argument('--processed_cds_pdata', default="temp_data/barcode_enrichment/aggregated_cds.barcode_enrichment.pdata.txt", help='pData table for processed_cds.')
parser$add_argument('--swap_rate', type='double', help='Simulate a swap rate prior to performing barcode enrichment.')
args = parser$parse_args()

cds.v1 = readRDS(args$cds)

cds.v1 = cds.v1[, with(pData(cds.v1), !is.na(proportion) & proportion >= 0.8 & guide_count == 1 & condition %in% c("mock_singles","dox_100nm_singles"))]

cds.v1.mock = cds.v1[, pData(cds.v1)$treatment == "mock"]
cds.v1.dox.100nm = cds.v1[, pData(cds.v1)$treatment == "dox_100nm"]

# Simulate swap rate
if (! is.null(args$swap_rate)) {
    cds.v1.mock = simulate_swap(cds.v1.mock, swap_rate=args$swap_rate, seed=0)
    cds.v1.dox.100nm = simulate_swap(cds.v1.dox.100nm, swap_rate=args$swap_rate, seed=0)
}

set.up.cds = function(cds) {
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds)
    cds = detectGenes(cds, 0.5)
    return(cds)
}

cds.v1.mock = set.up.cds(cds.v1.mock)
cds.v1.dox.100nm = set.up.cds(cds.v1.dox.100nm)

expressed_genes.v1.mock <-  row.names(fData(cds.v1.mock)[rowSums(exprs(cds.v1.mock) > 0) > 50 ,])
expressed_genes.v1.dox.100nm <-  row.names(fData(cds.v1.dox.100nm)[rowSums(exprs(cds.v1.dox.100nm) > 0) > 50 ,])

cds.v1.mock = reduceDimension(cds.v1.mock[expressed_genes.v1.mock], reduction_method = "tSNE",
    max_components = 3, norm_method = "log", num_dim = 20, verbose = T)

cds.v1.dox.100nm = reduceDimension(cds.v1.dox.100nm[expressed_genes.v1.dox.100nm], reduction_method = "tSNE",
    max_components = 3, norm_method = "log", num_dim = 20, verbose = T)

cds.v1.mock = clusterCells(cds.v1.mock)
cds.v1.dox.100nm = clusterCells(cds.v1.dox.100nm)

cds.v1.mock = clusterCells(cds.v1.mock,
    rho_threshold = 75, delta_threshold = 20, skip_rho_sigma = T)

plot_cell_clusters(cds.v1.mock, x = 1, y = 2, cell_size = 0.5, color_by = "Cluster") +
	theme_cfg() +
	ggsave('diagnostic_plots/mock.jose_enrichment.tsne.1_2.png')

plot_cell_clusters(cds.v1.mock, x = 1, y = 3, cell_size = 0.5, color_by = "Cluster") +
	theme_cfg() +
	ggsave('diagnostic_plots/mock.jose_enrichment.tsne.1_3.png')

plot_cell_clusters(cds.v1.mock, x = 2, y = 3, cell_size = 0.5, color_by = "Cluster") +
	theme_cfg() +
	ggsave('diagnostic_plots/mock.jose_enrichment.tsne.2_3.png')


cds.v1.dox.100nm = clusterCells(cds.v1.dox.100nm,
    rho_threshold = 75, delta_threshold = 20, skip_rho_sigma = T)

plot_cell_clusters(cds.v1.dox.100nm, x = 1, y = 2, cell_size = 0.5, color_by = "Cluster") +
        theme_cfg() +
        ggsave('diagnostic_plots/dox.jose_enrichment.tsne.1_2.png')

plot_cell_clusters(cds.v1.dox.100nm, x = 1, y = 3, cell_size = 0.5, color_by = "Cluster") +
        theme_cfg() +
        ggsave('diagnostic_plots/dox.jose_enrichment.tsne.1_3.png')

plot_cell_clusters(cds.v1.dox.100nm, x = 2, y =3, cell_size = 0.5, color_by = "Cluster") +
        theme_cfg() +
        ggsave('diagnostic_plots/dox.jose_enrichment.tsne.2_3.png')

# Guide level analysis
analysis.guides = list()

analysis.guides[["V1 mock"]] =
    (pData(cds.v1.mock) %>% filter(guide_count == 1, gene != "NONTARGETING") %>% group_by(gene, barcode) %>%
    summarize(n.guide.cells = n()) %>% group_by(gene) %>% mutate(n.target.cells = sum(n.guide.cells)) %>%
    filter(n.guide.cells >= 10) %>% ungroup())$barcode

analysis.guides[["V1 dox 100nm"]] =
    (pData(cds.v1.dox.100nm) %>% filter(guide_count == 1, gene != "NONTARGETING") %>% group_by(gene, barcode) %>%
    summarize(n.guide.cells = n()) %>% group_by(gene) %>% mutate(n.target.cells = sum(n.guide.cells)) %>%
    filter(n.guide.cells >= 10) %>% ungroup())$barcode

# Target level analysis
analysis.targets = list()

analysis.targets[["V1 mock"]] = as.data.frame(pData(cds.v1.mock) %>%
    group_by(gene) %>% summarize(
        n.cells = n(),
        n.guides = length(intersect(unique(barcode), analysis.guides[["V1 mock"]]))) %>%
    filter(n.cells >= 15, n.guides >= 1) %>% select(gene))[, 1]

analysis.targets[["V1 dox 100nm"]] = as.data.frame(pData(cds.v1.dox.100nm) %>%
    group_by(gene) %>% summarize(
        n.cells = n(),
        n.guides = length(intersect(unique(barcode), analysis.guides[["V1 dox 100nm"]]))) %>%
    filter(n.cells >= 15, n.guides >= 1) %>% select(gene))[, 1]

target.to.guide.map = list()

target.to.guide.map[["V1 mock"]] = list()
target.to.guide.map[["V1 dox 100nm"]] = list()

for (target in analysis.targets[["V1 mock"]]) {
    target.to.guide.map[["V1 mock"]][[target]] =
        sort(unique(as.data.frame(pData(cds.v1.mock) %>%
            filter(gene == target, barcode %in% analysis.guides[["V1 mock"]]) %>%
            select(barcode))[, 1]))
}


for (target in analysis.targets[["V1 dox 100nm"]]) {
    target.to.guide.map[["V1 dox 100nm"]][[target]] =
        sort(unique(as.data.frame(pData(cds.v1.dox.100nm) %>%
            filter(gene == target, barcode %in% analysis.guides[["V1 dox 100nm"]]) %>%
            select(barcode))[, 1]))
}


guide.to.target.map = list()

guide.to.target.map[["V1 mock"]] = list()
guide.to.target.map[["V1 dox 100nm"]] = list()

for (target in analysis.targets[["V1 mock"]]) {
    for (guide in target.to.guide.map[["V1 mock"]][[target]]) {
        guide.to.target.map[["V1 mock"]][[guide]] = target
    }
}

for (target in analysis.targets[["V1 dox 100nm"]]) {
    for (guide in target.to.guide.map[["V1 dox 100nm"]][[target]]) {
        guide.to.target.map[["V1 dox 100nm"]][[guide]] = target
    }
}


target.cluster.mat = list()

target.cluster.mat[["V1 mock"]] = acast(
    pData(cds.v1.mock) %>%
    filter(barcode %in% analysis.guides[["V1 mock"]] | gene == "NONTARGETING") %>%
    mutate(dummy = 1) %>% select(gene, Cluster, dummy),
    gene ~ Cluster, value.var = "dummy", fun.aggregate = sum, fill = 0)

target.cluster.mat[["V1 dox 100nm"]] = acast(
    pData(cds.v1.dox.100nm) %>%
    filter(barcode %in% analysis.guides[["V1 dox 100nm"]] | gene == "NONTARGETING") %>%
    mutate(dummy = 1) %>% select(gene, Cluster, dummy),
    gene ~ Cluster, value.var = "dummy", fun.aggregate = sum, fill = 0)


NTC.cluster.p = list()

NTC.cluster.p[["V1 mock"]] <- pData(cds.v1.mock)[pData(cds.v1.mock)$gene == "NONTARGETING",] %>%
    group_by(Cluster) %>%
    summarize(n = n()) %>%
    complete(Cluster, fill = list(n = 0.1))


NTC.cluster.p[["V1 dox 100nm"]] <- pData(cds.v1.dox.100nm)[pData(cds.v1.dox.100nm)$gene == "NONTARGETING",] %>%
    group_by(Cluster) %>%
    summarize(n = n()) %>%
    complete(Cluster, fill = list(n = 0.1))

guide.cluster.mat = list()

guide.cluster.mat[["V1 mock"]] = acast(
    pData(cds.v1.mock) %>% filter(barcode %in% analysis.guides[["V1 mock"]]) %>%
    mutate(dummy = 1) %>% select(barcode, Cluster, dummy),
    barcode ~ Cluster, value.var = "dummy", fun.aggregate = sum, fill = 0)

guide.cluster.mat[["V1 dox 100nm"]] = acast(
    pData(cds.v1.dox.100nm) %>% filter(barcode %in% analysis.guides[["V1 dox 100nm"]]) %>%
    mutate(dummy = 1) %>% select(barcode, Cluster, dummy),
    barcode ~ Cluster, value.var = "dummy", fun.aggregate = sum, fill = 0)

ntc.distribution = list()

ntc.distribution[["V1 mock"]] = target.cluster.mat[["V1 mock"]]["NONTARGETING",]
ntc.distribution[["V1 dox 100nm"]] = target.cluster.mat[["V1 dox 100nm"]]["NONTARGETING",]

ntc.distribution[["V1 mock"]] = ntc.distribution[["V1 mock"]] / sum(ntc.distribution[["V1 mock"]])
ntc.distribution[["V1 dox 100nm"]] = ntc.distribution[["V1 dox 100nm"]] / sum(ntc.distribution[["V1 dox 100nm"]])

initial.target.level.chisq.pval = list()

set.seed(42)
initial.target.level.chisq.pval[["V1 mock"]] = sapply(
    analysis.targets[["V1 mock"]], function(target) {

    message(target)
    chisq.test(
        target.cluster.mat[["V1 mock"]][target,],
        p = NTC.cluster.p[["V1 mock"]]$n,
        simulate.p.value = F, rescale.p = T, B = 20000)$p.value
})

set.seed(42)
initial.target.level.chisq.pval[["V1 dox 100nm"]] = sapply(
    analysis.targets[["V1 dox 100nm"]], function(target) {

    message(target)
    chisq.test(
        target.cluster.mat[["V1 dox 100nm"]][target,],
        p = NTC.cluster.p[["V1 dox 100nm"]]$n,
        simulate.p.value = F, rescale.p = T, B = 100000)$p.value
})


initial.guide.level.chisq.pval = list()

set.seed(42)
initial.guide.level.chisq.pval[["V1 mock"]] = sapply(
    analysis.guides[["V1 mock"]], function(guide) {

    message(guide)
    chisq.test(
        guide.cluster.mat[["V1 mock"]][guide,],
        p = NTC.cluster.p[["V1 mock"]]$n,
        simulate.p.value = F, rescale.p = T, B = 20000)$p.value
})

set.seed(42)
initial.guide.level.chisq.pval[["V1 dox 100nm"]] = sapply(
    analysis.guides[["V1 dox 100nm"]], function(guide) {

    message(guide)
    chisq.test(
        guide.cluster.mat[["V1 dox 100nm"]][guide,],
        p = NTC.cluster.p[["V1 dox 100nm"]]$n,
        simulate.p.value = F, rescale.p = T, B = 20000)$p.value
})

initial.target.level.chisq.qval <- list()

initial.target.level.chisq.qval[["V1 mock"]] <- p.adjust(initial.target.level.chisq.pval[["V1 mock"]], method = "BH")
initial.target.level.chisq.qval[["V1 dox 100nm"]] <- p.adjust(initial.target.level.chisq.pval[["V1 dox 100nm"]], method = "BH")

initial.guide.level.chisq.qval <- list()

initial.guide.level.chisq.qval[["V1 mock"]] <- p.adjust(initial.guide.level.chisq.pval[["V1 mock"]], method = "BH")
initial.guide.level.chisq.qval[["V1 dox 100nm"]] <- p.adjust(initial.guide.level.chisq.pval[["V1 dox 100nm"]], method = "BH")

initial.target.level.chisq.qval[["V1 mock"]] <- sapply(initial.target.level.chisq.qval[["V1 mock"]],
                                                       function(x){if(x < 1e-50){return(1e-50)}else{return(x)}})

initial.target.level.chisq.qval[["V1 dox 100nm"]] <- sapply(initial.target.level.chisq.qval[["V1 dox 100nm"]],
                                                       function(x){if(x < 1e-50){return(1e-50)}else{return(x)}})

initial.guide.level.chisq.qval[["V1 mock"]] <- sapply(initial.guide.level.chisq.qval[["V1 mock"]],
                                                       function(x){if(x < 1e-50){return(1e-50)}else{return(x)}})

initial.guide.level.chisq.qval[["V1 dox 100nm"]] <- sapply(initial.guide.level.chisq.qval[["V1 dox 100nm"]],
                                                       function(x){if(x < 1e-50){return(1e-50)}else{return(x)}})

saveRDS(initial.target.level.chisq.qval, args$target_level_chisq)
saveRDS(initial.guide.level.chisq.qval, args$guide_level_chisq)

pass.target.level.screen = list()

pass.target.level.screen[["V1 mock"]] =
    sort(names(which(initial.target.level.chisq.qval[["V1 mock"]] < 0.05)))

pass.target.level.screen[["V1 dox 100nm"]] =
    sort(names(which(initial.target.level.chisq.qval[["V1 dox 100nm"]] < 0.05)))


pass.guide.level.screen = list()

pass.guide.level.screen[["V1 mock"]] = sort(unlist(unique(sapply(
    names(which(initial.guide.level.chisq.qval[["V1 mock"]] < 0.05)), function(guide) {
        guide.to.target.map[["V1 mock"]][[guide]]
    }))))

pass.guide.level.screen[["V1 dox 100nm"]] = sort(unlist(unique(sapply(
    names(which(initial.guide.level.chisq.qval[["V1 dox 100nm"]] < 0.05)), function(guide) {
        guide.to.target.map[["V1 dox 100nm"]][[guide]]
    }))))


targets.passing.initial.screen = list()

targets.passing.initial.screen[["V1 mock"]] = sort(union(
    pass.target.level.screen[["V1 mock"]], pass.guide.level.screen[["V1 mock"]]))

targets.passing.initial.screen[["V1 dox 100nm"]] = sort(union(
    pass.target.level.screen[["V1 dox 100nm"]], pass.guide.level.screen[["V1 dox 100nm"]]))


#EM
get.guide.weights = function(mat, ntc.dist, n.iterations = 30) {
    n.guides = nrow(mat)
    n.cells = rowSums(mat)
    empirical.dist = sweep(mat, 1, n.cells, "/")

    lof.prop = rep(0.5, n.guides)
    expected.n.lof = n.cells * lof.prop

    for (i in 1:n.iterations) {
        lof.dist = sapply(1:n.guides, function(guide) {
            p = lof.prop[guide]
            (empirical.dist[guide,] - (1-p) * ntc.dist) / p
        })

        lof.dist = rowSums(sweep(lof.dist, 2, expected.n.lof / sum(expected.n.lof), "*"))
        lof.dist = ifelse(lof.dist < 0, 0, lof.dist)
        lof.dist = lof.dist / sum(lof.dist)

        lof.prop = sapply(1:n.guides, function(guide) {
            optimize(function(p) dmultinom(mat[guide,], prob = p * lof.dist + (1-p) * ntc.dist, log = T),
                c(0.0, 1.0), maximum = T)$maximum
        })

        expected.n.lof = n.cells * lof.prop
    }

    return(lof.prop)
}

weighted.target.cluster.mat = list()

for (condition in c("V1 mock", "V1 dox 100nm")) {
    weighted.target.cluster.mat[[condition]] = t(sapply(targets.passing.initial.screen[[condition]],
            function(target) {
        guides = target.to.guide.map[[condition]][[target]]
        if (length(guides) == 1) {
            return(target.cluster.mat[[condition]][target,])
        } else {
            mat = guide.cluster.mat[[condition]][guides,]
            guide.weights = get.guide.weights(mat, ntc.distribution[[condition]])
            guide.weights = guide.weights / max(guide.weights)

            #return(target.cluster.mat[[condition]][target,])
            return(round(colSums(sweep(mat, 1, guide.weights, "*"))))
        }
    }))
}


cluster.enrichment.df = list()

for (condition in c("V1 mock", "V1 dox 100nm")) {
    weighted.mat = weighted.target.cluster.mat[[condition]]
    ntc.counts = target.cluster.mat[[condition]]["NONTARGETING",]

    cluster.enrichment.df[[condition]] = do.call(rbind, lapply(rownames(weighted.mat), function(target) {
        do.call(rbind, lapply(1:ncol(weighted.mat), function(cluster) {
            test = fisher.test(cbind(
                c(weighted.mat[target, cluster], sum(weighted.mat[target, -cluster])),
                c(ntc.counts[cluster], sum(ntc.counts[-cluster]))))

            data.frame(
                target = target,
                cluster = cluster,
                odds.ratio = unname(test$estimate),
                p.value = test$p.value)
        }))
    }))

    cluster.enrichment.df[[condition]]$q.value = p.adjust(cluster.enrichment.df[[condition]]$p.value, "fdr")

    cluster.enrichment.df[[condition]]$log2.odds = with(cluster.enrichment.df[[condition]],
        ifelse(odds.ratio == 0, -5, log2(odds.ratio)))

    cluster.enrichment.df[[condition]]$heatmap.value = with(cluster.enrichment.df[[condition]],
        ifelse(q.value < 0.1, log2.odds, 0))
}

cluster.enrichment.df[["V1 mock"]]$id <- paste0(cluster.enrichment.df[["V1 mock"]]$target,"_",cluster.enrichment.df[["V1 mock"]]$cluster)

cluster.enrichment.df[["V1 dox 100nm"]]$id <- paste0(cluster.enrichment.df[["V1 dox 100nm"]]$target,"_",cluster.enrichment.df[["V1 dox 100nm"]]$cluster)


################### Return cells flagged as informative by Fishers exact test ######################

### have been using p.value < 0.05 and odds_ratio > 1 as cutoff

flag_functional_edits <- function(cds, enrichment_output){


    sig_enriched_subset <- enrichment_output %>% filter((q.value < 0.1 & log2.odds > 0) | (q.value < 0.1 & log2.odds == Inf))

    functional_edits <- list()

    enriched_targets_byCluster <- sig_enriched_subset$id


    for (enriched in enriched_targets_byCluster){
         print(enriched)
        edited_Cells_byCluster <- row.names(pData(cds)[pData(cds)$Cluster == sig_enriched_subset[sig_enriched_subset$id == enriched,]$cluster &
                                                       pData(cds)$gene == sig_enriched_subset[sig_enriched_subset$id == enriched,]$target,])

        functional_edits[[which(enriched_targets_byCluster == enriched)]] <- edited_Cells_byCluster

    }


        NTC_id <- row.names(pData(cds)[pData(cds)$gene == "NONTARGETING",])

        functional_edits_id <- union(unlist(functional_edits), NTC_id)



    return(functional_edits_id)

}

mock.informative <- flag_functional_edits(cds = cds.v1.mock, enrichment_output = cluster.enrichment.df[["V1 mock"]])

dox.100nm.informative <- flag_functional_edits(cds.v1.dox.100nm, cluster.enrichment.df[["V1 dox 100nm"]])

pData(cds.v1.mock)$informative <- row.names(pData(cds.v1.mock)) %in% mock.informative

pData(cds.v1.dox.100nm)$informative <- row.names(pData(cds.v1.dox.100nm)) %in% dox.100nm.informative

pData(cds.v1.mock)$used_in_enrichment <- (rep(TRUE, length(row.names(pData(cds.v1.mock)))))
pData(cds.v1.dox.100nm)$used_in_enrichment <- (rep(TRUE, length(row.names(pData(cds.v1.dox.100nm)))))

mock_used_in_enrichment <- row.names(pData(cds.v1.mock))
dox_used_in_enrichment <- row.names(pData(cds.v1.dox.100nm))

pData(cds.v1)$informative <- row.names(pData(cds.v1)) %in% c(mock.informative, dox.100nm.informative)

pData(cds.v1)$used_in_enrichment <- row.names(pData(cds.v1)) %in% c(mock_used_in_enrichment, dox_used_in_enrichment)

saveRDS(cds.v1, args$processed_cds)

write.table(pData(cds.v1), sep='\t', quote=F, row.names=F, file=args$processed_cds_pdata)

saveRDS(cds.v1.mock, args$cds_mock)
saveRDS(cds.v1.dox.100nm, args$cds_dox)
