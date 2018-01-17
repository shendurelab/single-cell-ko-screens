replace_gene_names_vec <- function(input_vec, name_vec, retain_inds = c(-1,-2)) {
  temp <- merge(name_vec, input_vec, by="row.names")
  temp2 <- temp[,retain_inds]
  names(temp2) <- temp[,2]
  return(temp2)
}

plot_gsa_hyper_heatmap <- function(gsa_results, significance=0.05)
{
  hyper_df <- ldply(gsa_results, function(gsa_res)
  {
    data.frame(gene_set = names(gsa_res$pvalues), pval = gsa_res$pvalues, qval = gsa_res$p.adj)
  })
  colnames(hyper_df)[1] <- "cluster_id"
  
  #hyper_df 
  
  hyper_df <- subset(hyper_df, qval <= significance)
  print (head(hyper_df))
  hyper_df <- merge(hyper_df, ddply(hyper_df, .(gene_set), function(x) { nrow(x) }), by="gene_set")
  #print (hyper_df)
  hyper_df$gene_set <- factor(hyper_df$gene_set, levels=unique(arrange(hyper_df, V1, cluster_id)$gene_set))
  
  qplot(cluster_id, gene_set, fill=-log10(qval), geom="tile", data=hyper_df) + scale_fill_gradientn(colours=rainbow(7)) + theme(text = element_text(size=4), axis.text.x = element_text(size=6), axis.title.x = element_text(size=6), axis.title.y = element_text(size=6), legend.title=element_text(size=6),legend.text=element_text(size=6))
}

plot_gsa_hyper_heatmap_cfg <- function(gsa_results, significance=0.05)
{
  hyper_df <- ldply(gsa_results, function(gsa_res)
  {
    data.frame(gene_set = names(gsa_res$pvalues), pval = gsa_res$pvalues, qval = gsa_res$p.adj)
  })
  colnames(hyper_df)[1] <- "Gene_ko"
  
  #hyper_df 
  
  hyper_df <- subset(hyper_df, pval <= significance)
  print (head(hyper_df))
  hyper_df <- merge(hyper_df, ddply(hyper_df, .(gene_set), function(x) { nrow(x) }), by="gene_set")
  #print (hyper_df)
  hyper_df$gene_set <- factor(hyper_df$gene_set, levels=unique(arrange(hyper_df, V1, Gene_ko)$gene_set))
  
  qplot(cluster_id, gene_set, fill=-log10(pval), geom="tile", data=hyper_df) + scale_fill_gradientn(colours=rainbow(7)) + theme(text = element_text(size=4), axis.text.x = element_text(size=10, angle = 90, hjust = 1), axis.title.x = element_text(size=6), axis.title.y = element_text(size=6), legend.title=element_text(size=6),legend.text=element_text(size=6))
}


collect_gsa_hyper_results_clusters <- function(genes_list, clusters, gsc)
{
  gene_universe <- unique(as.character(genes_list))
  gsa_results <- list()
  cluster_ids <- unique(clusters)
  for (i in (1:length(cluster_ids))) {
    cluster_genes <- unique(names(clusters[clusters == i]))
    gsaRes <- runGSAhyper(cluster_genes, gsc=gsc, universe=gene_universe, adjMethod = "BH")
    gsa_results[[length(gsa_results) + 1]] <- gsaRes
  }
  names(gsa_results) <- cluster_ids
  gsa_results
}

collect_gsa_hyper_results_All <- function(genes, universe, gsc)
{
  gene_universe <- unique(as.character(universe))
  gsa_results <- list()
  cluster_ids=c('diffExp')
  gsaRes <- runGSAhyper(as.character(genes), gsc=gsc, universe=gene_universe, adjMethod ="BH")
  gsa_results[[length(gsa_results) + 1]] <- gsaRes
  names(gsa_results) <- cluster_ids
  gsa_results
}

