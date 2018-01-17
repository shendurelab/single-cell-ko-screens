library(monocle)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(matrixStats)
library(viridis)
library(tidyr)
library(stringr)
library(glmnet)
library(parallel)

library(reshape2)
library(scales)

.tenx_to_cds = function(pipeline_dirs, genome="hg19", filtered=TRUE) {
  # Takes a list of 10X pipeline output directories and generates a cellDataSet containing all cells in these experiments
  #
  # Args:
  #pipeline_dirs: Directory name or list of directory names of the top level 10X output directory for an experiment(s)
  # genome: String with genome name specified for 10X run (such as hg19)
  # filtered: bool indicating whether you want to load the filtered matrix or not (unfiltered)
  #
  # Returns:
  # A cellDataSet object containing data from all experiments.

  # Outer scope variables for collecting data
  expression_matrices = list()
  metadata_dfs = list()
  gene_table = NULL # will keep first gene table so can match ordering for each dataset

  lapply(pipeline_dirs, function(pipeline_dir) {
    # Check initial user input
    if( ! file.exists(pipeline_dir) ) { stop(paste("Specified 10X output directory does not exist:", pipeline_dir)) }

    # Construct paths for pipeline output files and check that they exist
    if (filtered) {
      base_path = file.path(pipeline_dir, "outs", "filtered_gene_bc_matrices_mex", genome)
      base_path.unaggregated = file.path(pipeline_dir, "outs", "filtered_gene_bc_matrices", genome)
      } else {
        base_path = file.path(pipeline_dir, "outs", "raw_gene_bc_matrices_mex", genome)
        base_path.unaggregated = file.path(pipeline_dir, "outs", "raw_gene_bc_matrices", genome)
      }

    # Get the file path
    if(file.exists(base_path)) {
      base_path = base_path
      } else if(file.exists(base_path.unaggregated)) {
      # This is an aggregated run, change the base path
      base_path = base_path.unaggregated
      } else {
      # No expected directories were found
      stop(paste("Specified genome does not appear in 10X output:", base_path, ' or ', base_path.unaggregated))
    }

    matrix_path = file.path(base_path, "matrix.mtx")
    genes_path = file.path(base_path, "genes.tsv")
    barcodes_path = file.path(base_path, "barcodes.tsv")
    analysis_path = file.path(pipeline_dir, "outs", "analysis")

    if( ! file.exists(matrix_path) ) { stop(paste("Expression matrix not found in 10X output:", matrix_path)) }
    if( ! file.exists(genes_path) ) { stop(paste("Genes file not found in 10X output:", genes_path)) }
    if( ! file.exists(barcodes_path) ) { stop(paste("Barcodes file not found in 10X output:", barcodes_path)) }
    if( ! file.exists(analysis_path) ) { stop(paste("Analysis  path not found in 10X output:", analysis_path)) }

    # All files exist, read them in
    matrix = Matrix::readMM(matrix_path)
    barcodes = read.table(barcodes_path, header=F, as.is=T)[,1]
    current_gene_table = read.table(genes_path, header=F, as.is=T) ## saves for later

    if ( is.null(gene_table) ) { gene_table <<- current_gene_table } ## store the first gene table so can match ordering between all experiments
    genes = current_gene_table[, 1]

    # Add gene and sample names to expression matrix (adding dataset post-fix in case barcodes appear in multiple samples)
    sample = basename(pipeline_dir)
    row.names(matrix) = genes
    colnames(matrix) = paste(barcodes, sample, sep="_") ## adds dataset post-fix
    matrix = matrix[gene_table[, 1], ] ## ensures order of genes matches between experiments

    # Construct metadata table that includes directory samples come from and other stats
    total_umis = colSums(matrix)

    metadata_df = data.frame(
      cell = colnames(matrix),
      total_umis = total_umis,
      sample = sample
      )

    # Add both matrices to the running list
    expression_matrices[[length(expression_matrices) + 1]] <<- matrix
    metadata_dfs[[length(metadata_dfs) + 1]] <<- metadata_df
    })

  # Now combine all the dataframes into one and make CDS
  combined_expression_matrix = do.call(cBind, expression_matrices)
  row.names(combined_expression_matrix) = gene_table[, 1]

  combined_metadata_df = do.call(rbind, metadata_dfs)
  row.names(combined_metadata_df) = combined_metadata_df$cell

  colnames(gene_table) = c("id", "gene_short_name")
  row.names(gene_table) = gene_table$id

  pd = new("AnnotatedDataFrame", data = combined_metadata_df)
  fd = new("AnnotatedDataFrame", data = gene_table)
  cds = newCellDataSet(combined_expression_matrix,
   phenoData=pd,
   featureData=fd,
   expressionFamily=negbinomial.size(),
   lowerDetectionLimit=0.5)
  return(cds)
}

### Not in function
'%!in%' <- function(x,y)!('%in%'(x,y))

# Helper function for GSEA stuff from Jose
replace_gene_names_vec <- function(input_vec, name_vec, retain_inds = c(-1,-2)) {
  temp <- merge(name_vec, input_vec, by="row.names")
  temp2 <- temp[,retain_inds]
  names(temp2) <- temp[,2]
  return(temp2)
}

#EM approach for getting the guide weights to go into final enrichment testing
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

# Wrapper for getting information about enrichment in clusters
flag_enriched_clusters = function(cds, cell_detection_threshold, rho_threshold, delta_threshold, rho_delta_plot, pca_plot, tsne_plot, max_components=12, guide_count_min=10, target_count_min=15, qval_threshold=0.1) {
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds)
    cds = detectGenes(cds, 0.5)

    expressed_genes = row.names(fData(cds)[rowSums(exprs(cds) > 0) > cell_detection_threshold ,])

    # PCA variance explained plot
    plot_pc_variance_explained(cds[expressed_genes], norm_method = "log", max_components=50, pseudo_expr = 1, return_all = F) +
      geom_vline(xintercept=max_components, color="red", linetype="dashed", size=0.8) +
      theme(axis.text = element_text(face="bold"), axis.title = element_text(face="bold")) +
      ylab('variance explained by component') +
      xlab('component') +
      ggsave(pca_plot, height=3, width=3)

    cds = reduceDimension(cds[expressed_genes], reduction_method = "tSNE",
    max_components = 3, norm_method = "log", num_dim = max_components, verbose = T)

    # Cluster cells (required for rho delta)
    cds = clusterCells(cds,
      verbose = T,
      method = "densityPeak")

    # Produce rho delta plot prior to clustering
    plot_rho_delta(cds, rho_threshold = 1000, delta_threshold = 9) +
      geom_point(aes(colour=rho > rho_threshold & delta > delta_threshold)) +
      geom_vline(xintercept=rho_threshold, color="red", linetype="dashed", size=0.8) +
      geom_hline(yintercept=delta_threshold, color="red", linetype="dashed", size=0.8) +
      theme(axis.text = element_text(face="bold"), axis.title = element_text(face="bold")) +
      ylab('density peak delta') +
      xlab('density peak rho') +
      scale_colour_manual(values=c("TRUE"="red", "FALSE"="#d3d3d3")) +
      guides(color=FALSE) +
      ggsave(rho_delta_plot, height=3, width=3)

    cds = clusterCells(cds, rho_threshold=rho_threshold, delta_threshold=delta_threshold, skip_rho_sigma=TRUE)

    # Plot TSNE in 3D
    tsne_plot1 = plot_cell_clusters(cds, x = 1, y = 2, cell_size = 0.5) + theme_cfg()
    tsne_plot2 = plot_cell_clusters(cds, x = 1, y = 3, cell_size = 0.5) + theme_cfg()
    tsne_plot3 = plot_cell_clusters(cds, x = 2, y = 3, cell_size = 0.5) + theme_cfg()
    ggsave(tsne_plot, arrangeGrob(tsne_plot1, tsne_plot2, tsne_plot3, ncol=3), height=3.5, width=9)

    analysis.guides =
        (pData(cds) %>%
            dplyr::filter(gene != "NONTARGETING") %>%
            dplyr::group_by(gene, barcode) %>%
                dplyr::summarize(n.guide.cells = n()) %>%
            dplyr::group_by(gene) %>%
                dplyr::mutate(n.target.cells = sum(n.guide.cells)) %>%
                dplyr::filter(n.guide.cells >= guide_count_min) %>%
                dplyr::ungroup())$barcode


    analysis.targets = as.data.frame(pData(cds) %>%
                            dplyr::group_by(gene) %>% dplyr::summarize(n.cells = n(),
                                                         n.guides = length(intersect(unique(barcode), analysis.guides))) %>%
                            dplyr::filter(n.cells >= target_count_min, n.guides >= 1) %>%
                            dplyr::select(gene))[, 1]

    target.to.guide.map = list()

    for (target in analysis.targets) {
        target.to.guide.map[[target]] =
            sort(unique(as.data.frame(pData(cds) %>%
                dplyr::filter(gene == target & barcode %in% analysis.guides) %>%
                dplyr::select(barcode))[, 1]))
    }


    guide.to.target.map = list()

    for (target in analysis.targets) {
        for (guide in target.to.guide.map[[target]]) {
            guide.to.target.map[[guide]] = target
        }
    }


    target.cluster.mat = acast(
        pData(cds) %>%
        dplyr::filter(barcode %in% analysis.guides | gene == "NONTARGETING") %>%
        dplyr::mutate(dummy = 1) %>%
        dplyr::select(gene, Cluster, dummy),
               gene ~ Cluster,
               value.var = "dummy",
               fun.aggregate = sum,
               fill = 0)

    NTC.cluster.p <- pData(cds)[pData(cds)$gene == "NONTARGETING",] %>%
        dplyr::group_by(Cluster) %>%
        dplyr::summarize(n = n()) %>%
        tidyr::complete(Cluster, fill = list(n = 0.1))


    guide.cluster.mat = acast(
        pData(cds) %>%
        dplyr::filter(barcode %in% analysis.guides) %>%
        dplyr::mutate(dummy = 1) %>%
        dplyr::select(barcode, Cluster, dummy),
            barcode ~ Cluster,
            value.var = "dummy",
            fun.aggregate = sum,
            fill = 0)

    ntc.distribution = target.cluster.mat["NONTARGETING",] / sum(target.cluster.mat["NONTARGETING",])

    set.seed(42)
    initial.target.level.chisq.pval = sapply(
        analysis.targets, function(target) {

        message(target)
        chisq.test(
            target.cluster.mat[target,],
            p = NTC.cluster.p$n,
            simulate.p.value = T,
            rescale.p = T,
            B = 20000)$p.value
    })


    set.seed(42)
    initial.guide.level.chisq.pval = sapply(
        analysis.guides, function(guide) {

        message(guide)
        chisq.test(
            guide.cluster.mat[guide,],
            p = NTC.cluster.p$n,
            simulate.p.value = T,
            rescale.p = T,
            B = 20000)$p.value
    })

    pass.target.level.screen =
        sort(names(which(initial.target.level.chisq.pval < 0.05 /
            length(initial.target.level.chisq.pval))))


    pass.guide.level.screen = sort(unlist(unique(sapply(
        names(which(initial.guide.level.chisq.pval < 0.05 /
            length(initial.guide.level.chisq.pval))), function(guide) {
            guide.to.target.map[[guide]]
        }))))

    targets.passing.initial.screen = sort(union(
        pass.target.level.screen, pass.guide.level.screen))

    # Get the guide weights
    all_guide_weights = data.frame(guide=c(), weight=c())

    weighted.target.cluster.mat = t(sapply(targets.passing.initial.screen,
            function(target) {
                guides = target.to.guide.map[[target]]
                if (length(guides) == 1) {
                    return(target.cluster.mat[target,])
                } else {
                    mat = guide.cluster.mat[guides,]
                    guide.weights = get.guide.weights(mat, ntc.distribution)
                    guide.weights = guide.weights / max(guide.weights)
                    all_guide_weights <<- rbind(all_guide_weights, data.frame(guide=guides, weight=guide.weights))
                    return(round(colSums(sweep(mat, 1, guide.weights, "*"))))
                }
            }))

    ntc.counts = target.cluster.mat["NONTARGETING",]

    cluster.enrichment.df = do.call(rbind, lapply(rownames(weighted.target.cluster.mat), function(target) {
        do.call(rbind, lapply(1:ncol(weighted.target.cluster.mat), function(cluster) {
            test = fisher.test(cbind(
                c(weighted.target.cluster.mat[target, cluster], sum(weighted.target.cluster.mat[target, -cluster])),
                c(ntc.counts[cluster], sum(ntc.counts[-cluster]))))

            data.frame(
                target = target,
                cluster = cluster,
                odds.ratio = unname(test$estimate),
                p.value = test$p.value)
        }))
    }))

    cluster.enrichment.df$q.value = p.adjust(cluster.enrichment.df$p.value, "fdr")

    cluster.enrichment.df$log2.odds = with(cluster.enrichment.df,
        ifelse(odds.ratio == 0, -5, log2(odds.ratio)))

    # Determine a list of cells that are informative in this CDS given qvalue threshold
    informative_targets = subset(cluster.enrichment.df, q.value <= qval_threshold & log2.odds > 0)
    informative_targets = with(informative_targets, paste(target, cluster))

    cds_id = with(pData(cds), paste(gene, Cluster))
    informative_cells = colnames(cds[, cds_id %in% informative_targets])

    return(list('informative_cells'=informative_cells, 'barcode_enrichment'=cluster.enrichment.df, 'processed_cds'=cds, 'guide_weights'=all_guide_weights))
}


# Helper function for DE
diff_fold_change_pseudo = function(X, id_to_average_over, normalize_to, pseudo_exprs=0.01){
  Group_subset = list()
  Grouping_list = unique(pData(X)[,id_to_average_over])
  Grouping_list = setdiff(Grouping_list, normalize_to)

  Norm_cds = X[,pData(X)[,id_to_average_over] == normalize_to]
  Norm_Size_Factor = as.matrix(exprs(Norm_cds))/pData(Norm_cds)$Size_Factor

  for (Group in Grouping_list) {

    temp_cds = X[,pData(X)[,id_to_average_over] == Group]
    temp_Size_Factor = as.matrix(exprs(temp_cds))/pData(temp_cds)$Size_Factor
    temp_foldChange = as.matrix(log2(rowMeans(temp_Size_Factor)/rowMeans(Norm_Size_Factor)))

    Group_subset[[which(Grouping_list == Group)]] = temp_foldChange
    names(Group_subset) = paste0("log2FC_",Grouping_list[seq_along(Group_subset)],sep ="")

    print(Group)
  }
  return(Group_subset)
}

# Theme for plots
theme_cfg <- function(base_size=12, font=NA, grid_lines=TRUE){
  result = theme(
    strip.background = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(face="bold")
    )
  if (grid_lines) {
    result = result + theme(
      panel.grid.major = element_line(colour = '#d3d3d3', size = 0.1, linetype = 'dotted'),
      panel.grid.minor = element_line(colour = '#d3d3d3', size = 0.1, linetype = 'dotted')
    )
  }

  return(result)
}

###### This function is based off of plot_genes_jitter from Monocle ######
plot_genes_violin = function (cds_subset, grouping = "State", min_expr = NULL, cell_size = 0.75,
  nrow = NULL, ncol = 1, panel_order = NULL, color_by = NULL,
  plot_trend = FALSE, label_by_short_name = TRUE, relative_expr = TRUE, log_scale = FALSE, show_combined=NULL)
{
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial",
    "negbinomial.size")) {
    integer_expression = TRUE
  }
  else {
    integer_expression = FALSE
    relative_expr = TRUE
  }
  if (integer_expression) {
    cds_exprs = exprs(cds_subset)
    if (relative_expr) {
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs = Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    }
        #cds_exprs = reshape2::melt(round(as.matrix(cds_exprs)))
        cds_exprs = reshape2::melt(as.matrix(cds_exprs))
      }
      else {
        cds_exprs = exprs(cds_subset)
        cds_exprs = reshape2::melt(as.matrix(cds_exprs))
      }
      if (is.null(min_expr)) {
        min_expr = cds_subset@lowerDetectionLimit
      }
      colnames(cds_exprs) = c("f_id", "Cell", "expression")
      cds_exprs$expression[cds_exprs$expression < min_expr] = min_expr
      cds_pData = pData(cds_subset)


    # Custom bit for adding in a group for
    if(! is.null(show_combined)) {
      for(combine_gene in show_combined) {
        cds_pData_all = subset(cds_pData, gene == combine_gene)
        cds_pData_all[, grouping] = paste("All", combine_gene)
        cds_pData = rbind(cds_pData, cds_pData_all)
      }
    }

    cds_fData = fData(cds_subset)
    cds_exprs = merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "id")
    cds_exprs = merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "cell")
    cds_exprs$adjusted_expression = log10(cds_exprs$expression)




    if (label_by_short_name == TRUE) {
      if (is.null(cds_exprs$gene_short_name) == FALSE) {
        cds_exprs$feature_label = cds_exprs$gene_short_name
        cds_exprs$feature_label[is.na(cds_exprs$feature_label)] = cds_exprs$f_id
      }
      else {
        cds_exprs$feature_label = cds_exprs$f_id
      }
    }
    else {
      cds_exprs$feature_label = cds_exprs$f_id
    }
    if (is.null(panel_order) == FALSE) {
      cds_exprs$feature_label = factor(cds_exprs$feature_label,
        levels = panel_order)
    }

	
    # This prints stats for each category
    # Required for Nature methods
    print(head(cds_exprs))
    print(tapply(cds_exprs$expression, paste(cds_exprs$feature_label, cds_exprs[, grouping]), summary))

    q = ggplot(aes_string(x = grouping, y = "expression"), data = cds_exprs)
    if (is.null(color_by) == FALSE) {
      q = q + geom_violin(aes_string(fill = color_by))
    }
    else {
      q = q + geom_violin()
    }
    if (plot_trend == TRUE) {
      q = q + stat_summary(fun.data = "mean_cl_boot",
        size = 0.2)
      q = q + stat_summary(aes_string(x = grouping, y = "expression",
        group = color_by), fun.data = "mean_cl_boot",
      size = 0.2, geom = "line")
    }
    q = q + facet_wrap(~feature_label, nrow = nrow,
      ncol = ncol, scales = "free_y")
    if (min_expr < 1) {
      q = q + expand_limits(y = c(min_expr, 1))
    }


    q = q + ylab("Expression") + xlab(grouping)

    if (log_scale == TRUE){

      q = q + scale_y_log10()
    }
    q
  }


################################################
# Other functions from Jose
################################################

diff_foldChange = function(X, id_to_average_over, normalize_to){
  Group_subset = list()
  Grouping_list = unique(pData(X)[,id_to_average_over])
  Grouping_list = setdiff(Grouping_list, normalize_to)

  Norm_cds = X[,pData(X)[,id_to_average_over] == normalize_to]
  Norm_Size_Factor = exprs(Norm_cds)/pData(Norm_cds)$Size_Factor

  for (Group in Grouping_list) {

    temp_cds = X[,pData(X)[,id_to_average_over] == Group]
    temp_Size_Factor = exprs(temp_cds)/pData(temp_cds)$Size_Factor
    temp_foldChange = as.matrix(log2(rowMeans(temp_Size_Factor)/rowMeans(Norm_Size_Factor)))

    Group_subset[[which(Grouping_list == Group)]] = temp_foldChange
    names(Group_subset) = paste0("log2FC_",Grouping_list[seq_along(Group_subset)],sep ="")

    print(Group)
  }
  return(Group_subset)
}

expressed_genes_clusterCutoff = function(X, gene_list, cutoff) {

  if (is.null(pData(X)$Cluster)) {
    stop("Error: to call this function you must assign genes to clusters using clusterCells_Density_Peak() first")
  }

  Cluster_list = unique(pData(X)$Cluster)
  expressed_in_cluster = list()

  for (cluster in Cluster_list) {

    expressed_over_cutoff = rowSums(exprs(X[gene_list,pData(X)$Cluster == cluster])) >
    (nrow(subset(pData(X),Cluster == cluster)) * cutoff)

    expressed_in_cluster[[which(Cluster_list == cluster)]] = expressed_over_cutoff

    print(cluster)
  }

  expressed_over_cutoff_AllClusters = do.call(cbind,expressed_in_cluster)
  expressed_over_cutoff_by_gene = rowSums(expressed_over_cutoff_AllClusters) > 1

  return(expressed_over_cutoff_by_gene)
}




plot_gsea_go<-function(gsa_res,
                       mode="distinct",
                       top_g_gene_sets=15,
                       q_thresh=0.05,
                       fill_color="black",
                       whitelist=NULL,
                       plot_q_vals=FALSE)
{
  resTab <- GSAsummaryTable(gsa_res)

  res_df <- resTab[,c("Name", "Genes (tot)", "Genes (up)", "Genes (down)")]

  if (mode == "distinct"){
    resTab <- resTab[with(resTab, order(-log10(resTab[,"p (dist.dir.up)"]))), ]
    res_df$stat_up <- resTab[,"Stat (dist.dir.up)"]
    res_df$stat_dn <- -resTab[,"Stat (dist.dir.dn)"]
    res_df$p_up <- resTab[,"p (dist.dir.up)"]
    res_df$p_dn <- resTab[,"p (dist.dir.dn)"]
    res_df$q_up <- resTab[,"p adj (dist.dir.up)"]
    res_df$q_dn <- resTab[,"p adj (dist.dir.dn)"]
  }
  else if (mode == "mixed"){
    resTab <- resTab[with(resTab, order(-log10(resTab[,"p (mix.dir.up)"]))), ]
    res_df$stat_up <- resTab[,"Stat (mix.dir.up)"]
    res_df$stat_dn <- -resTab[,"Stat (mix.dir.dn)"]
    res_df$p_up <- resTab[,"p (mix.dir.up)"]
    res_df$p_dn <- resTab[,"p (mix.dir.dn)"]
    res_df$q_up <- resTab[,"p adj (mix.dir.up)"]
    res_df$q_dn <- resTab[,"p adj (mix.dir.dn)"]
  }


  #if (is.null(whitelist)){
  to_draw_df_dn <- arrange(subset(res_df, q_dn < 1), desc(q_dn))
  to_draw_df_dn <- to_draw_df_dn[,c("Name", "Genes (tot)", "Genes (up)", "Genes (down)", "stat_dn", "p_dn", "q_dn")]
  colnames(to_draw_df_dn) <- c("Name", "Genes (tot)", "Genes (up)", "Genes (down)", "stat", "p_val", "q_val")
  to_draw_df_dn$q_val <- log10(to_draw_df_dn$q_val + 1e-3)

  to_draw_df_up <- arrange(subset(res_df, q_up < 1), desc(q_up))
  to_draw_df_up <- to_draw_df_up[,c("Name", "Genes (tot)", "Genes (up)", "Genes (down)", "stat_up", "p_up", "q_up")]
  colnames(to_draw_df_up) <- c("Name", "Genes (tot)", "Genes (up)", "Genes (down)", "stat", "p_val", "q_val")
  to_draw_df_up$q_val <- -log10(to_draw_df_up$q_val + 1e-3)

  res_df <- rbind(to_draw_df_up, to_draw_df_dn)

  res_df <- arrange(res_df, desc(q_val))
  res_df <- transform(res_df, Name=reorder(Name, q_val) )

  if (is.null(whitelist) == FALSE){
    g <- qplot(Name, q_val, data=res_df, size=I(0.65), color = Name %in% whitelist)

    up_gene_sets <- subset(res_df, stat > 0 & Name %in% whitelist)
    if (nrow(up_gene_sets) > 0){
      g <- g + geom_text(aes(Name, q_val, label=Name, hjust="right", check_overlap=TRUE), data=up_gene_sets)
    }
    down_gene_sets <- subset(res_df, stat < 0 & Name %in% whitelist)
    if (nrow(down_gene_sets) > 0){
      g <- g + geom_text(aes(Name, q_val, label=Name, hjust="left", check_overlap=TRUE), data=down_gene_sets)
    }

    print (subset(res_df, Name %in% whitelist))
  }
  else{
    g <- qplot(Name, q_val, data=res_df)
  }
  g <- g + theme(axis.text.x = element_blank()) +
    theme(axis.ticks.x = element_blank())
  g <- g + geom_hline(yintercept=-log10(q_thresh))
  g <- g + geom_hline(yintercept=log10(q_thresh))
  g <- g + scale_color_manual(values=c("black", "red"))
  g <- g + theme(legend.position="none")
  g <- g + theme(panel.grid.minor=element_line(colour=NA), panel.grid.major=element_line(colour=NA))

  g
}


loadGSCSafe <- function (file, type = "auto", addInfo, sep="\t", encoding="latin1")
{
  if (missing(addInfo)) {
    addUserInfo <- "skip"
    addInfo <- "none"
  }
  else {
    addUserInfo <- "yes"
  }
  tmp <- try(type <- match.arg(type, c("auto", "gmt", "sbml",
                                       "sif", "data.frame"), several.ok = FALSE), silent = TRUE)
  if (class(tmp) == "try-error") {
    stop("argument type set to unknown value")
  }
  if (type == "auto") {
    if (class(file) == "character") {
      tmp <- unlist(strsplit(file, "\\."))
      type <- tolower(tmp[length(tmp)])
      if (!type %in% c("gmt", "sif", "sbml", "xml"))
        stop(paste("can not handle .", type, " file extension, read manually using e.g. read.delim() and load as data.frame",
                   sep = ""))
    }
    else {
      type <- "data.frame"
    }
  }
  if (type == "gmt") {
    con <- file(file, encoding=encoding)
    tmp <- try(suppressWarnings(open(con)), silent = TRUE)
    if (class(tmp) == "try-error")
      stop("file could not be read")
    if (addUserInfo == "skip")
      addInfo <- vector()
    gscList <- list()
    i <- 1
    tmp <- try(suppressWarnings(while (length(l <- scan(con,
                                                        nlines = 1, what = "character", quiet = T, sep=sep)) > 0) {
      if (addUserInfo == "skip")
        addInfo <- rbind(addInfo, l[1:2])
      tmp <- l[3:length(l)]
      gscList[[l[1]]] <- unique(tmp[tmp != "" & tmp !=
                                      " " & !is.na(tmp)])
      i <- i + 1
    }), silent = TRUE)
    if (class(tmp) == "try-error")
      stop("file could not be read")
    close(con)
    gsc <- gscList[!duplicated(names(gscList))]
    if (addUserInfo == "skip")
      addInfo <- unique(addInfo)
  }
  else if (type %in% c("sbml", "xml")) {
    require(rsbml)
    tmp <- try(sbml <- rsbml_read(file))
    if (class(tmp) == "try-error") {
      stop("file could not be read by rsbml_read()")
    }
    gsc <- list()
    for (iReaction in 1:length(reactions(model(sbml)))) {
      metIDs <- names(c(reactants(reactions(model(sbml))[[iReaction]]),
                        products(reactions(model(sbml))[[iReaction]])))
      geneIDs <- names(modifiers(reactions(model(sbml))[[iReaction]]))
      if (length(geneIDs) > 0) {
        geneNames <- rep(NA, length(geneIDs))
        for (iGene in 1:length(geneIDs)) {
          geneNames[iGene] <- name(species(model(sbml))[[geneIDs[iGene]]])
        }
        for (iMet in 1:length(metIDs)) {
          gsc[[metIDs[iMet]]] <- c(gsc[[metIDs[iMet]]],
                                   geneNames)
        }
      }
    }
    if (length(gsc) == 0) {
      stop("no gene association found")
    }
    else {
      for (iMet in 1:length(gsc)) {
        tmp1 <- name(species(model(sbml))[[names(gsc)[iMet]]])
        tmp2 <- compartment(species(model(sbml))[[names(gsc)[iMet]]])
        names(gsc)[iMet] <- paste(tmp1, " (", tmp2, ")",
                                  sep = "")
      }
    }
  }
  else if (type == "sif") {
    tmp <- try(gsc <- as.data.frame(read.delim(file, header = FALSE,
                                               quote = "", as.is = TRUE), stringsAsFactors = FALSE),
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("argument file could not be read and converted into a data.frame")
    }
    if (ncol(gsc) != 3) {
      stop("sif file should contain three columns")
    }
    if (addUserInfo == "skip")
      addInfo <- gsc[, c(1, 2)]
    gsc <- gsc[, c(3, 1)]
    tmp <- nrow(gsc)
    gsc <- unique(gsc)
    geneSets <- unique(gsc[, 2])
    gscList <- list()
    for (iGeneSet in 1:length(geneSets)) {
      gscList[[iGeneSet]] <- gsc[gsc[, 2] == geneSets[iGeneSet],
                                 1]
    }
    names(gscList) <- geneSets
    gsc <- gscList
  }
  else if (type == "data.frame") {
    tmp <- try(gsc <- as.data.frame(file, stringsAsFactors = FALSE),
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("argument file could not be converted into a data.frame")
    }
    for (i in 1:ncol(gsc)) {
      gsc[, i] <- as.character(gsc[, i])
    }
    if (ncol(gsc) != 2) {
      stop("argument file has to contain exactly two columns")
    }
    tmp <- nrow(gsc)
    gsc <- unique(gsc)
    geneSets <- unique(gsc[, 2])
    gscList <- list()
    for (iGeneSet in 1:length(geneSets)) {
      gscList[[iGeneSet]] <- gsc[gsc[, 2] == geneSets[iGeneSet],
                                 1]
    }
    names(gscList) <- geneSets
    gsc <- gscList
  }
  if (addUserInfo == "yes") {
    tmp <- try(addInfo <- as.data.frame(addInfo, stringsAsFactors = FALSE),
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("failed to convert additional info in argument 'addInfo' into a data.frame")
    }
  }
  if (class(addInfo) == "data.frame") {
    if (ncol(addInfo) != 2)
      stop("additional info in argument 'file' or 'addInfo' has to contain 2 columns")
    tmp <- nrow(addInfo)
    addInfo <- unique(addInfo[addInfo[, 1] %in% names(gsc),
                              ])
  }
  else {
  }
  res <- list(gsc, addInfo)
  names(res) <- c("gsc", "addInfo")
  class(res) <- "GSC"
  return(res)
}

# Functions used for collecting and plotting GSEA
collect_gsa_hyper_results <- function(genes_list, testing_list, deg_results, qval, gsc, grouping_variable="target"){

    target_list <- testing_list
    gene_universe <- genes_list
    gsa_results <- list()

    for (target in target_list){

        deg_subset <- deg_results[deg_results[, grouping_variable] == target &
                                  deg_results$qval < qval,]
        deg_subset_genes <- unique(as.character(deg_subset$gene_short_name))

        gsaRes <- runGSAhyper(deg_subset_genes, gsc=gsc, universe=gene_universe, adjMethod = "BH")
        gsa_results[[which(target_list == target)]] <- gsaRes


        print(target)
    }
    names(gsa_results) <- target_list
    return(gsa_results)
}

plot_gsa_hyper_heatmap <- function(gsa_results, significance=0.05)
{
  hyper_df <- ldply(gsa_results, function(gsa_res)
  {
    data.frame(gene_set = names(gsa_res$pvalues), pval = gsa_res$pvalues, qval = gsa_res$p.adj)
  })
  colnames(hyper_df)[1] <- "Gene_KO"

  hyper_df <- subset(hyper_df, qval <= significance)
  hyper_df <- merge(hyper_df, ddply(hyper_df, .(gene_set), function(x) { nrow(x) }), by="gene_set")

  hyper_df$gene_set = str_split_fixed(hyper_df$gene_set, '%', n=3)[, 1]


  hyper_df$gene_set <- factor(hyper_df$gene_set, levels=unique(arrange(hyper_df, V1, Gene_KO)$gene_set))

  plot_object = ggplot(hyper_df, aes(Gene_KO, gene_set)) +
    geom_tile(aes(fill=-log10(qval)), color="black", size=0.5) +
    scale_fill_viridis(option="magma", name="-log10(qval)", na.value="white") +
    xlab('Target') +
    ylab('GO Term') +
    theme_cfg(grid_lines=FALSE) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  return(plot_object)

}

get_gsea_results_wrapper = function(pairwise_deg_results, column_value, gsc_object, direction=c("all", "positive", "negative"), qval_threshold=0.05, grouping_variable="target") {
  # Use list of all expressed genes as gene universe for GO analysis
  All_Ensembl_GSAlist = as.matrix(unique(pairwise_deg_results$gene_short_name))
  All_Ensembl_GSAlist = All_Ensembl_GSAlist[,1]
  All_Ensembl_GSAlist = unique(toupper(All_Ensembl_GSAlist))

  # Get the DEGs for specified column value
  significant_degs = pairwise_deg_results[pairwise_deg_results$qval < 0.05 & pairwise_deg_results$column_value == column_value,]
  significant_degs_ensembl_gsalist  = as.matrix(unique(significant_degs$gene_short_name))
  significant_degs_ensembl_gsalist = significant_degs_ensembl_gsalist[,1]
  significant_degs_ensembl_gsalist = toupper(significant_degs_ensembl_gsalist)

  # Do the test including genes that have a log fold change in the specified direction
  if (direction == "all") {
    # Lump everything that is up and down together
    deg_df = significant_degs
    testing_list = as.character(unique(deg_df[, grouping_variable]))
  } else if (direction == "positive") {
    # Look at positive and negative log fold changes separately
    deg_df = significant_degs[significant_degs$log2_fc > 0, ]
    testing_list = as.character(unique(deg_df[, grouping_variable]))
  } else if (direction == "negative") {
    deg_df = significant_degs[significant_degs$log2_fc < 0, ]
    testing_list = as.character(unique(deg_df[, grouping_variable]))
  }

  testing_list = toupper(testing_list)

  gsa_hyper = collect_gsa_hyper_results(genes_list = All_Ensembl_GSAlist, testing_list = testing_list, deg_results = deg_df, qval = qval_threshold, gsc = gsc_object, grouping_variable=grouping_variable)

  return(gsa_hyper)
}

