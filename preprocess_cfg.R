library(dplyr)
library(Matrix)
library(monocle)
library(stringr)
library(argparse)
library(ggplot2)
library(gridExtra)
library(scales)
library(viridis)
library(data.table)
library(readr)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')

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

    if( ! file.exists(matrix_path) ) { stop(paste("Expression matrix not found in 10X output:", matrix_path)) }
    if( ! file.exists(genes_path) ) { stop(paste("Genes file not found in 10X output:", genes_path)) }
    if( ! file.exists(barcodes_path) ) { stop(paste("Barcodes file not found in 10X output:", barcodes_path)) }

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
    total_umis = Matrix::colSums(matrix)

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
  combined_expression_matrix = do.call(cbind, expression_matrices)
  row.names(combined_expression_matrix) = gene_table[, 1]

  combined_metadata_df = dplyr::bind_rows(metadata_dfs)
  row.names(combined_metadata_df) = combined_metadata_df$cell

  colnames(gene_table) = c("id", "gene_short_name")
  row.names(gene_table) = gene_table$id

  pd = new("AnnotatedDataFrame", data = as.data.frame(combined_metadata_df))
  fd = new("AnnotatedDataFrame", data = gene_table)
  cds = newCellDataSet(combined_expression_matrix,
   phenoData=pd,
   featureData=fd,
   expressionFamily=VGAM::negbinomial.size(),
   lowerDetectionLimit=0.5)
  return(cds)
}

# Functions for preprocessing of KO barcode data
# Can use with UMIs (umi = TRUE) or reads (UMI=FALSE)
process_barcodes = function(barcodes, umis=FALSE) {
	if (umis) {
		barcodes$count_column = barcodes$umi_count
	} else {
		barcodes$count_column = barcodes$read_count
	}

  # Prefiltering
  barcodes = dplyr::as_tibble(barcodes)
  barcodes = barcodes %>% dplyr::filter(count_column > 0)
  barcodes = barcodes %>% filter(!grepl('unprocessed', barcode))

	barcodes = barcodes %>%
              dplyr::group_by(cell) %>%
              dplyr::mutate(proportion=count_column / sum(count_column)) %>%
              dplyr::ungroup()

	# Make sure case matches
	barcodes = barcodes %>% mutate(barcode=toupper(barcode))
	return(barcodes)
}

get_filtered_barcodes = function(barcodes, guide_metadata, reads_threshold=2, proportion_threshold=0.05, umis=FALSE) {
  if (umis) {
          barcodes$count_column = barcodes$umi_count
  } else {
          barcodes$count_column = barcodes$read_count
  }

  barcodes.initial = barcodes %>%
    dplyr::filter(proportion > proportion_threshold & count_column > reads_threshold)

  guide_metadata = guide_metadata %>% mutate(guide=toupper(guide))
  barcodes.initial = dplyr::inner_join(barcodes.initial, guide_metadata, by=c("barcode"="guide"))

  barcodes.filtered = barcodes.initial %>%
    dplyr::arrange(gene, barcode) %>%
    dplyr::group_by(cell) %>%
    dplyr::summarize(gene = paste0(unique(gene), collapse = "_"), all_gene=paste0(gene, collapse='_'), barcode = paste0(barcode, collapse = "_"), read_count = sum(read_count), umi_count = sum(umi_count), proportion = sum(proportion), guide_count = n())

  # Make sure that the original metadata is preserved
  barcodes.filtered = dplyr::inner_join(barcodes.filtered, unique(barcodes.initial[, !colnames(barcodes.initial) %in% c("gene", "all_gene", "barcode", "read_count", "umi_count", "proportion", "guide_count", "count_column")]))
  return(barcodes.filtered)
}

## helper function to take initial gene column in pdata and expand to a set of columns (one for each KO)
## that indicates whether or not that cell had a given KO
## useful for some DE tests
expand_genotypes_to_indicator = function(pdata, guide_metadata) {
	all_kos = unique(guide_metadata$gene)

	tx = lapply(all_kos, function(x) {
    grepl(pattern = paste0('(^|_)', x, '(_|$)'), x = pdata$gene)
    })

	tx_merged = as.data.frame(do.call(cbind, tx))
	colnames(tx_merged) = all_kos

	new_pData = cbind(pdata, tx_merged)
}

flag_low_size_factor_clusters = function(cds, cell_detection_threshold, log2_size_factor_threshold) {
  # Preprocess
  expressed_genes = row.names(fData(cds)[rowSums(exprs(cds) > 0) > cell_detection_threshold,])
  cds = cds[expressed_genes, ]

  cds = detectGenes(cds, min_expr = 1)
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)

  # Reduce the dimensionality with a reasonable number of components
  cds = reduceDimension(cds[expressed_genes, ], max_components=2, norm_method = 'log', num_dim = 15 , reduction_method = 'tSNE', verbose = T)

  # Cluster with lax density peak parameters
  cds = clusterCells(cds, method = "densityPeak")

  # Flag clusters with low average size factor
  size_factor_score_by_cluster = pData(cds) %>% dplyr::group_by(Cluster) %>% dplyr::mutate(score = mean(log2(Size_Factor))) %>% dplyr::ungroup()
  pData(cds)$low_size_factor_cell_score = size_factor_score_by_cluster$score
  pData(cds)$low_size_factor_cell = pData(cds)$low_size_factor_cell_score < log2_size_factor_threshold

  low_quality_cells = row.names(pData(cds))[size_factor_score_by_cluster$score < log2_size_factor_threshold]

  tsne_size_factors_plot = plot_cell_clusters(cds, color_by = 'log2(Size_Factor)', cell_size = 0.5) +
  theme(axis.text = element_text(face="bold"), axis.title = element_text(face="bold"), panel.margin = unit(1, "lines")) +
  scale_color_viridis(option='magma')

  cds.scores_by_cluster = pData(cds) %>% dplyr::group_by(Cluster) %>% dplyr::summarize(score = mean(low_size_factor_cell_score))

  score_by_cluster_plot = ggplot(cds.scores_by_cluster, aes(Cluster, score, color=score < log2_size_factor_threshold)) +
  geom_point(size=2) +
  geom_linerange(ymin=0, aes(ymax=score)) +
  geom_hline(yintercept=log2_size_factor_threshold, color="red", size=1, linetype="dashed") +
  theme_classic() +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="#d3d3d3")) +
  theme(axis.text = element_text(face="bold"), axis.title = element_text(face="bold"), panel.margin = unit(1, "lines")) +
  ylab('log2(average size factor)') +
  guides(color=FALSE)

  tsne_by_flag_plot = plot_cell_clusters(cds, color_by = 'low_size_factor_cell', cell_size = 0.5) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="#d3d3d3")) +
  theme(axis.text = element_text(face="bold"), axis.title = element_text(face="bold"), panel.margin = unit(1, "lines")) +
  guides(color=FALSE)

  return(list("low_size_factor_cells"=low_quality_cells, "tsne_size_factors_plot"=tsne_size_factors_plot, "score_by_cluster_plot"=score_by_cluster_plot, "tsne_by_flag_plot"=tsne_by_flag_plot))
}

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

parser = argparse::ArgumentParser(description="Script to preprocess CFG data into a CDS object.")
parser$add_argument('sample_metadata', help='Table with sample metadata.')
parser$add_argument('output_cds', help='Output with CDS.')
parser$add_argument('output_pdata', help='Output with just pData of CDS.')
parser$add_argument('--genome', default='hg19', help='Genome used with cellranger to process data. Default is hg19')
parser$add_argument('--barcode_enrichment_qc_plot', required=FALSE, help='A file to output QC plot of barcode enrichment libraries.')
parser$add_argument('--guide_metadata', required=TRUE, help='A table of each guide and the gene it is associated with and any other metadata, must have "gene" and "guide" columns.')
parser$add_argument('--aggregated', action='store_true', help='Flag to indicate that the run directory specified in sample_directory column is an aggregated run. Must be same for all samples, samples must be in same order as aggregate samplesheet.')
parser$add_argument('--umis', action='store_true', help='Base thresholds and plots on UMI counts rather than the default read counts when assigning barcodes to cells.')
parser$add_argument('--ko_assignment_reads_threshold', type='double', default=10, help='The number of reads a guide must have to be considered a valid observation.')
parser$add_argument('--ko_assignment_proportion_threshold', type='double', default=0.075, help='The total proportion of data from a cell a guide must make up to be considered a valid guide observation.')
parser$add_argument('--ko_barcode_total_proportion_threshold', type='double', default=0.8, help='Total proportion all guides must make up to be considered a valid assignment.')
parser$add_argument('--cell_detection_threshold', type='integer', default=50, help='The number of cells a gene must be detected in to be considered for dim reduction in flagging low size factor cells.')
parser$add_argument('--log2_size_factor_threshold', type='double', default=-0.85, help='The lower threshold on average size factor for any cluster if --size_factor_filter is turned on.')
parser$add_argument('--size_factor_filter', action='store_true', help='Turns on filtering of clusters that have low size factors on average.')
parser$add_argument('--genotype_indicator_columns', action='store_true', help='This generates a TRUE/FALSE column for every individual genotype from being generated.')

args = parser$parse_args()

# Load in metadata (sample_directory, ko_barcode_file, ...)
cat('Loading metadata...\n')
sample_metadata = utils::read.delim(args$sample_metadata, header=T, sep='\t')
if (! all(c('sample_directory', 'ko_barcode_file') %in% colnames(sample_metadata))) {
  stop('Sample metadata file must have sample_directory and ko_barcode_file columns.')
}

sample_metadata$sample = unlist(lapply(sample_metadata$sample_directory, basename))
sample_directories = sample_metadata$sample_directory


# Load in data about mutants
guide_metadata = read.delim(args$guide_metadata, header=F, col.names=c('gene', 'guide'))
guide_metadata$gene[guide_metadata$gene == "NonTargetingControlGuideForHuman"] = "NONTARGETING"

duplicated_guides = duplicated(guide_metadata$guide)
duplicated_rows = duplicated(guide_metadata)

if (! all(duplicated_guides == duplicated_rows)) {
    stop('Your guide_metadata sheet should have all unique rows, but does not. Amongst the duplicate sequences, at least one has discrepancies, so there are conflicting entries that you need to fix. Please correct conflicts and rerun.')
} else if (sum(duplicated_rows) > 0) {
    cat("WARNING: guide_metadata sheet should have all unique rows but does not. However, in this csae the entire row was duplicated so there are no conflicts. Running with the duplicates removed safely.\n")
    guide_metadata = guide_metadata[!duplicated_rows, ]
}

# Load in barcodes
cat('Loading get_barcodes output...\n')
ko_barcodes = lapply(1:nrow(sample_metadata), function(i) {
  filename = sample_metadata$ko_barcode_file[i]
  cat(paste0('    Loading file ', i, ' of ', nrow(sample_metadata), '...\n'))
  return(readr::read_delim(filename, delim='\t'))
})

# Now process the barcodes
# Now process the barcodes
cat('Processing get_barcodes output...\n')
guide_metadata = dplyr::as_tibble(guide_metadata)

ko_barcodes = lapply(1:length(ko_barcodes), function(i) {
  cat(paste0('    processing barcode set ', i, ' of ', length(ko_barcodes), '...\n'))
  x = ko_barcodes[[i]]
  process_barcodes(x, umis=args$umis)
})

## Add on sample prefix to make cell names match the CDS
if (! args$aggregated) {
  cat('Prefixing cell IDs...\n')
  ko_barcodes = lapply(1:length(ko_barcodes), function(i) {
    x = ko_barcodes[[i]]
    x = x %>% mutate(cell = paste(cell, sample_metadata[i, "sample"], sep='_'),
                     sample = sample_metadata[i, "sample"])
    return(x)
  })

} else {
  # Make everything match the aggregated output
  sample_metadata$sample = paste(sample_metadata$sample, '-', 1:nrow(sample_metadata), sep='')
  ko_barcodes = lapply(1:length(ko_barcodes), function(i) {
    x = ko_barcodes[[i]]
    x = x %>% mutate(cell = stringr::str_replace(cell, '-1', ''),
                     cell = paste(cell, '-', i, sep=''),
                     sample = sample_metadata[i, "sample"])
    return(x)
  })
}

# Load in the RNA-seq data
cat('Loading RNA-seq data...\n')
if (! args$aggregated) {
  aggregated_cds = .tenx_to_cds(sample_directories, genome=args$genome)
} else {
  # User should have only specified one...
  sample_directories = unique(sample_directories)
  if (length(sample_directories) > 1) { stop('More than one sample directory specified for aggregated run, not allowed...')}
  aggregated_cds = .tenx_to_cds(unique(sample_directories), genome=args$genome)

  # Match cell names from KO barcodes
  new_cells = stringr::str_replace(colnames(aggregated_cds), paste('_', basename(sample_directories), sep=''), '')
  colnames(aggregated_cds) = new_cells
  row.names(pData(aggregated_cds)) = new_cells
  aggregated_cds$cell = new_cells
  aggregated_cds$sample = paste(aggregated_cds$sample, str_extract(aggregated_cds$cell, '-[0-9]+'), sep='')
}


## Generate QC plot for barcode enrichment
cat('Generating QC plot...\n')
if (! is.null(args$barcode_enrichment_qc_plot)) {
    ## Generate QC plot for barcode enrichment
    cat('    Combining KO barcodes across samples...\n')
    combined_barcodes_plot_df = dplyr::bind_rows(ko_barcodes)

    if (args$umis) {
        combined_barcodes_plot_df$count_column = combined_barcodes_plot_df$umi_count
    } else {
        combined_barcodes_plot_df$count_column = combined_barcodes_plot_df$read_count
    }

    combined_barcodes_plot_df = combined_barcodes_plot_df %>%
    		dplyr::group_by(sample, cell) %>%
				dplyr::mutate(valid=proportion > args$ko_assignment_proportion_threshold & count_column > args$ko_assignment_reads_threshold) %>%
				dplyr::mutate(guide_count = sum(valid)) %>%
				ungroup()

    combined_barcodes_plot_df$category = 'single guide'
    combined_barcodes_plot_df$category[!combined_barcodes_plot_df$cell %in% aggregated_cds$cell] = 'background barcode'
    combined_barcodes_plot_df$category[combined_barcodes_plot_df$guide_count < 1 & combined_barcodes_plot_df$category != 'background barcode'] = 'unassigned'
    combined_barcodes_plot_df$category[combined_barcodes_plot_df$guide_count > 1 & combined_barcodes_plot_df$category != 'background barcode'] = 'multiple guides'

    cat('    Plotting...\n')
    # Plot, excluding singletons if doesn't conflict with the user set threshold
    ggplot(combined_barcodes_plot_df %>% filter(count_column > min(args$ko_assignment_reads_threshold, 1)), aes(x=proportion, y=log10(count_column))) +
        geom_point(data=subset(combined_barcodes_plot_df, category == 'background_barcode'), size=0.15, alpha=0.5, aes(color=category)) +
        geom_point(data=subset(combined_barcodes_plot_df, category != 'background_barcode'), size=0.15, alpha=0.5, aes(color=category)) +
        geom_vline(xintercept=args$ko_assignment_proportion_threshold, color="red", size=1) +
        geom_hline(yintercept=log10(args$ko_assignment_reads_threshold), color="red", size=1, alpha=0.5) +
        facet_wrap(~sample, ncol=4) +
        xlab('barcode proportion') +
        ylab(paste('log10(', ifelse(args$umis, 'umis', 'reads'), ' from barcode)', sep='')) +
        theme_cfg() +
        theme(panel.margin = unit(1, "lines")) +
        scale_color_manual(values=c('background barcode'='#d3d3d3', 'multiple guides'='#51A7F9', 'single guide'='black', 'unassigned'='#F39019'), guide=FALSE) +
        ggsave(args$barcode_enrichment_qc_plot, height=5, width=8)
}

cat('Making barcode assignments...\n')
# Make KO assignments
ko_barcodes.filtered = lapply(ko_barcodes, get_filtered_barcodes, guide_metadata, reads_threshold = args$ko_assignment_reads_threshold, proportion_threshold = args$ko_assignment_proportion_threshold, umis=args$umis)

## Also make sure any metadata not calculated in filtered set is present in final dataframe
#ko_barcodes.filtered = lapply(1:length(ko_barcodes.filtered), function(i) { merge(ko_barcodes.filtered[[i]], unique(ko_barcodes[[i]][, !colnames(ko_barcodes[[i]]) %in% c("gene", "all_gene", "barcode", "read_count", "umi_count", "proportion", "guide_count")])) })

# Make a final set of assignments
all_assignments = dplyr::bind_rows(ko_barcodes.filtered)

aggregate_samples_pdata = merge(pData(aggregated_cds), all_assignments, all.x=T)
aggregate_samples_pdata = merge(aggregate_samples_pdata, sample_metadata)


# Expand pdata to include a TRUE/FALSE column for every individual KO
if (args$genotype_indicator_columns) {
    aggregate_samples_pdata = expand_genotypes_to_indicator(aggregate_samples_pdata, guide_metadata)
}

# Replace original pData
cat('Saving final output...\n')
row.names(aggregate_samples_pdata) = aggregate_samples_pdata$cell
pData(aggregated_cds) = aggregate_samples_pdata[row.names(pData(aggregated_cds)), ]

##################################################################
# Remove clusters of cells with much lower average size factors
##################################################################
# Subset out each sample
if ( args$size_factor_filter) {
  cat('Finding and removing low size factor clusters...\n')
	low_size_factor_results = lapply(unique(pData(aggregated_cds)$sample), function(x) {
		cds_subset = aggregated_cds[, pData(aggregated_cds)$sample == x]
		flag_low_size_factor_clusters(cds_subset, args$cell_detection_threshold, args$log2_size_factor_threshold)
		})

	# Get the names of these cells to filter final dataset
	low_size_factor_cells = unlist(lapply(low_size_factor_results, function(x) x$low_size_factor_cells))
} else {
	low_size_factor_cells = c()
}

##################################################################
# Merge info from enrichment with final CDS and save
##################################################################
# Filter the low size-factor cluster out of the final CDS entirely
final_cds = aggregated_cds[, ! row.names(pData(aggregated_cds)) %in% low_size_factor_cells]

# Save filtered and flagged CDS
saveRDS(final_cds, file=args$output_cds)
write.table(pData(final_cds), quote=F, row.names=F, sep="\t", file=args$output_pdata)
