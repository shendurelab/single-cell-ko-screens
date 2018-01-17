library(stringr)
library(monocle)
library(plyr)
library(dplyr)
library(reshape2)
library(argparse)
source('/net/trapnell/vol1/ajh24/proj/2017emt_cfg/bin/cfg_helper_functions.R')

parser = argparse::ArgumentParser(description="Script to generate and save DEG results between two gene KOs for a given sample condition.")
parser$add_argument('cds', help="CDS to load (RDS format).")
parser$add_argument('gene_1', help="First gene KO to test.")
parser$add_argument('gene_2', help="Second gene KO to test.")
parser$add_argument('output_file', help="Output file to save DEG table.")
parser$add_argument('--column_name', help="Column used to partion samples.")
parser$add_argument('--column_value', help="Test cells from matching this value for the given column_name.")
parser$add_argument('--informative', action="store_true", help="Flag to filter only to the set of cells labeled as informative.")
args = parser$parse_args()

getDEG<-function(cds, column_name, column_value, gene_1, gene_2, cell_detection_threshold, informative=FALSE){
      # Get the sample
      if (is.null(column_name) | is.null(column_value)) {
        de_df = cds
      } else {
        de_df = cds[, pData(cds)[, column_name] == column_value]
      }

      # Estimate size factors and dispersions using all the cells from that sample
      de_df = detectGenes(de_df, min_expr = 1)
      de_df = estimateSizeFactors(de_df)

      # Now subset the cells with the genes we care about
      genes_to_include = c(gene_1, gene_2)
      reduced_model = "~1"
      full_model = paste("~gene")

      # If requested, filter only to the informative set of genes
      if (informative) {
        de_df = de_df[, pData(de_df)$gene %in% genes_to_include & pData(de_df)$informative]
      } else {
        de_df = de_df[, pData(de_df)$gene %in% genes_to_include]
      }

      de_df = detectGenes(de_df, min_expr = 1)
      expressed_genes = row.names(fData(de_df)[rowSums(exprs(de_df) > 0) > cell_detection_threshold,])
      de_df = estimateDispersions(de_df)
      
      closeAllConnections()
      differential_test = differentialGeneTest(de_df[expressed_genes,], 
                                               fullModelFormulaStr=full_model,
                                               reducedModelFormulaStr=reduced_model,
                                               cores=4)

      differential_test$full_model = full_model
      differential_test$reduced_model = reduced_model
      differential_test$column = column_name
      differential_test$column_value = column_value
      differential_test$target = gene_1
      differential_test$versus = gene_2

      differential_test$log2_fc=unlist(diff_fold_change_pseudo(de_df[row.names(differential_test), ], "gene", gene_2))

      return(differential_test)
}

print('Reading data...')
aggregated_cds = readRDS(args$cds)
aggregated_cds = aggregated_cds[, (!is.na(pData(aggregated_cds)$gene)) & pData(aggregated_cds)$guide_count == 1 & pData(aggregated_cds)$proportion > 0.8]

print('Getting DEG...')
NTC_target_diff_test = getDEG(aggregated_cds, args$column_name, args$column_value, args$gene_1, args$gene_2, cell_detection_threshold=50, informative=args$informative)

print('DEG complete.')

print('Writing file...')
write.table(NTC_target_diff_test, quote=F, row.names=F, file=args$output_file, sep="\t")
