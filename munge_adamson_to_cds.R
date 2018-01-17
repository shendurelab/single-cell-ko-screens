library(argparse)
library(Matrix)
library(stringr)
library(monocle)

.adamson_to_cds = function(barcodes, genes, matrix, cell_identities, treatment) {
  # Quick helper function to munge data from Adamson et al. into a CDS object that
  # allows for use with some of our scripts.
  matrix_path = matrix
  genes_path = genes
  barcodes_path = barcodes

  if( ! file.exists(matrix_path) ) { stop(paste("Expression matrix not found in 10X output:", matrix_path)) }
  if( ! file.exists(genes_path) ) { stop(paste("Genes file not found in 10X output:", genes_path)) }
  if( ! file.exists(barcodes_path) ) { stop(paste("Barcodes file not found in 10X output:", barcodes_path)) }

  # All files exist, read them in
  matrix = Matrix::readMM(matrix_path)
  barcodes = read.table(barcodes_path, header=F, as.is=T)[,1]
  gene_table = read.table(genes_path, header=F, as.is=T)

  genes = gene_table[, 1]

  # Add cell names to matrix
  colnames(matrix) = barcodes

  # Read in assignments and filter to the confident set
  assignments = read.delim(cell_identities, sep=',', header=T)
  assignments$gene = str_split_fixed(assignments$guide.identity, '_', n=2)[, 1]
  assignments$cell = assignments$cell.BC
  assignments = subset(assignments, good.coverage == TRUE & number.of.cells == 1)
  rownames(assignments) = assignments$cell

  # Construct metadata table that includes directory samples come from and other stats
  total_umis = colSums(matrix)

  metadata_df = data.frame(
    cell = colnames(matrix),
    total_umis = total_umis
    )

  metadata_df = merge(metadata_df, assignments, by="cell")
  row.names(metadata_df) = metadata_df$cell

  labeled_cells = intersect(row.names(metadata_df), colnames(matrix))
  matrix = matrix[, labeled_cells]
  metadata_df = metadata_df[labeled_cells, ]

  colnames(gene_table) = c("id", "gene_short_name")
  row.names(gene_table) = gene_table$id

  pd = new("AnnotatedDataFrame", data = metadata_df)
  fd = new("AnnotatedDataFrame", data = gene_table)
  cds = newCellDataSet(matrix,
   phenoData=pd,
   featureData=fd,
   expressionFamily=negbinomial.size(),
   lowerDetectionLimit=0.5)

  pData(cds)$treatment = treatment
  pData(cds)$condition = treatment
  pData(cds)$guide_count = 1
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  return(cds)
}

parser = argparse::ArgumentParser(description='Script to munge Adamson data into a monocle cds.')
parser$add_argument('--matrix', help='MTX file')
parser$add_argument('--barcodes', help='Barcodes file')
parser$add_argument('--genes', help='Genes file')
parser$add_argument('--assignments', help='File with assignments')
parser$add_argument('--output_file', help='RDS formatted file.')
parser$add_argument('--treatment', help='Value to provide for treatment/condition column.')
args = parser$parse_args()

cds = .adamson_to_cds(args$barcodes, args$genes, args$matrix, cell_identities=args$assignments, args$treatment)
saveRDS(cds, file=args$output_file)
