#!/usr/bin/env Rscript
##
## Generates a SCE object with the relevant TSO/WTA data
##
## Izaskun Mallona
## Started Oct 18th 2023

suppressPackageStartupMessages( {
    library(SingleCellExperiment)
    library(argparse)
    library(Matrix)
})

parser <- ArgumentParser(description='Builds a rock and roi SingleCellExperiment object for a given sample.')

parser$add_argument('--sample',
                    type = "character",
                    help = 'Sample identifier')

parser$add_argument('--run_mode', 
                    type = 'character',
                    help = 'Run mode')

parser$add_argument('--working_dir', 
                    type = 'character',
                    help = 'Working directory')

parser$add_argument('--output_fn', 
                    type = 'character',
                    help = 'Output SCE filename (path)')

args <- parser$parse_args()

read_matrix <- function(mtx, cells, features, cell.column = 1, feature.column = 1) {
  cell.barcodes <- read.table(
    file = cells,
    header = FALSE,
    row.names = cell.column)


  feature.names <- read.table(
    file = features,
    header = FALSE,
    row.names = feature.column)

  d <- readMM(mtx)

  colnames(d) <- rownames(cell.barcodes)
  rownames(d) <- rownames(feature.names)
  ## d <- as(d, "CSparseMatrix")

  return(d)
}


## wta start

wd <- args$working_dir
id <- args$sample

wta <- read_matrix(mtx = file.path(wd, 'align_wta', id,  'Solo.out', 'Gene', 'filtered', 'matrix.mtx'),
                   cells = file.path(wd, 'align_wta', id,  'Solo.out', 'Gene', 'filtered', 'barcodes.tsv'),
                   features = file.path(wd, 'align_wta', id,  'Solo.out', 'Gene', 'filtered', 'features.tsv'),
                   cell.column = 1,
                   feature.column = 1)

wta_feat <- read.table(file.path(wd, 'align_wta', id,  'Solo.out', 'Gene', 'filtered', 'features.tsv'),
                       row.names = 1,
                       header = FALSE)
      
## wta end


(sce <- SingleCellExperiment(assays = list(wta = wta),
                            mainExpName = id,
                            rowData = wta_feat))

saveRDS(object = sce, file = args$output_fn)
