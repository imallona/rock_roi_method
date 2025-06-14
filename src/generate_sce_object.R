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

parser$add_argument('--captured_gtf', 
                    type = 'character',
                    help = 'Captured features GTF (path)')

args <- parser$parse_args()

get_captured_gene_ids <- function(gtf) {
    fd <- read.table(gtf, sep = '\t')
    if (ncol(fd) == 1)
        stop('Wrong GTF format')
    return(sapply(strsplit(fd$V9, split = ';'), function(x) return(gsub('gene_id ', '', x[1]))))
}

read_matrix <- function(mtx, cells, features, cell.column = 1, feature.column = 1, modality, wta_whitelist) {
  cell.barcodes <- read.table(
    file = cells,
    header = FALSE,
    row.names = cell.column)


  feature.names <- read.table(
    file = features,
    header = FALSE,
    row.names = feature.column)

  d <- as(readMM(mtx), 'CsparseMatrix')

  ## if (modality == 'wta') {
  ## colnames(d) <- gsub('_', '', rownames(cell.barcodes))
  ## } else if (modality == 'tso') {
  ##     ## remove the fixed parts of the TSO CBs
  ##     colnames(d) <- paste0(
  ##         substr(rownames(cell.barcodes), 1, 9),
  ##         substr(rownames(cell.barcodes), 9+4+1, 9+4+9),
  ##         substr(rownames(cell.barcodes), 9+4+9+4+1, 9+4+9+4+9))
  ## }

  colnames(d) <- rownames(cell.barcodes)
  rownames(d) <- rownames(feature.names)
  
  if (modality == 'tso') {
      d <- d[,wta_whitelist]
  }

  return(d)
}

read_featurecounts <- function(fn, wta_whitelist) {
    mtso <- read.table(fn, header = TRUE, sep = '\t')
    if (all(!duplicated(mtso$Geneid))) {
        rownames(mtso) <- mtso$Geneid
    } else {
        rownames(mtso) <- sprintf('%s_%s:%s-%s', mtso$Geneid, mtso$Chr, mtso$Start, mtso$End)
    }
    
    colnames(mtso) <- sapply(strsplit(split = '.bam.', colnames(mtso)), function(x) return(x[2]))

    mtso <- mtso[,grep('tso', colnames(mtso))]
    colnames(mtso) <- gsub('tso_', '', colnames(mtso))

    ## append 0s instead of NAs for non-detected CBs
    if (!all(wta_whitelist %in% colnames(mtso))) {
        empty <- matrix(nrow = nrow(mtso), ncol = length(setdiff(wta_whitelist, colnames(mtso))), data = 0)
        colnames(empty) <- setdiff(wta_whitelist, colnames(mtso))
        rownames(empty) <- rownames(mtso)
    
        mtso <- as(as.matrix(cbind(mtso, empty))[,wta_whitelist], 'CsparseMatrix')
    } else {
        mtso <- as(as.matrix(mtso)[,wta_whitelist], 'CsparseMatrix')
    }
    return(mtso)
}

## wta start

wd <- args$working_dir
id <- args$sample
gtf <- args$captured_gtf

wta <- read_matrix(mtx = file.path(wd, 'align_wta', id,  'Solo.out', 'Gene', 'filtered', 'matrix.mtx'),
                   cells = file.path(wd, 'align_wta', id,  'Solo.out', 'Gene', 'filtered', 'barcodes.tsv'),
                   features = file.path(wd, 'align_wta', id,  'Solo.out', 'Gene', 'filtered', 'features.tsv'),
                   cell.column = 1,
                   feature.column = 1,
                   modality = 'wta')

wta_feat <- read.table(file.path(wd, 'align_wta', id,  'Solo.out', 'Gene', 'filtered', 'features.tsv'),
                       row.names = 1,
                       header = FALSE)

captured <- get_captured_gene_ids(gtf)

colnames(wta_feat) <- c("name", "type", "value")

wta_feat$captured <- ifelse(wta_feat$name %in% captured, yes = 'captured', no = 'not_captured')

## wta end

if (args$run_mode == 'tso ontarget multi') {
    mtso <- read_featurecounts(file.path(wd, 'multimodal', id, 'tso_featurecounted'),
                               wta_whitelist = colnames(wta))
    mwta <- read_featurecounts(file.path(wd, 'multimodal', id, 'wta_featurecounted'),
                               wta_whitelist = colnames(wta))
    
    (sce <- SingleCellExperiment(assays = list(wta = wta),
                                 altExps = list(wta_ontarget_multi = SummarizedExperiment(mwta),
                                                tso_ontarget_multi = SummarizedExperiment(mtso)),
                                 mainExpName = id,
                                 rowData = wta_feat))

} else if (args$run_mode == 'tso off- and ontarget unique') {
    utso <- read_matrix(mtx = file.path(wd, 'align_tso', id,  'Solo.out', 'Gene', 'raw', 'matrix.mtx'),
                   cells = file.path(wd, 'align_tso', id,  'Solo.out', 'Gene', 'raw', 'barcodes.tsv'),
                   features = file.path(wd, 'align_tso', id,  'Solo.out', 'Gene', 'raw', 'features.tsv'),
                   cell.column = 1,
                   feature.column = 1,
                   modality = 'tso',
                   wta_whitelist = colnames(wta))


    (sce <- SingleCellExperiment(assays = list(wta = wta, tso_off_and_ontarget_unique = utso),
                            mainExpName = id,
                            rowData = wta_feat))
} else if (args$run_mode == 'all') {

    mtso <- read_featurecounts(file.path(wd, 'multimodal', id, 'tso_featurecounted'),
                               wta_whitelist = colnames(wta))
    mwta <- read_featurecounts(file.path(wd, 'multimodal', id, 'wta_featurecounted'),
                               wta_whitelist = colnames(wta))
    
    utso <- read_matrix(mtx = file.path(wd, 'align_tso', id,  'Solo.out', 'Gene', 'raw', 'matrix.mtx'),
                   cells = file.path(wd, 'align_tso', id,  'Solo.out', 'Gene', 'raw', 'barcodes.tsv'),
                   features = file.path(wd, 'align_tso', id,  'Solo.out', 'Gene', 'raw', 'features.tsv'),
                   cell.column = 1,
                   feature.column = 1,
                   modality = 'tso',
                   wta_whitelist = colnames(wta))

    (sce <- SingleCellExperiment(assays = list(wta = wta, tso_off_and_ontarget_unique = utso),
                                 altExps = list(wta_ontarget_multi = SummarizedExperiment(mwta),
                                                tso_ontarget_multi = SummarizedExperiment(mtso)),
                                 mainExpName = id,
                                 rowData = wta_feat))
}

saveRDS(object = sce, file = args$output_fn)
