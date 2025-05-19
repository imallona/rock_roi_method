#!/usr/bin/env R
##
## Generate a plain text unstructured report with relevant alignment metrics, per experiment (all samples)
##
## R --no-echo --no-restore --file=generate_mapping_report.R --args --path {path}
## Izaskun Mallona
## GPLv3
##
## started 1st Dec 2023


suppressPackageStartupMessages({
    library("argparse")
})

parser <- ArgumentParser()

parser$add_argument("-p", "--path", type="character", default='none',
                    help="Path to the output folder",
                    metavar="path")

tryCatch({
    args <- parser$parse_args()
}, error = function(x) parser$print_help())

## parse wta logs here

## args <- commandArgs(trailingOnly=TRUE)

## p <- '/home/imallona/ebrunner_spectral/mixing'
## p <- args[1]
p <- args$path

samples <- list.files(file.path(p, 'align_wta'))

wta_fns <- list.files(file.path(p, 'align_wta'), pattern = 'Summary.csv', recursive = TRUE)

wta <- lapply(wta_fns, function(x) read.csv(file.path(p, 'align_wta', x), header = FALSE))
names(wta) <- samples

## parse tso logs here

tso_fns <- list.files(file.path(p, 'align_tso'), pattern = 'Summary.csv', recursive = TRUE)

tso <- lapply(tso_fns, function(x) read.csv(file.path(p, 'align_tso', x), header = FALSE))
names(tso) <- samples

## parse self featurecounts here tso

wta_fc_fns <- list.files(file.path(p, 'multimodal'), pattern =  'tso_featurecounted.summary', recursive = TRUE)

if (length(wta_fc_fns) > 0) {
    wta_fc <- lapply(wta_fc_fns, function(x) read.table(file.path(p, 'multimodal', x),
                                                    header = TRUE, row.names = 1))
    names(wta_fc) <- dirname(wta_fc_fns)
}

tso_fc_fns <- list.files(file.path(p, 'multimodal'), pattern =  'tso_featurecounted.summary', recursive = TRUE)

if (length(tso_fc_fns) > 0) {
    tso_fc <- lapply(tso_fc_fns, function(x) read.table(file.path(p, 'multimodal', x),
                                                    header = TRUE, row.names = 1))
names(tso_fc) <- dirname(tso_fc_fns)
}

## assemble (rather dirty...)

for (sample in samples) {
    cat('# ', sample, '\n')
    
    cat('## STARsolo summary- WTA modality\n')
    print(wta[[sample]])

    cat('\n\n## STARsolo summary - TSO modality\n')
    print(tso[[sample]])

    if (length(wta_fc_fns) > 0) {
        cat('\n\n## Custom counting summary - WTA modality\n')
        cat('### Number of features\n')
        nrow(wta_fc)
        cat('### Alignment counting summary\n')
        print(apply(wta_fc[[sample]], 1, sum))
        cat('### Number of features\n')
        nrow(wta_fc[[sample]])
        cat('### Number of cells with at least one read assigned\n')   
        sum(wta_fc[[sample]]['Assigned',] > 0)
    }

    if (length(tso_fc_fns) > 0) {
        cat('\n\n## Custom counting summary - TSO modality\n')
        cat('### Number of features\n')
        nrow(tso_fc)
        cat('### Alignment counting summary\n')
        print(apply(tso_fc[[sample]], 1, sum))
        cat('### Number of features\n')
        nrow(tso_fc[[sample]])
        cat('### Number of cells with at least one read assigned\n')   
        print(sum(tso_fc[[sample]]['Assigned',] > 0))
    }
    cat('\n\n')
}
