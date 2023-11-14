#!/usr/bin/env R

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos =  "http://cran.us.r-project.org")


for (pkg in c("devtools", "argparse", "Seurat",  "scater", "scuttle", "ggplot2")){
    if (! pkg %in% installed.packages()) {
        
        BiocManager::install(pkg, update = FALSE)

    }
}

## https://community.rstudio.com/t/seurat-functions-throwing-error-messages-due-to-conflicting-packages/176936/3
devtools::install_version("Matrix",version = "1.6.1.1")
