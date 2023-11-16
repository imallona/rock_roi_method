#!/usr/bin/env R

options(repos=c(CRAN="http://cloud.r-project.org"))

print(R.Version())
print(.libPaths())

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", contriburl =  "http://cloud.r-project.org")

BiocManager::install('devtools', update = FALSE)

## https://community.rstudio.com/t/seurat-functions-throwing-error-messages-due-to-conflicting-packages/176936/3
## devtools::install_version("Matrix",version = "1.6.1.1")
# devtools::install_version('Seurat', version = '4.4.0')

for (pkg in c("argparse", #"Seurat",
              "scater", "scuttle", "ggplot2")){
    if (! pkg  %in% installed.packages()) {
        
        BiocManager::install(pkg, update = FALSE) 

    }
}

