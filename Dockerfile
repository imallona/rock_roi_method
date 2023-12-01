## FROM ubuntu:22.04
FROM bioconductor/bioconductor:RELEASE_3_17

LABEL maintainer="Izaskun Mallona izaskun.mallona@gmail.com"
LABEL version="v0.1"

RUN useradd -m rock

RUN apt-get update && \
    apt install -y python3-pip python-is-python3 samtools wget less git nano bedtools

WORKDIR /home/rock
COPY ./data data/
COPY ./ .

RUN pip install snakemake pandas

# STAR, subread and R packages installs

RUN mkdir -p /home/rock/soft/star && \
    cd /home/rock/soft/star && \
    wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz && \
    tar -xzf 2.7.10b.tar.gz && \
    cd STAR-2.7.10b/source && \
    make && \
    mkdir -p /home/rock/soft/subread && \
    cd /home/rock/soft/subread && \
    wget 'https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-source.tar.gz/download' && \
    tar xzvf download && \
    cd subread-2.0.6-source/src && \
    make -f Makefile.Linux && \
    ln -s  /home/rock/soft/subread/subread-2.0.6-source/bin/featureCounts /usr/local/bin/featureCounts && \
    R -e 'BiocManager::install(c("argparse", "SingleCellExperiment", "Matrix", "Seurat", \
          "scater", "scuttle", "ggplot2"), update = FALSE)' && \
    chown -R rock:rock /home/rock && \
    cd /usr/local/bin && \
    wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/faSize && \
    chmod +x faSize && \
    wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bedGraphToBigWig && \
    chmod +x bedGraphToBigWig
    
USER rock

ENV R_LIBS_USER=/home/rock/Rlibs:/usr/local/lib/R/site-library:/usr/local/lib/R/library
ENV R_LIBS=$R_LIBS_USER

RUN export R_LIBS && \
    mkdir -p /home/rock/Rlibs

ENV PATH=/home/rock/soft/star/STAR-2.7.10b/source:$PATH
RUN echo "export PATH=$PATH" >> ~/.bashrc

## showcase the method with simulated data and using snakemake (without conda)

#RUN cd /home/rock && \
#   snakemake -s Snakefile --cores 1 --configfile config.yaml -p
