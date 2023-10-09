FROM ubuntu:22.04

LABEL maintainer="izaskun.mallona@gmail.com"

RUN useradd -m rock

RUN apt-get update && \
    apt install -y python3-pip python-is-python3 samtools wget less
     
COPY ./main/module.tar.gz /home/rock/
COPY ./data/alien_*bam /home/rock/

## rockandroi python module

RUN cd /home/rock/ && \
    tar xzvf module.tar.gz && \
    cd module && \
    pip install -r rock_n_roi_requirements.txt && \
    mkdir -p data && \
    ln -s /home/rock/alien_wta.bam ./data/wta_alien.bam && \
    ln -s /home/rock/alien_tso.bam ./data/tso_alien.bam

# STAR

RUN mkdir -p /home/rock/soft/star && \
    cd /home/rock/soft/star && \
    wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz && \
    tar -xzf 2.7.10b.tar.gz && \
    cd STAR-2.7.10b/source && \
    make

     
# subread

RUN mkdir -p /home/rock/soft/subread && \
    cd /home/rock/soft/subread && \
    wget 'https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-source.tar.gz/download' && \
    tar xzvf download && \
    cd subread-2.0.6-source/src && \
    make -f Makefile.Linux && \
    ln -s  /home/rock/soft/subread/subread-2.0.6-source/bin/featureCounts /usr/local/bin/featureCounts

USER rock

ENV PATH=/home/rock/soft/star/STAR-2.7.10b/source:$PATH
RUN echo "export PATH=$PATH" >> ~/.bashrc

RUN echo $PATH
# RUN STAR --version


