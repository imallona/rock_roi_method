FROM ubuntu:22.04

RUN useradd -m rock

RUN apt-get update && \
    apt install -y python3-pip python-is-python3 samtools wget

     
COPY ./main/module.tar.gz /home/rock/

## python deps

RUN cd /home/rock/ && \
    tar xzvf module.tar.gz && \
    cd module && \
    pip install -r rock_n_roi_requirements.txt

# STAR

RUN mkdir -p /home/rock/soft/star && \
    cd $_ && \
    wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz && \
    tar -xzf 2.7.10b.tar.gz && \
    cd STAR-2.7.10b && \
    cd source && \
    make && \
    make install && \
    export PATH=/home/rock/soft/star/STAR-2.7.10b/source:$PATH

USER rock

RUN STAR --version


