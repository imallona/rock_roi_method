FROM ubuntu:22.04

RUN useradd -m rock
USER rock

COPY ./main/module $HOME/ 

RUN cd $HOME/module
RUN python --version
RUN pip --version


RUN uname -a

