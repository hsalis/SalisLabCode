FROM ubuntu:22.04

LABEL Description="This Docker image is used to Python 3 code, mainly for machine learning models" Version="1.0"

ENV DEBIAN_FRONTEND noninteractive
ENV HOME /home/ubuntu/
WORKDIR $HOME

RUN apt-get update && apt-get install -y python3 python3-pip build-essential zip wget pypy graphviz nano \
                                         apt-transport-https ca-certificates curl software-properties-common \
                                         intel-mkl git-core cmake \
                                         pandoc memcached \
                                         gcc gfortran g++ openmpi-bin openmpi-common libopenmpi-dev libopenblas-base liblapack3 liblapack-dev make \
                                         && rm -rf /var/lib/apt/lists/*

RUN pip3 install numpy pandas matplotlib seaborn scipy sklearn lightgbm biopython tensorflow tensorflow-probability dm-sonnet sonnet graphs mpi4py openpyxl zipfile

RUN git clone https://github.com/biointec/blamm.git && \
    cd blamm && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install

RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_5_x/ViennaRNA-2.5.0.tar.gz && \
    tar xvfz ViennaRNA-2.5.0.tar.gz && \
    cd ViennaRNA-2.5.0 && \
    ./configure && \
    make && \
    make install

RUN pip3 install jupyterlab
