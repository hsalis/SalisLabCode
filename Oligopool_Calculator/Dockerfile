FROM ubuntu:latest

RUN apt-get update \
    && apt-get install -y python3 python3-pip build-essential zip wget pypy graphviz nano apt-transport-https ca-certificates curl software-properties-common intel-mkl git-core cmake pandoc memcached mpich mpich-doc libgfortran5 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install numpy pandas matplotlib seaborn scipy scikit-learn biopython dm-sonnet sonnet graphs mpi4py openpyxl leveldb==0.201 networkx nrpcalc statsmodels msgpack pyfastx edlib dinopy bounter numba primer3-py

 