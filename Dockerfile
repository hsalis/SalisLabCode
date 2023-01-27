FROM ubuntu:18.04
LABEL Description="This Docker image is used to run software developed by the Salis Lab" Version="1.0"
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y ubuntu-server
RUN apt-get update && apt-get install -y apt-utils
RUN apt-get update && apt-get install -y python-pip python-dev build-essential zip wget build-essential pypy graphviz nano \
                                         apt-transport-https ca-certificates curl software-properties-common \
                                         libreadline-gplv2-dev libncursesw5-dev libssl-dev tk-dev libgdbm-dev libc6-dev libbz2-dev libssl-dev libffi-dev pandoc  \
                                         libcr-dev mpich mpich-doc \
                                         && rm -rf /var/lib/apt/lists/*

RUN pip2 install --upgrade pip
RUN pip2 install --upgrade virtualenv
RUN pip2 install mpmath==1.1.0
RUN pip2 install biopython==1.76 boto certifi configparser cycler python-dateutil numpy openpyxl pandas pexpect pickleshare prompt-toolkit==1.0.15 ptyprocess pyExcelerator deap decorator entrypoints enum34 functools32 ipykernel ipython pathos pygments pypandoc pyparsing cssselect==1.1.0 pyquery pysvg python_dateutil pytz pyzmq jinja jsonschema scipy matplotlib mpi4py sympy slacker
RUN pip2 install nbconvert nose xlrd xlwt terminado tornado wcwidth nbformat mistune==0.8.4 gprof2dot pybloom_live==3.1.0 psutil sortedcontainers sqlitedict bayesian-optimization networkx edlib npm requests boto3 python-memcached dnaplotlib Pillow treeinterpreter
RUN pip2 install --upgrade hyperloglog scipy numpy

RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.11.tar.gz && \
    tar xvfz ViennaRNA-2.4.11.tar.gz && \
    cd ViennaRNA-2.4.11 && \
    ./configure PYTHON2_DIR='/usr/local/lib/python2.7/dist-packages' PYTHON2_EXECDIR='/usr/local/lib/python2.7/dist-packages' && \
    make && \
    make install
