FROM debian:jessie

RUN apt-get update && \
	apt-get -y install build-essential gcc g++ gfortran swig wget bzip2 ca-certificates libmysqlclient-dev && \
	rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.continuum.io/archive/Anaconda-2.3.0-Linux-x86_64.sh && \
    /bin/bash Anaconda-2.3.0-Linux-x86_64.sh -b -p /opt/conda && \                                 
    rm Anaconda-2.3.0-Linux-x86_64.sh && \                                                         
    /opt/conda/bin/conda install --yes conda==3.10.1

ENV PATH /opt/conda/bin:$PATH

ADD . .

RUN python setup.py install --prefix /opt/spower

ENV PATH /opt/spower/bin:$PATH
ENV PYTHONPATH /opt/spower/lib/python2.7/site-packages:$PYTHONPATH
