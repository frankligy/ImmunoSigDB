# syntax=docker/dockerfile:1

FROM debian:bullseye
USER root


RUN apt-get update \
    && apt-get install -y wget bzip2 parallel git perl build-essential zip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


# install the main python environment including ncbi orffinder
RUN wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh -O ~/anaconda.sh \
    && /bin/bash ~/anaconda.sh -b -p /opt/conda \
    && rm ~/anaconda.sh

ENV PATH /opt/conda/bin:$PATH

WORKDIR /
RUN conda create -n python_env python=3.10
RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN /opt/conda/envs/python_env/bin/pip install gseapy biopython dash gunicorn

# copy the file
COPY ./static /static
COPY ./assets /assets

# copy the code
COPY web.py /

# modify permission
RUN chmod -R 777 /static
RUN chmod -R 777 /assets
RUN chmod 777 /web.py

# web
CMD /opt/conda/envs/python_env/bin/gunicorn --bind 0.0.0.0:$PORT web:server


