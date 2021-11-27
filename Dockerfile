FROM ubuntu:20.04
WORKDIR /root
RUN apt-get update && apt-get -y install git python3-dev python3-pip libxrender1 vim
RUN pip install Django social-auth-app-django rdkit-pypi pandas annoy
RUN git clone https://github.com/n-yoshikawa/molecule-sanitizer
