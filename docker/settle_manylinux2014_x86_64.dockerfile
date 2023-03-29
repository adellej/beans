FROM quay.io/pypa/manylinux2014_x86_64

LABEL maintainer="martin.cupack@curtin.edu.au"

# this fails as the mirrors are no longer up
# RUN yum -y update
# RUN yum clean all

ARG PY_VER_MAJOR=3
ARG PY_VER_MINOR=8

RUN ln -s /opt/python/cp${PY_VER_MAJOR}${PY_VER_MINOR}-cp${PY_VER_MAJOR}${PY_VER_MINOR}/bin/python3 /usr/local/bin/python3; \
    ln -s /opt/python/cp${PY_VER_MAJOR}${PY_VER_MINOR}-cp${PY_VER_MAJOR}${PY_VER_MINOR}/bin/pip3 /usr/local/bin/pip3
    
RUN which python3; \
    python3 --version

RUN python3 -m pip install --upgrade build twine

# clone settle and build it
ARG BEANS_REPO="https://github.com/ADACS-Australia/beans.git"
RUN cd /usr/src; \
    git clone ${BEANS_REPO}

RUN cd /usr/src/beans/settle; \
    python3 -m build; \
    python3 -m pip install .

RUN cd /usr/src/beans/settle; \
    python3 tests/test_settle_sft.py

# this is interactive. run in na shell.
# 
#RUN python3 -m twine upload --verbose --repository testpypi dist/*
    
    
