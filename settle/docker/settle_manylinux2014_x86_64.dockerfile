FROM quay.io/pypa/manylinux2014_x86_64

LABEL maintainer="martin.cupack@curtin.edu.au"

# select python version - works with 3.6 - 3.11
ARG PY_VER_MAJOR=3
ARG PY_VER_MINOR=8

RUN ln -s /opt/python/cp${PY_VER_MAJOR}${PY_VER_MINOR}-cp${PY_VER_MAJOR}${PY_VER_MINOR}*/bin/python3 /usr/local/bin/python3; \
    ln -s /opt/python/cp${PY_VER_MAJOR}${PY_VER_MINOR}-cp${PY_VER_MAJOR}${PY_VER_MINOR}*/bin/pip3 /usr/local/bin/pip3
    
RUN which python3; \
    python3 --version

RUN python3 -m pip install --upgrade build twine

# clone settle and build it
ARG GIT_REPO="https://github.com/ADACS-Australia/beans.git"
ARG BASE_DIR="/usr/src"
ARG SETTLE_DIR=${BASE_DIR}"/beans/settle"

RUN cd ${BASE_DIR}; \
    git clone ${GIT_REPO}

RUN cd ${SETTLE_DIR}; \
    python3 -m build; \
    python3 -m pip install .

RUN cd ${SETTLE_DIR}; \
    python3 tests/test_settle_sft.py

RUN cd ${SETTLE_DIR}; \
    auditwheel repair dist/*.whl

# NOTE ... INSTRUCTIONS ...
# NOTE this last step is interactive, to enter username and apssword.
# NOTE Run in na shell inside the container:
# docker run -it settle_manylinux2014_x86_64:latest bash
# cd ${SETTLE_DIR} # (cd /usr/src/beans/settle/)
# NOTE auditwheel puts the fixed wheel into a wheelhouse.
# NOTE remove "--repository testpypi" to install on real PyPI
# python3 -m twine upload --verbose --repository testpypi wheelhouse/*
    
    
