#
# Good references:
# https://denibertovic.com/posts/handling-permissions-with-docker-volumes/
# https://www.fromlatest.io
#
# ~AngryMaciek

##### BASE IMAGE #####
FROM bitnami/minideb:bullseye
# INFO: https://github.com/bitnami/minideb

##### METADATA #####
LABEL base.image="bitnami/minideb:bullseye"
LABEL version="1.0.0"
LABEL software="hypercomplex"
LABEL software.description="Header-only C++ template library for lattice-based cryptosystems in high-dimensional algebras"
LABEL software.documentation="https://angrymaciek.github.io/hypercomplex"
LABEL software.website="https://github.com/AngryMaciek/hypercomplex"
LABEL software.license="Apache-2.0"
LABEL software.tags="Mathematics"
LABEL maintainer="Maciek Bak"
LABEL maintainer.email="wsciekly.maciek@gmail.com"

##### INSTALL SYSTEM-LEVEL DEPENDENCIES #####
RUN install_packages \
  ca-certificates doxygen g++ gnupg2 git gosu libmpfr-dev make python3 python3-pip valgrind \
  && pip install --no-cache-dir cpplint pre-commit

##### SET ENVIROMENTAL VARIABLES #####
ENV LANG C.UTF-8

##### PREPARE WORKING DIRECTORY #####
VOLUME /hypercomplex
WORKDIR /hypercomplex

##### SETUP ENTRYPOINT #####
COPY entrypoint.sh /bin/entrypoint.sh
RUN /bin/bash -c "chmod +x /bin/entrypoint.sh"
ENTRYPOINT ["/bin/entrypoint.sh"]
CMD ["/bin/bash"]
