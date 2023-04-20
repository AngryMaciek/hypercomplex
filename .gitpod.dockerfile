FROM condaforge/mambaforge:latest

RUN \
apt-get update -y \
&& apt-get install libmpfr-dev -y \
&& apt-get clean

# Package conda-build needs to be installed in base env.
# https://github.com/conda/conda-build/issues/3813
# However, installing through YAML configuration of Gitpod
# results in permission denied error; this needs to be executed
# on the conteiner build level.
# ~AngryMaciek
RUN conda install -c conda-forge \
  "conda-build >= 3.23.3" "boa >= 0.14.0" "conda-verify >= 3.1.1"
